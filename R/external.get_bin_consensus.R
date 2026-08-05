#' Compute consensus BIN taxonomy
#'
#' @description
#' Computes and returns consensus taxonomic identifications for each BIN in search results or a BCDM data frame.
#'
#' @details
#' Consensus is defined as any name that exceeds the specified `threshold`, expressed as a proportion of
#' records with a concordant identification (i.e., same name, same rank). The function steps backwards
#' (i.e., from subspecies to kingdom) through the eligible `ranks` to determine the lowest available concordant
#' identification that meets the criteria specified by `threshold`, `min.ids`, and `enforce.scientific`.
#' Different thresholds can be supplied for each rank, if desired (either as a vector of equal length to ranks
#' or as a named list). The function can also be applied to any other grouping variable by modifying `groups`.
#'
#' The provided `bold.search.res` input can be a search result object from \code{\link{bold_parquet_search}}
#' or a BCDM data frame. Alternatively, it can be any data frame or data table minimally containing `bin_uri`
#' (or other grouping variable) and taxonomic identifications for all available records.
#'
#' **Important Note**: This function performs operations on the input data and may be slow when applied to very large datasets or run on systems with limited      #' resources. Before using this function, check the size of the `bold.search.res`object using \code{\link{get_concise_summary}} and proceed with caution.
#'
#' @param bold.search.res A `tbl_sql` object obtained from \code{\link{bold_parquet_search}} or a data frame or data table in BCDM format.
#' @param ranks A character vector of ranks to consider for consensus identifications. Defaults to the standard BOLD ranks.
#' @param threshold Numeric value(s) between 0 and 1 indicating the minimum proportion of records in a BIN that must must share the same taxonomic assignment to establish a consensus. Supply as a single value, a vector with length equal to the number of taxonomic ranks considered, or a named list with names corresponding to specific ranks. If supplied as a named list, an optional "default" value can be set for any ranks that are not explicitly specified (e.g., \code{threshold = list(species = 0.95, default = 0.75)}). The default is 1.0 requiring complete agreement among records at all ranks.
#' @param min.ids Numeric value(s) indicating the minimum number of identifications needed to establish a consensus (names with fewer identifications are still included when calculating proportions). Supply as a single value, a vector with length equal to the number of taxonomic ranks in consideration, or a named list with names corresponding to specific ranks. If supplied as a named list, an optional "default" value can be set for any ranks that are not explicitly specified (e.g., \code{min.ids = list(family = 1, default = 2)}). The default is 2, requiring a minimum of two identifications at any rank.
#' @param enforce.scientific A logical value indicating whether non-scientific, provisional names should be ignored when determining consensus. Default value is TRUE, meaning non-scientific names are ignored.
#' @param groups Grouping variable. Default value is "bin_uri".
#' @param discord.format String indicating the desired output format for the `discordant_ids` column. Can be one of "text" or "list". If "text" (the default), the output is a string column with comma-separated values in the format "Taxon (proportion)". If "list", the output is a list column with names indicating competing identifications and values indicating proportions of discordant identifications for each taxon.
#'
#' @returns A table of consensus identifications for each BIN (or other grouping variable), with the following columns:
#'    `bin_uri`, `member_count`, `concordant_rank`, `concordant_id`, `discordant_rank`, `discordant_ids`.
#'
#' @importFrom data.table as.data.table fcase setnames set copy
#'
#' @examples
#' \dontrun{
#'
#' # Search BOLD data package
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   taxonomy = "Coleoptera",
#'   geography = "Canada"
#' )
#'
#' # Compute strict consensus identifications for BINs in searched data
#' strict_consensus <- get_bin_consensus(
#'   bold.search.res = bold_search,
#'   threshold = 1.0
#' )
#'
#' # Compute identifications concordant among at least 75% of BIN members,
#' # as long as at least three records carry the majority identifications
#' # in each BIN.
#' bin_ids <- get_bin_consensus(
#'   bold.search.res = bold_search,
#'   threshold = 0.75,
#'   min.ids = 3
#' )
#'
#' # Include non-scientific names (i.e., interim taxonomy or placeholder names)
#' # in consideration of BIN consensus.
#' bin_consensus <- get_bin_consensus(
#'   bold.search.res = bold_search,
#'   threshold = 0.9,
#'   enforce.scientific = FALSE
#' )
#' }
#'
#' @export
get_bin_consensus <- function(
  bold.search.res,
  ranks = c("kingdom", "phylum", "class", "order", "family", "subfamily", "tribe", "genus", "species", "subspecies"),
  threshold = 1.0,
  min.ids = 1,
  enforce.scientific = TRUE,
  groups = "bin_uri",
  discord.format = c("text", "list")
) {
  # Check input format
  is_tbl_sql <- isTRUE(try(check.tbl.sql(bold.search.res), silent = TRUE))
  if (!is_tbl_sql && !is.data.frame(bold.search.res)) stop("`bold.search.res` must be either a bold_parquet_search output (tbl_sql / dbplyr table) or a data frame / data table.")
  # Check parameters
  stopifnot(
    "One or more provided `ranks` is/are missing from `bold.df`." = all(ranks %in% colnames(bold.search.res)),
    "Provided `groups` column is missing from `bold.df`." = (groups %in% colnames(bold.search.res)),
    "`threshold` value(s) must be one or more real numbers (i.e., doubles) between 0 and 1." = is.double(unlist(threshold)) & all(unlist(threshold) >= 0) & all(unlist(threshold) <= 1),
    "`threshold` must be either a single number, a vector of unnamed numbers equal in length to `ranks`, or a named list or vector of numbers with names corresponding to ranks." = ((length(threshold) == 1) | (length(threshold) == length(ranks)) | (!is.null(names(threshold)))),
    "`min.ids` value(s) must be one or more whole numbers greater than zero." = is.numeric(unlist(min.ids)) & all(unlist(min.ids) > 0) & all(unlist(min.ids) %% 1 == 0),
    "`min.ids` must be either a single number, a vector of unnamed numbers equal in length to `ranks`, or a named list or vector of numbers with names corresponding to ranks." = ((length(min.ids) == 1) | (length(min.ids) == length(ranks)) | (!is.null(names(min.ids)))),
    '`discord.format` must be one of "list" or "text".' = all(discord.format %in% c("list", "text"))
  )
  # Collect into data.frame if needed, then create a data.table copy of the data
  dt <- if (is_tbl_sql) {
    bold.search.res %>%
      dplyr::filter((!is.na(bin_uri)) & (bin_uri != "")) %>%
      dplyr::select(all_of(c(groups, ranks))) %>%
      collect() %>%
      as.data.table() %>%
      copy()
  } else {
    as.data.table(copy(bold.search.res))
  }
  # Parse threshold & min.ids parameters and align them with ranks
  parse_param_vector <- function(param) {
    if ((length(param) != 1) | !is.null(names(param))) {
      if (is.null(names(param))) {
        param <- unlist(unname(param))
      } else {
        named <- as.list(param[(names(param) %in% ranks) & (!duplicated(param))])
        default <- ifelse("default" %in% names(param), param[["default"]], max(unlist(param)))
        if ((!length(named) %in% c(0, length(ranks))) & (!"default" %in% names(param))) {
          warning(paste0("Only some ranks found among `", substitute(param), "` values, with no default given; highest value applied to all unspecified ranks."))
        }
        param <- rep(default, length(ranks))
        for (r in names(named)) param[match(r, ranks)] <- named[[r]]
      }
    } else {
      param <- rep(unlist(param), length(ranks))
    }
    return(param)
  }
  threshold <- parse_param_vector(threshold)
  min.ids <- parse_param_vector(min.ids)
  # Replace NA in taxonomy columns with empty values (if ignoring non-scientific names, replace those too)
  dt[, (ranks) := lapply(.SD, function(x) {
    fcase(enforce.scientific & grepl(re_int, x, perl = TRUE), "",
      is.na(x), "",
      grepl("^\\s$", x), "",
      default = as.character(x)
    )
  }), .SDcols = ranks]
  # Convert data table to matrix for faster row access
  mat <- as.matrix(dt[, c(groups, ranks), with = FALSE])
  # Replace trailing "" with NA so that they are not counted as alternative names
  for (i in seq_len(nrow(mat))) {
    row_vals <- mat[i, ]
    non_blank_idx <- which(row_vals != "")
    if (length(non_blank_idx) > 0) {
      last <- max(non_blank_idx)
      if (last < ncol(mat)) {
        mat[i, (last + 1):ncol(mat)] <- NA_character_
      }
    } else {
      mat[i, ] <- NA_character_ # Entire row is blank
    }
  }
  # Convert back to data.table and restore column names
  dt <- as.data.table(mat)
  setnames(dt, c(groups, ranks))
  # Core consensus logic
  get_consistent_taxon <- function(sub_dt, ranks, threshold, min.ids) {
    id_hier <- sapply(ranks, function(x) NULL)
    concordant <- FALSE
    rank_set <- ranks
    # Ensure min.ids does not exceed group size
    if (any(min.ids > nrow(sub_dt))) {
      for (i in seq_along(min.ids)) min.ids[[i]] <- nrow(sub_dt)
    }
    # Expand threshold and min.ids parameters into full vectors if applicable
    if (length(threshold) == 1) {
      threshold <- rep(threshold, length(ranks))
    }
    if (length(min.ids) == 1) {
      min.ids <- rep(min.ids, length(ranks))
    }
    result <- list(
      member_count = nrow(sub_dt),
      concordant_rank = NA_character_,
      concordant_id = NA_character_,
      concordant_id_count = 0L,
      discordant_rank = NA_character_,
      discordant_ids = list(),
      discordant_id_count = 0L
    )
    for (rank_col in rev(ranks)) { # Step backwards through ranks
      rank_threshold <- threshold[which(ranks == rank_col)]
      rank_min_ids <- min.ids[which(ranks == rank_col)]
      vals <- sub_dt[[rank_col]]
      filtered <- table(vals[!is.na(vals)])
      name_vals <- proportions(filtered)
      props <- proportions(filtered)[(proportions(filtered) >= rank_threshold) & (filtered >= rank_min_ids)]
      names(name_vals) <- sub("^$", "<None>", names(name_vals))
      if ((length(props) != 1) & (length(unique(filtered)) > 0)) {
        concordant <- FALSE
        if (id_hier[rank_col] != "") {
          rank_set <- ranks[0:(which(ranks == rank_col) - 1)]
          id_hier <- id_hier[rank_set]
        }
        if (length(name_vals) > 1) {
          result$discordant_rank <- rank_col
          result$discordant_ids <- list(stats::setNames(as.vector(name_vals), names(name_vals)))
          result$discordant_id_count <- sum(filtered)
        }
      } else if ((length(props) == 1) && (names(props)[1] != "")) {
        if (!concordant) {
          result$concordant_rank <- rank_col
          result$concordant_id <- names(props)[1]
          result$concordant_id_count <- unname(filtered[names(props)[1]])
        }
        concordant <- TRUE
        if (is.null(id_hier[[rank_col]]) || is.na(id_hier[[rank_col]]) || is.na(names(props)[1]) || (names(props)[1] != id_hier[[rank_col]])) {
          rank_set <- ranks[0:which(ranks == rank_col)]
          id_hier <- as.list(sub_dt[get(rank_col) == names(props)[1], .SD, .SDcols = rank_set][1])
        }
      }
    }
    for (r in setdiff(ranks, names(id_hier))) {
      id_hier[r] <- NA_character_
    }
    result[ranks] <- id_hier
    return(result)
  }
  # Generate summary of consensus by BIN
  consensus <- dt[!is.na(get(groups)), do.call(get_consistent_taxon, list(.SD, ranks, threshold, min.ids)), by = groups, .SDcols = ranks]
  # Convert discordant_ids to text if appropriate
  if (discord.format[1] == "text") {
    data.table::set(consensus,
      j = "discordant_ids",
      value = sapply(consensus[["discordant_ids"]], function(x) {
        if (length(x) == 0) {
          return("")
        }
        sort(x, decreasing = TRUE)
        pairs <- paste0(names(x), " (", formatC(as.numeric(x), format = "f", digits = 2), ")")
        paste(pairs, collapse = ", ")
      })
    )
  }
  return(consensus)
}
