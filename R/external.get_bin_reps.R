#' Select representative records by BIN
#'
#' @description
#' Obtain one or more representative record(s) from each BIN in search results or a BCDM data frame.
#' BIN representatives can also be selected for each unique BIN-taxon combination. The chosen
#' selection criteria are applied in sequence to select representatives, thus criteria should be given
#' in order of priority. If multiple records are tied after all criteria have been applied, the tie is
#' broken at random by default. Alternatively, it is possible to obtain reproducible results for the
#' same input values by setting a random seed.
#'
#' @details
#' # Input
#' The provided `bold.search.res` input can be a search result object from \code{\link{bold_parquet_search}}
#' or a BCDM data frame. Alternatively, it can be any data frame or data table minimally containing
#' `bin_uri`, unique record identifiers (e.g. `processid` or `sampleid`), taxonomic identifications
#' for all available records, and any fields relevant to the provided selection `criteria` (see below).
#'
#' **Important Note**: This function is not optimized for very large `tbl_sql` search results,
#' particularly those with many unique BINs. On the other hand, input provided as a data frame or data table
#' can be processed much more efficiently. Therefore, if you intend to keep the full search results in addition
#' to BIN representatives, it is recommended that you first collect the results into a data frame using
#' \code{\link{bold_search_collect}}.
#'
#' # Criteria
#' Selection `criteria` must be listed in order of priority. Each one acts upon data from a particular BCDM
#' field, which must be present in the input data, as indicated below. Available criteria include the following:
#' \describe{
#'   \item{**vouchered**}{If `TRUE`, prioritize records with known voucher repositories over those mined from databases like GenBank. Setting this to `FALSE` will prioritize records *without* vouchers. To exclude this criterion, simply omit it. (Required field: `inst`.)}
#'   \item{**seq_length**}{Can be used either to specify target barcode sequence length as an integer, to preferentially select longer or shorter sequences ("longest", "shortest"), or to select sequences that match the modal* barcode length for each BIN ("COI_auto"). If a target length is provided, sequences closest to that target are prioritized. *Modal barcode length ("COI_auto" option) is determined after first rounding sequence lengths to the nearest full codon; in the event of multiple modes, the one closest to 658bp is chosen. (Required field: `nuc_basecount`.)}
#'   \item{**id_method**}{A character vector listing preferred identification methods in order of priority. The function definition lists all available values per the BCDM specification. (Required field: `identification_method`.)}
#'   \item{**inst**}{A character vector listing preferred voucher specimen repositories in order of priority. Values not present in the input data are ignored. (Required field: `inst`.)}
#'   \item{**coll_date**}{Prioritize recently collected specimens ("latest") or those collected longest ago ("oldest"). Records without dates are selected last in either case. (Required field: `collection_date_start`.)}
#'   \item{**seq_date**}{Prioritize recently uploaded sequences ("latest") or those with the earliest upload date ("oldest"). (Required field: `sequence_upload_date`.)}
#'  }
#'
#' Default selection criteria are as follows:
#' ```R
#' criteria = list(vouchered = TRUE,
#'                 seq_length = "COI_auto",
#'                 id_method = c("Morphology", "Morphology and sequence based",
#'                               "Image based", "Image and sequence based",
#'                               "Tree based", "BIN based", "BOLD ID Engine",
#'                               "Other sequence based approach", "Other"),
#'                 inst = "Centre for Biodiversity Genomics",
#'                 coll_date = "latest",
#'                 seq_date = "latest")
#' ```
#'
#' # Sampling representatives by taxon
#' When \code{by.taxon = TRUE}, representatives are selected for each unique combination of `bin_uri` and
#' `identification`. For example: Sampling single representatives from a BIN containing records identified as
#' "Agrilinae", "Agrilus", and "Agrilus VVG_sp.42" will yield three records—one for each name. Two
#' additional parameters can be used to tune this behaviour: `enforce.scientific`
#' and `non.redundant.taxa`.
#'
#' When \code{enforce.scientific = TRUE}, interim / provisional names are ignored when considering unique
#' identifications. For example: "Agrilus VVG_sp.42" and "Agrilus" will both be treated as "Agrilus". For the
#' BIN from the previous example, this will yield two representatives: one each for "Agrilus" and "Agrilinae".
#' Note that records with such names may still be selected as representatives if they are prioritized
#' according to the provided criteria, or if they are the only available representatives in a BIN.
#'
#' When \code{non.redundant.taxa = TRUE}, the identifications in a BIN are compared in terms of their
#' full taxonomic classification to determine the most specific name available for each distinct taxonomic
#' lineage, and any records thus identified will be selected first. For example: "Agrilinae", "Agrilus",
#' and "Agrilus VVG_sp.42" will all be treated as a single lineage when selecting representatives by taxon.
#' This taxonomic de-duplication behaviour is also affected by the previous parameter: when `enforce.scientific`
#' and `non.redundant.taxa` are both `TRUE`, the names "Agrilinae", "Agrilus", "Agrilus VVG_sp.42", and
#' "Agrilus crataegi" all collapse to the single taxon "Agrilus crataegi". When `non.redundant.taxa` is `TRUE`,
#' and `enforce.scientific` is `FALSE`, the same four names collapse to "Agrilus VVG_sp.42" and "Agrilus crataegi".
#' Note that this option prioritizes records with the lowest available taxonomic information, and in some cases
#' can thus override the selection criteria, e.g., if otherwise preferred records have unresolved identifications.
#' Note also that the function will always return at least `Nreps` records per BIN where possible; if there are
#' fewer than `Nreps` records bearing the lowest non-redundant identification(s) in a BIN, additional records
#' will be selected as back-fill beginning with the next-most specific identification, working upwards.
#'
#' # Reproducibility
#' It is possible to consistently obtain the same representatives using the same input parameters by
#' supplying a random `seed`. However, the same data provided either as a `tbl_sql` object or a
#' data frame may yield different representative records for a given random seed due to differences
#' in sorting and tie-breaking behaviours across backends (i.e. DuckDB vs. data.table).
#'
#' @param bold.search.res A `tbl_sql` object obtained from \code{\link{bold_parquet_search}} or a data frame or data table in BCDM format.
#' @param Nreps Integer indicating the maximum number of representatives to select for each BIN (or BIN-taxon combination).
#' @param by.taxon Logical value indicating whether to select representatives for each unique combination of BIN and taxonomic identification. If `TRUE`, the additional parameters `non_redundant_taxa` and `enforce_scientific` are also applied. (See 'Sampling representatives by taxon' for more details.)
#' @param enforce.scientific Logical value indicating whether to ignore non-scientific, provisional names for the purposes of sampling representatives by taxon. Ignored if `by.taxon` is `FALSE`. (See 'Sampling representatives by taxon' for more details.)
#' @param non.redundant.taxa Logical value, when sampling by taxon (by.taxon = TRUE), representatives are selected from the lowest available taxonomic rank within each lineage. For example, the identifications "Apidae", "Bombus", and "Bombus terrestris" are considered to belong to a single taxonomic lineage, and records assigned to "Bombus terrestris" will be selected first. Ignored if `by.taxon` is `FALSE`. (See 'Sampling representatives by taxon' for more details.)
#' @param criteria Named list of selection criteria to apply when sampling representatives, given in priority order. See 'Criteria' for more information and default values.
#' @param seed Optional positive integer to use as a random seed for reproducible tie-breaking. If `NULL` (the default), ties are broken randomly and selected records may differ between runs.
#'
#' @returns A data frame of selected representatives.
#'
#' @usage
#' get_bin_reps(
#'   bold.search.res,
#'   Nreps = 1,
#'   by.taxon = FALSE,
#'   enforce.scientific = FALSE,
#'   non.redundant.taxa = FALSE,
#'   criteria = list(vouchered = TRUE,
#'                   seq_length = c("COI_auto", "longest", "shortest", 658),
#'                   id_method = c("Morphology", "Morphology and sequence based",
#'                                 "Image based", "Image and sequence based",
#'                                 "Tree based", "BIN based", "BOLD ID Engine",
#'                                 "Other sequence based approach", "Other"),
#'                  inst = "Centre for Biodiversity Genomics",
#'                  coll_date = c("latest", "oldest"),
#'                  seq_date = c("latest", "oldest")),
#'   seed = NULL
#' )
#'
#' @examples
#' \dontrun{
#'
#' # Search BOLD data package
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   taxonomy = "Araneae",
#'   geography = "Canada"
#' )
#'
#' # Select three representatives per BIN from the searched data,
#' # prioritizing those with a morphological ID and with vouchers
#' # deposited at the CBG (e.g. in case vouchers need to be examined)
#' bin_reps <- get_bin_reps(
#'   bold.search.res = bold_search,
#'   Nreps = 3,
#'   criteria = list(
#'     inst = "Centre for Biodiversity Genomics",
#'     id_method = c(
#'       "Morphology",
#'       "Morphology and sequence based"
#'     )
#'   )
#' )
#'
#' # Select one representative for each combination of BIN and taxonomic
#' # lineage (scientific names only), with preference for 658-bp barcodes
#' # (e.g. for building a sequence tree of all known taxa)
#' bin_tax_reps <- get_bin_reps(
#'   bold.search.res = bold_search,
#'   Nreps = 1,
#'   by.taxon = TRUE,
#'   enforce.scientific = TRUE,
#'   non.redundant.taxa = TRUE,
#'   criteria = list(seq_length = 658)
#' )
#' }
#'
#' @export
get_bin_reps <- function(
  bold.search.res,
  Nreps = 1,
  by.taxon = FALSE,
  enforce.scientific = FALSE,
  non.redundant.taxa = FALSE,
  criteria = list(
    vouchered = TRUE,
    seq_length = c("COI_auto", "longest", "shortest", 658),
    id_method = c(
      "Morphology", "Morphology and sequence based",
      "Image based", "Image and sequence based",
      "Tree based", "BIN based", "BOLD ID Engine",
      "Other sequence based approach", "Other"
    ),
    inst = "Centre for Biodiversity Genomics",
    coll_date = c("latest", "oldest"),
    seq_date = c("latest", "oldest")
  ),
  seed = NULL
) {
  # Check input format
  is_tbl_sql <- isTRUE(try(check.tbl.sql(bold.search.res), silent = TRUE))
  if (!is_tbl_sql && !is.data.frame(bold.search.res)) stop("`bold.search.res` must be either a bold_parquet_search output (tbl_sql / dbplyr table) or a data frame / data table.")
  # Check parameters
  stopifnot(
    "'Nreps' must be a single numeric value." = is.numeric(Nreps) && length(Nreps) == 1,
    "'by.taxon' must be a single logical value." = is.logical(by.taxon) && length(by.taxon) == 1,
    "'non.redundant.taxa' must be a single logical value." = is.logical(non.redundant.taxa) && length(non.redundant.taxa) == 1,
    "'enforce.scientific' must be a single logical value." = is.logical(enforce.scientific) && length(enforce.scientific) == 1
  )
  if (!missing(criteria)) {
    if (length(criteria$seq_length) > 1) {
      stop("'seq_length' criterion must be a single value.")
    } else if (length(criteria$seq_length) == 1) {
      if (!any(
        (criteria$seq_length %in% c("COI_auto", "longest", "shortest")),
        is.numeric(criteria$seq_length)
      )) {
        stop("'seq_length' criterion must be one of 'COI_auto', 'longest', 'shortest', or a numeric value.")
      }
    }
    stopifnot(
      "'vouchered' criterion must be a single logical value." =
        length(criteria$vouchered) == 0 || is.logical(criteria$vouchered) && length(criteria$vouchered) == 1
    )
    if (length(criteria$coll_date) > 1) stop("'coll_date' criterion must be one of 'latest' or 'oldest'.")
    if (length(criteria$seq_date) > 1) stop("'seq_date' criterion must be one of 'latest' or 'oldest'.")
  } else {
    # Supply default criteria if param is omitted
    criteria <- list(
      vouchered = TRUE,
      seq_length = "COI_auto",
      id_method = c(
        "Morphology", "Morphology and sequence based",
        "Image based", "Image and sequence based",
        "Tree based", "BIN based", "BOLD ID Engine",
        "Other sequence based approach", "Other"
      ),
      inst = "Centre for Biodiversity Genomics",
      coll_date = "latest",
      seq_date = "latest"
    )
  }
  # Default BOLD ranks
  ranks <- c("kingdom", "phylum", "class", "order", "family", "subfamily", "tribe", "genus", "species", "subspecies")
  # Define grouping columns
  select_by <- if (by.taxon) {
    c("bin_uri", "identification")
  } else {
    "bin_uri"
  }
  # Ensure Nreps is a whole number
  Nreps <- round(Nreps, 0)
  # Deduce specimen identifier
  id_col <- intersect(
    c("processid", "sampleid", "fieldid", "museumid", "record_id", "specimenid"),
    colnames(bold.search.res)
  )[[1]]
  # Proceed according to input object type
  bin_reps <- if (is_tbl_sql) {
    # Shuffle data to randomize representative order, or set pre-determined order using random seed
    df_shuffle <- if (!is.null(seed)) {
      bold.search.res %>% dplyr::mutate(.rand = md5(paste0(!!sym(id_col), !!as.character(seed))))
    } else {
      bold.search.res %>% dplyr::mutate(.rand = random())
    }
    df_shuffle <- df_shuffle %>%
      dplyr::filter(!is.na(bin_uri)) %>%
      dplyr::filter(marker_code == "COI-5P") %>%
      dplyr::distinct(!!sym(id_col), .keep_all = TRUE)
    # Build sequence of sort keys to use in query
    sort_sequence <- get_sql_sort(df_shuffle, criteria)
    # Assign temporary column names for sorting
    temp_names <- paste0(".sort_", names(criteria))
    # Assemble expressions for temporary sort columns
    sort_exprs <- setNames(
      lapply(sort_sequence, function(step) step$key),
      temp_names
    )
    # Assemble sequential sort expressions
    arrange_exprs <- c(
      lapply(seq_along(temp_names), function(i) {
        nm <- sym(temp_names[[i]])
        if (isTRUE(sort_sequence[[i]]$desc)) dplyr::expr(desc(!!nm)) else nm
      }),
      list(sym(".rand"))
    )
    # Apply any preliminary transformations (e.g. get modal sequence length per BIN)
    df_prepped <- Reduce(function(df, step) {
      if (!is.null(step$prep)) step$prep(df) else df
    }, sort_sequence, init = df_shuffle)
    # Apply query expressions and collect representatives
    df_reps <- df_prepped %>%
      dplyr::mutate(!!!sort_exprs) %>%
      dplyr::group_by(across(all_of(select_by))) %>%
      dbplyr::window_order(!!!arrange_exprs) %>%
      dplyr::mutate(.row_rank = dplyr::row_number()) %>%
      dplyr::filter(.row_rank <= Nreps) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(!!!arrange_exprs) %>%
      dplyr::collect() %>%
      dplyr::select(-all_of(temp_names), -.row_rank, -.rand, -any_of(c(".bin_mode")))
    # For simplicity when by.taxon == TRUE all unique combinations of `bin_uri` and `identification`
    # are collected initially, then filtered further if applicable
    if (by.taxon) {
      # Move interim names to bottom (may still be selected if they are the only reps that fulfill the criteria)
      if (enforce.scientific) df_reps <- df_reps %>% dplyr::arrange(grepl(re_int, identification, perl = TRUE))
      # Build a taxonomy lookup table with matching indices, which can be safely modified
      id_lookup <- df_reps %>%
        dplyr::select(all_of(c("bin_uri", "identification", "identification_rank", ranks))) %>%
        dplyr::mutate(across(all_of(ranks), ~ case_when(
          enforce.scientific & grepl(re_int, .x, perl = TRUE) ~ NA_character_,
          .x == "" ~ NA_character_,
          grepl("^\\s$", .x) ~ NA_character_,
          .default = as.character(.x)
        )))
      if (enforce.scientific) {
        # Interim names were already suppressed in the previous step; now we re-compute identification & rank
        vals <- as.matrix(id_lookup[ranks])
        mask <- !is.na(vals)
        idx <- max.col(mask, ties.method = "last")
        all_na <- rowSums(mask) == 0L
        id_lookup <- id_lookup %>%
          dplyr::mutate(
            identification = ifelse(all_na, NA_character_, vals[cbind(dplyr::row_number(), idx)]),
            identification_rank = ifelse(all_na, NA_character_, ranks[idx])
          )
        # Group by cleaned identifications and select first record(s) (scientific, if available, due to earlier sort step)
        keep_rows <- id_lookup %>%
          dplyr::mutate(.row = dplyr::row_number()) %>%
          dplyr::group_by(bin_uri, identification) %>%
          dplyr::slice_head(n = Nreps) %>%
          dplyr::ungroup() %>%
          dplyr::pull(.row)
        # Discard unwanted records in both the reps and lookup tables
        id_lookup <- id_lookup %>% dplyr::slice(keep_rows)
        df_reps <- df_reps %>% dplyr::slice(keep_rows)
      }
      if (non.redundant.taxa) {
        id_lookup_by_bin <- split(id_lookup, id_lookup$bin_uri)
        # Group-wise function to find lowest non-redundant identification for each BIN
        filter_ids <- function(df, ...) {
          key <- list(...)[[1]]
          group_lookup <- id_lookup_by_bin[[as.character(key$bin_uri)]]
          group_lookup$.row <- seq_len(nrow(group_lookup))
          ranks_to_check <- ranks[sapply(ranks, function(rank) !all(is.na(group_lookup[[rank]]) | group_lookup[[rank]] == ""))]
          id_counts <- group_lookup %>%
            dplyr::distinct(identification, identification_rank, across(all_of(ranks_to_check)), .keep_all = TRUE)
          id_counts$count <- mapply(
            function(rank, val) sum(id_counts[[rank]] == val, na.rm = TRUE),
            id_counts[["identification_rank"]], id_counts[["identification"]]
          )
          id_counts <- id_counts %>% dplyr::arrange("count")
          unique_ids <- id_counts %>%
            dplyr::filter(count == 1) %>%
            dplyr::pull("identification") %>%
            as.character()
          next_best <- id_counts %>%
            dplyr::filter(count != 1) %>%
            dplyr::pull("identification") %>%
            as.character()
          # Select Nreps for each unique taxon
          keep <- group_lookup %>%
            dplyr::filter(identification %in% unique_ids) %>%
            dplyr::group_by(identification) %>%
            dplyr::slice_head(n = Nreps) %>%
            dplyr::ungroup() %>%
            dplyr::pull(".row")
          # If less than `Nreps` have lowest non-redundant ID, fill out the rest with the next-best IDs
          if (length(keep) < Nreps) {
            keep <- c(
              keep,
              group_lookup %>%
                dplyr::filter(identification %in% next_best) %>%
                dplyr::mutate(identification = factor(identification, levels = next_best)) %>%
                dplyr::arrange(identification) %>%
                dplyr::slice_head(n = Nreps - length(keep)) %>%
                dplyr::pull(".row")
            )
          }
          df %>% dplyr::slice(keep)
        }
        # Apply group-wise filtering while preserving original column order
        col_order <- names(df_reps)
        df_reps <- df_reps %>%
          dplyr::group_by(bin_uri) %>%
          dplyr::group_modify(filter_ids) %>%
          dplyr::ungroup() %>%
          dplyr::relocate(all_of(col_order))
      }
    }
    df_reps
  } else { # If local data is supplied, select BIN reps from there
    # Randomize representative order (using random seed for pre-determined order if provided)
    if (!is.null(seed)) set.seed(seed)
    data <- as.data.table(bold.search.res)[sample(nrow(bold.search.res))]
    data <- unique(data[!is.na(bin_uri) & (bin_uri != "") & (marker_code == "COI-5P")], by = id_col)
    # Filter by identification if applicable
    if (by.taxon) {
      # Move interim names to the bottom (they may still be selected if they match other criteria)
      if (enforce.scientific) data <- data[order(bin_uri, grepl(re_int, identification, perl = TRUE))]
      if (enforce.scientific || non.redundant.taxa) {
        # Build a taxonomy lookup table with matching indices to be mutated (suppressing interim names if applicable)
        id_lookup <- data[, .SD, .SDcols = c("bin_uri", "identification", "identification_rank", ranks)]
        id_lookup[, c(ranks, if (enforce.scientific) c("identification", "identification_rank")) := {
          cleaned <- lapply(.SD, function(x) {
            fcase(enforce.scientific & grepl(re_int, x, perl = TRUE), NA_character_,
              x == "", NA_character_,
              grepl("^\\s$", x), NA_character_,
              default = as.character(x)
            )
          })
          if (enforce.scientific) {
            vals <- do.call(cbind, cleaned)
            mask <- !is.na(vals)
            idx <- max.col(mask, ties.method = "last")
            all_na <- rowSums(mask) == 0L
            c(cleaned, list(
              ifelse(all_na, NA_character_, vals[cbind(seq_len(.N), idx)]),
              ifelse(all_na, NA_character_, ranks[idx])
            ))
          } else {
            cleaned
          }
        }, .SDcols = ranks]
      }
      if (enforce.scientific) {
        # Substitute cleaned identification for later rep selection

        data[, "id_clean" := id_lookup$identification]
        select_by <- c("bin_uri", "id_clean")
      }
    }
    # Build sequence of sort keys to apply to data table
    sort_sequence <- get_dt_sort(data, criteria)
    # Apply sort keys in sequence to select and return representatives
    rep_idx <- data[do.call("order", sort_sequence), .I[seq_len(min(Nreps, .N))], by = select_by]$V1
    data <- data[rep_idx, ]
    if (by.taxon) {
      if (enforce.scientific || non.redundant.taxa) id_lookup <- id_lookup[rep_idx, ]
      if (non.redundant.taxa) {
        id_lookup_by_bin <- split(id_lookup, id_lookup$bin_uri)
        # Find and apply the indices of records with the lowest non-redundant identification in each BIN
        data <- data[id_lookup[, .I[({
          bin_taxa <- id_lookup_by_bin[[.BY$bin_uri]]
          cols_to_use <- c("identification", "identification_rank", ranks[sapply(ranks, function(rank) !all(is.na(bin_taxa[[rank]]) | bin_taxa[[rank]] == ""))])
          bin_taxa <- bin_taxa[, .SD, .SDcols = cols_to_use]
          ids_ranked <- unique(bin_taxa)
          ids_ranked[, id_count := mapply(
            function(rank, id) .SD[.SD[[rank]] == id, .N],
            as.character(identification_rank),
            as.character(identification)
          )]
          data.table::setorder(ids_ranked, id_count)
          keep <- bin_taxa[identification %in% ids_ranked[id_count == 1, identification], which = TRUE]
          # If less than `Nreps` have lowest non-redundant ID, fill out the rest with the next-best IDs
          if (length(keep) < Nreps) {
            next_best <- unique(bin_taxa[identification %in% ids_ranked[id_count != 1, identification], identification])
            add <- bin_taxa[, .(identification, row = .I)][identification %in% next_best][
              , identification := factor(identification, levels = next_best)
            ][
              order(identification)
            ][
              , head(row, Nreps - length(keep))
            ]
            keep <- c(keep, add)
          }
          bin_taxa[, .I] %in% keep
        })], by = "bin_uri"]$V1]
      }
      # Remove cleaned ID column if added previously
      if (enforce.scientific) data[, "id_clean" := NULL]
    }
    data
  }
  return(as.data.frame(bin_reps))
}
