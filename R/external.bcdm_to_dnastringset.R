#' Convert the BOLD parquet search into a DNAStringSet object
#'
#' @description Converts the sequence data from the search result into a DNAStringSet object for downstream multiple sequence alignment with customized headers.
#'
#' @details This function transforms the search results from `bold_parquet_search` into a `DNAStringSet` object suitable for downstream data analyses. The `cols_for_seq_names` argument lets users create custom headers using the BCDM column names.
#'
#' @param bold.search.res A `tbl_sql` object containing BOLD search results.
#' @param marker Character vector specifying the genetic marker.
#' @param cols_for_seq_names Character vector of field names to include in the header.
#'
#' @return A `DNAStringSet` object.
#'
#' @importFrom dplyr select filter mutate compute collect all_of %>%
#' @importFrom dbplyr sql remote_con
#' @importFrom stats setNames
#' @importFrom DBI dbQuoteIdentifier
#' @examples
#' \dontrun{
#'
#' # Search the BOLD data package
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   taxonomy = "Coleoptera",
#'   geography = "Canada",
#'   basecount = c(500, 660)
#' )
#'
#' # Get the DNAStringset object (library Biostrings needs to be imported beforehand)
#'
#' bold.dnastringset <- bcdm_to_dnastringset(bold_search,
#'   marker = "COI-5P",
#'   cols_for_seq_names = c("processid", "family")
#' )
#' }
#'
#' @export
bcdm_to_dnastringset <- function(bold.search.res,
                                 marker = NULL,
                                 cols_for_seq_names) {
  # Check if input is a tbl_sql (helper function you already have)
  check.tbl.sql(bold.search.res)
  # Check for Biostrings package
  check_biostrings <- function() {
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
      stop(
        "Error: The 'Biostrings' package is required for this function. ",
        "Please install it with: BiocManager::install('Biostrings')"
      )
    }
  }
  # Ensure required columns exist
  required_cols <- c("nuc", "marker_code", cols_for_seq_names)
  missing_cols <- setdiff(required_cols, colnames(bold.search.res))
  # To check if minimum columns needed for creating a DNAstringset are available
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from the input: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  # Start processing tbl_sql
  seq.data <- bold.search.res %>%
    select(
      .data$nuc,
      .data$marker_code,
      all_of(cols_for_seq_names)
    ) %>%
    filter(!is.na(.data$nuc)) %>%
    filter(.data$nuc != "") %>%
    mutate(nuc = sql("REGEXP_REPLACE(nuc, '-', '')"))
  # Apply marker filter separately
  if (!is.null(marker)) {
    seq.data <- seq.data %>%
      filter(.data$marker_code %in% marker)
  }
  con <- dbplyr::remote_con(seq.data)
  # Quote column names
  quoted_cols <- DBI::dbQuoteIdentifier(con, cols_for_seq_names)
  # Build sequence names from specified columns
  obtain.seq.from.data <- seq.data %>%
    select(nuc, all_of(cols_for_seq_names)) %>%
    mutate(
      # Quote each column name to avoid SQL keyword issues
      msa.seq.name = sql(
        paste0(
          paste0('"', cols_for_seq_names, '"', collapse = " || '|' || ")
        )
      )
    ) %>%
    select(msa.seq.name, nuc)
  # Pull into R and convert to named character vector
  seq.from.data <- obtain.seq.from.data %>%
    collect() %>%
    {
      setNames(as.character(.$nuc), .$msa.seq.name)
    }
  # Convert to DNAStringSet
  dna.4.align <- DNAStringSet(seq.from.data)
  return(dna.4.align)
}
