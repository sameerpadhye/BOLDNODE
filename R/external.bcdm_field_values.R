#' Extract unique values from BCDM fields in the BOLD data package
#'
#' @description Extracts distinct values of specified field(s) from the BOLD parquet data.
#'
#' @details This function extracts unique values from one or more specified columns from the BOLD parquet data. It handles both the parquet file and `tbl_sql` objects from `bold_parquet_search` as input. The results can be saved to disk as an `.rds` file for later use.
#'
#' @param input.data Path to the input parquet file or the `bold_parquet_search` result.
#' @param specific.cols Name of the column to extract unique values from.
#' @param save.data Logical value indicating whether to save the results to disk as a .rds file (default: FALSE).
#' @param output.file Path (without extension) for saving results as .rds file (required if save.data = TRUE).
#'
#' @return A list containing unique values from the specified column. If `save.data` = T, a `.rds` file is exported locally.
#'
#' @importFrom dplyr filter distinct collect %>%
#' @importFrom rlang .data
#'
#' @examples
#' \donttest{
#'
#' # Import the parquet file (This is a test parquet file composed of
#' # records of Cerambycidae beetles from Canada)
#' parquet_file <- system.file(
#' "extdata",
#' "test_data.parquet",
#' package = "BOLDNODE"
#' )
#'
#' # Search the BOLD data package
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   marker = "COI-5P"
#' )
#'
#' # Get the field values#'
#' vocab.data <- bcdm_field_values(bold_search,
#'   specific.cols = c("inst", "identified_by")
#' )
#' }
#' @export
#'
#'
bcdm_field_values <- function(
  input.data,
  specific.cols,
  save.data = FALSE,
  output.file = NULL
) {
  # Allow both parquet path OR tbl_sql input
  if (inherits(input.data, "tbl_sql")) {
    bold_parquet_data <- input.data
  } else {
    bold_parquet_data <- import_parquet_data(input.data)
  }
  # Get unique values per column separately
  terms_list <- lapply(specific.cols, function(col) {
    bold_parquet_data %>%
      dplyr::filter(!is.na(.data[[col]]) & .data[[col]] != "") %>%
      dplyr::distinct(.data[[col]]) %>%
      dplyr::collect() %>%
      dplyr::pull(.data[[col]])
  })
  # Name the elements of the list as per the column names specified in the specific.cols argument
  names(terms_list) <- specific.cols
  # Save RDS file locally
  if (save.data) {
    if (is.null(output.file)) {
      stop("output.file must be provided when save.data = TRUE")
    }
    saveRDS(terms_list, paste0(output.file, ".rds"))
  }
  return(terms_list)
}
