#' Generate a concise summary of the search results
#'
#' @description Creates a summary statistics table from the `bold_parquet_search` `tb_sql` object.
#'
#' @details The function provides a concise summary of the search obtained by `bold_parquet_search` that includes: total records, unique BINs, unique institutes, unique markers and amplicon size range.
#'
#' @param bold.search.res A `tbl_sql` object containing `bold_parquet_search` results.
#'
#' @return A data frame with the summary statistics.
#'
#' @importFrom dplyr summarise n n_distinct case_when pull mutate across select collect %>% everything
#' @importFrom tidyr pivot_longer
#' @importFrom DBI dbExecute
#' @importFrom dbplyr remote_con
#'
#' @examples
#' \dontrun{
#'
#' # Search the BOLD data package
#'
#' parquet_file <- "user-defined path to parquet file"
#'
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   taxonomy = "Hemiptera",
#'   geography = "India",
#'   marker = "COI-5P"
#' )
#' # Get the concise summary
#' bold_summary <- get_concise_summary(bold_search)
#' }
#' @export

get_concise_summary <- function(bold.search.res) {
  # checking the data type
  check.tbl.sql(bold.search.res)
  bold.search.res.cols <- bold.search.res %>% colnames()
  # Getting all the BCDM field column names
  bold_field_data <- bcdm_field_names(print.output = F) %>%
    dplyr::select(field)
  # Check to see if all BCDM fields are present in the input data
  if (!all(bold_field_data$field %in% bold.search.res.cols)) {
    stop("Error: Concise summary requires all BCDM fields. Please re-check the search criteria.")
  }
  con <- dbplyr::remote_con(bold.search.res)
  # A SQL like query to get the concise summary
  concise_summary <- bold.search.res %>%
    summarise(
      Total_records = n(),
      Total_records_w_sequences = sum(!is.na(nuc)),
      Unique_species = n_distinct(species, na.rm = TRUE),
      Unique_species_w_BINs = n_distinct(
        case_when(!is.na(bin_uri) ~ species),
        na.rm = TRUE
      ),
      Unique_BINs = n_distinct(bin_uri, na.rm = TRUE),
      Unique_countries = n_distinct(country.ocean, na.rm = TRUE),
      Unique_institutes = n_distinct(inst, na.rm = TRUE),
      Unique_identified_by = n_distinct(identified_by, na.rm = TRUE),
      Unique_specimen_depositories = n_distinct(sequence_run_site, na.rm = TRUE),
      min_amplicon = min(nuc_basecount, na.rm = TRUE),
      max_amplicon = max(nuc_basecount, na.rm = TRUE)
    ) %>%
    collect()
  DBI::dbExecute(con, "PRAGMA disable_progress_bar;")
  # Further edits in the collected data
  concise_summary <- concise_summary %>%
    mutate(
      Unique_markers = paste(unique(bold.search.res %>%
        distinct(marker_code) %>%
        collect() %>%
        pull()), collapse = ","),
      Amplicon_length_range = paste(min_amplicon, "-", max_amplicon),
      across(everything(), as.character)
    ) %>%
    select(-min_amplicon, -max_amplicon) %>%
    tidyr::pivot_longer(
      everything(),
      names_to = "Category",
      values_to = "Value"
    )
  DBI::dbExecute(con, "PRAGMA enable_progress_bar;")
  return(concise_summary)
}
