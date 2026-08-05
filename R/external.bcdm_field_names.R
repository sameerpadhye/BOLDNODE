#' Retrieve BCDM field names and descriptions
#'
#' @description Provides information on the field (column) names and their respective data type, all of which are compliant with the Barcode Core Data Model (BCDM).
#'
#' @param print.output Whether the output should be printed in the console. Default is FALSE.
#'
#' @details The function downloads the latest field (column) metadata (file type and brief description) for the Barcode Core Data Model (BCDM) from <https://github.com/boldsystems-central/BCDM/blob/main/field_definitions.tsv>; output = TRUE will print the information in the console.
#' \emph{Important Note}: Two field names 'country/ocean' and 'province/state' have been modified to 'country.ocean' and 'province.state' to match with BOLDconnectR output and for operational ease.
#'
#' @returns A data frame containing definitions for all fields (columns).
#'
#' @examples
#'
#' bold.field.data <- bcdm_field_names()
#'
#' head(bold.field.data, 10)
#'
#' @importFrom dplyr matches case_when select mutate %>%
#' @importFrom data.table fread
#'
#' @export
#'
bcdm_field_names <- function(print.output = FALSE) {
  bold.fields.data <- suppressMessages(data.table::fread("https://raw.githubusercontent.com/boldsystems-central/BCDM/refs/heads/main/field_definitions.tsv",
    sep = "\t",
    quote = "",
    check.names = FALSE,
    verbose = FALSE,
    showProgress = FALSE,
    data.table = FALSE,
    fill = TRUE,
    tmpdir = tempdir()
  )) %>%
    dplyr::select(
      dplyr::matches("field", ignore.case = TRUE),
      dplyr::matches("definition", ignore.case = TRUE),
      dplyr::matches("data_type", ignore.case = TRUE)
    ) %>%
    dplyr::mutate(R_field_types = dplyr::case_when(
      data_type == "string" ~ "character",
      data_type %in% c("char", "array") ~ "character",
      data_type == "float" ~ "numeric",
      data_type == "number" ~ "numeric",
      data_type == "integer" ~ "integer",
      data_type == "string:date" ~ "Date"
    )) %>%
    dplyr::select(-dplyr::matches("data_type", ignore.case = TRUE)) %>%
    dplyr::mutate(field = case_when(
      field == "country/ocean" ~ "country.ocean",
      field == "province/state" ~ "province.state",
      TRUE ~ field
    ))
  # Print the fields info in console
  if (print.output == TRUE) {
    return(bold.fields.data)
  } else {
    # This is so that the whole output is not printed in the console
    invisible(bold.fields.data)
  }
}
