#' Convert the BOLD parquet search in to a Darwin Core (DwC) format data model
#'
#' @description Converts `bold_parquet_search` results from BCDM format to Darwin Core Standard format.
#'
#' @details This function maps BCDM (Barcode Core Data Model) fields to their Darwin Core equivalents (<https://gbif.github.io/dwc-dp/qrg/>)
#'
#' \emph{Important Note}: All fields should be available in the `bold_parquet_search` `tbl_sql` object otherwise, the function will throw an error.
#'
#' @param bold.search.res A `tbl_sql` object containing BOLD search results.
#'
#' @return A data frame with columns mapped to Darwin Core equivalent fields.
#'
#' @importFrom dplyr select mutate case_when filter rename collect any_of
#' @importFrom data.table fread
#' @importFrom dbplyr sql
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
#'   taxonomy = "Lepturinae",
#'   marker = "COI-5P",
#'   basecount = c(550, 660)
#' )
#'
#' # Get the DwC object#'
#' bold.dwc <- bcdm_to_dwc(bold_search)
#' }
#' @export

bcdm_to_dwc <- function(bold.search.res) {
  # Check if the input is a 'tbl_sql' input
  check.tbl.sql(bold.search.res)
  # Check if all BCDM is available
  bcdm_fields <- bcdm_field_names()
  if (length(setdiff(
    bcdm_fields$field,
    colnames(bold.search.res)
  )) > 0) {
    stop("Please check if all the BCDM fields are present in the data")
  }
  # BOLD to BCDM mapping with certain updates for some field name mappings
  bcdm_path <- "https://raw.githubusercontent.com/boldsystems-central/BCDM/refs/heads/main/mapping_BCDM_to_DWC.tsv"

  bold.2.bcdm <- suppressMessages(data.table::fread(bcdm_path,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    verbose = FALSE,
    showProgress = FALSE,
    data.table = FALSE,
    fill = TRUE,
    tmpdir = tempdir()
  )) %>%
    select(
      bcdm_field,
      dwc_field,
      dwc_type
    ) %>%
    # change the country and province names for ease of use ('/' against '.')
    dplyr::mutate(bcdm_field = dplyr::case_when(
      bcdm_field == "country/ocean" ~ "country.ocean",
      bcdm_field == "province/state" ~ "province.state",
      TRUE ~ bcdm_field
    ))
  # Only the dwc fields are retained
  bold.2.bcdm <- bold.2.bcdm %>%
    dplyr::mutate(dwc_field = dplyr::case_when(
      bcdm_field == "processid" ~ "recordNumber",
      bcdm_field == "sampleid" ~ "materialEntityID",
      bcdm_field == "taxid" ~ "taxonID",
      bcdm_field == "species" ~ "acceptedNameUsage",
      bcdm_field == "identification" ~ "verbatimIdentification",
      bcdm_field == "identification_rank" ~ "verbatimTaxonRank",
      bcdm_field == "voucher type" ~ "disposition",
      bcdm_field == "tissue type" ~ "materialEntityType",
      bcdm_field == "specimen_linkout" ~ "associatedMedia",
      bcdm_field == "collection_notes" ~ "eventRemarks",
      bcdm_field == "geoid" ~ " locationId",
      bcdm_field == "site" ~ "locality",
      bcdm_field == "site_code" ~ "locationID",
      bcdm_field == "nuc" ~ "measurementValue",
      bcdm_field == "insdc_acs" ~ "materialSampleID",
      bcdm_field == "funding_src" ~ "fundingAttribution",
      bcdm_field == "marker_code" ~ "target_subfragment",
      bcdm_field == "primers_forward" ~ "pcr_primer_forward",
      bcdm_field == "primers_reverse" ~ "pcr_primer_reverse",
      bcdm_field == "sovereign_list" ~ "rightsHolder",
      bcdm_field == "tissue_type" ~ "materialEntityType",
      bcdm_field == "voucher_type" ~ "disposition",
      TRUE ~ dwc_field
    )) %>%
    mutate(bcdm_field = case_when(
      bcdm_field == "recordset_code_arr" ~ "bold_recordset_code_arr",
      TRUE ~ bcdm_field
    )) %>%
    dplyr::filter(dwc_field != "")
  # Adding the latitude, longitude data from BCDM's coord field
  dwc_output <- bold.search.res %>%
    select(any_of(bold.2.bcdm$bcdm_field)) %>%
    mutate(
      coord_clean = regexp_replace(coord, "\\[|\\]|\\s", ""),
      decimalLatitude = sql("split_part(coord_clean, ',', 1)"),
      # Extract longitude = second part
      decimalLongitude = sql("replace(split_part(coord_clean, ',', 2), ']', '')")
    ) %>%
    select(-c(coord)) %>%
    collect() %>%
    select(-coord_clean)
  # match the names
  rename_vec <- setNames(bold.2.bcdm$bcdm_field, bold.2.bcdm$dwc_field)
  # Remove the coord name from the rename vector as lat and lon are already created above
  coord_col <- which(rename_vec == "coord")
  # Create a vector of renamed names
  rename_vec <- rename_vec[-c(coord_col[[1]])]
  # Apply the column name changes to the data
  dwc_output <- dwc_output %>%
    rename(!!!rename_vec) %>%
    rename("country" = "country/waterBody")
  return(dwc_output)
}
