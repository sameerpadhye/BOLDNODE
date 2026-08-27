#' Search records within a BOLD parquet data package
#'
#' @description Query records in BOLD parquet data packages using taxonomic, geographic, and other search criteria.
#'
#' @details This function loads the BOLD public data package parquet files (<https://boldsystems.org/data/data-packages/>) via DuckDB and applies filters based on the provided parameters. It supports filtering by ids (sampleid, processid), taxonomy (using a combination of scope.taxonomy and taxonomy), geography (using a combination of scope.geography and geography), biogeography (from biome to ecoregion level), BINs, institutes, identifiers, sequence sources, genetic markers, nucleotide base counts, dataset or projects, spatial bounding boxes, and ambiguous base percent cutoffs. The taxonomy and geography filters use a two-level scoping system that allows users to define `scope.taxonomy` and `scope.geography` to control how records are filtered. This is particularly useful when the same name exists at multiple geographic or taxonomic levels. For example, the name Azerbaijan may refer to a country, but a region with the same name also exists in Iran. If `scope.geography` is set to "any", records from both geographic levels will be returned. However, setting `scope.geography` to "country.ocean" restricts the search to country records only, returning records associated with the country of Azerbaijan. A similar situation can occur with taxonomic names where identical names are assigned to different groups. For example, the genus **Iris** occurs in both plants and animals. Using `scope.taxonomy` allows users to specify the desired taxonomic context and avoid ambiguity between groups. Users can also specify particular columns to return using the `specific.cols` parameter (column names can be checked using the `bcdm_field_names` function). The `tbl_sql` object can then be used by any  `bcdm_to_*` function for data transformations or `bold_search_collect` to load the query results in memory.
#'
#' @param input.parquet Path to the input parquet file.
#' @param ids Vector of process IDs or sample IDs used to filter records.
#' @param bins Vector of BIN numbers (i.e., URIs) used to filter records.
#' @param scope.taxonomy Character value specifying the kingdom. Values include "all", "Animalia", "Plantae", "Protista", "Fungi", "Bacteria" (default: NULL).
#' @param taxonomy Vector of taxonomic names to filter by. Values include "kingdom", "phylum", "class", "order", "family", "subfamily","genus", "species" (default: NULL).
#' @param scope.geography A character string specifying the geographic hierarchy level for a search. Values include "any", "country.ocean", "province.state", "region", "sector", "site" (default: "any").
#' @param geography Vector of geographic locations to filter by based on the `scope.geography`.
#' @param institutes Vector of institute codes used to filter records.
#' @param identified.by Vector of identifiers used to filter records.
#' @param seq.source Vector of sequence run sites used to filter records.
#' @param marker Vector of marker codes used to filter records.
#' @param basecount Nucleotide base count filter - either a single value or a vector of two values for a range.
#' @param biogeo.cat Character vector of biogeographic or ecological categories for filter by, such as a biome, realm, or ecoregion.
#' @param dataset.projects Vector of dataset/project codes used to filter records.
#' @param bounding.box Numeric vector of length 4: c(min_lon, max_lon, min_lat, max_lat).
#' @param ambi.base.cutoff Character value for filtering data based proportion of ambiguous bases (IUPAC codes). Valid values are "<1%", "1-5%", and ">5%" (default: NULL).
#' @param specific.cols Optional character vector of specific columns to return.
#'
#' @return A `tbl_sql` object containing the filtered data. The total number of records matching the search criteria is printed to the console.
#'
#' @importFrom dplyr filter select
#' @importFrom tools file_ext
#' @importFrom rlang !!
#' @importFrom utils head
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
#'
#' # Taxonomy
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   taxonomy = "Clytus"
#' )
#'
#' # Geography
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   geography = "Ontario"
#' )
#'
#' # Combination of many search criteria
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   taxonomy = "Trachysida",
#'   geography = "British Columbia",
#'   marker = "COI-5P",
#'   basecount = c(500, 660)
#' )
#' }
#'
#' @export
bold_parquet_search <- function(input.parquet,
                                ids = NULL,
                                bins = NULL,
                                scope.taxonomy = NULL,
                                taxonomy = NULL,
                                scope.geography = 'any',
                                geography = NULL,
                                institutes = NULL,
                                identified.by = NULL,
                                seq.source = NULL,
                                marker = NULL,
                                basecount = NULL,
                                biogeo.cat = NULL,
                                dataset.projects = NULL,
                                bounding.box = NULL,
                                ambi.base.cutoff = NULL,
                                specific.cols = NULL) {
  # checking the data type
  if (tolower(tools::file_ext(input.parquet)) != "parquet") {
    stop("Error: Input file must be a parquet")
  }
  # Import the parquet data
  parquet_data <- tryCatch(
    import_parquet_data(input.parquet),
    error = function(e) {
      stop(
        "Error: Failed to import parquet file '", input.parquet, "'. ",
        "Details: ", conditionMessage(e)
      )
    }
  )
  # Empty parquet data check
  if (nrow(parquet_data %>% head(1) %>% collect()) == 0) stop("Error: Parquet data is empty")
  # mapping the fetch.filter arguments with the bold.search function arguments
  bold.search.args <- list(
    ids = ids,
    bins = bins,
    taxonomy.kingdom = scope.taxonomy,
    taxonomy.names = taxonomy,
    geography.level = scope.geography,
    geography.names = geography,
    institutes = institutes,
    identified.by = identified.by,
    seq.source = seq.source,
    marker = marker,
    basecount = basecount,
    biogeo_cat = biogeo.cat,
    dataset.projects = dataset.projects,
    bounding.box = bounding.box,
    ambi.base.cutoff = ambi.base.cutoff
  )
  # Filter out NULL values and get their values
  null_args <- sapply(
    bold.search.args,
    is.null
  )
  # Selecting the non null arguments
  non_null_args <- bold.search.args[!null_args]
  # add the parquet data to the list
  non_null_args$bold.df <- parquet_data
  # Apply all the search filters on the non null arguments
  result <- do.call(
    bold.search.filters,
    non_null_args
  )
  # If specific columns are provided
  if (!is.null(specific.cols)) {
    bcdm_fields <- bcdm_field_names()
    # To check if the column names are BCDM
    if (any(!specific.cols %in% bcdm_fields$field)) stop("Please re-check the column names")
    result <- result %>%
      dplyr::select(all_of(specific.cols))
    result
  } else {
    result
  }
  # To provide the total records available in the search
  tot_records <- result %>%
    summarise(Total_records = n()) %>%
    collect()
  message(paste("The search has", tot_records, "records in the dataset", sep = " "))
  return(invisible(result))
}
