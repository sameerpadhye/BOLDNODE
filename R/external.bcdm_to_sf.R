#' Convert the BOLD parquet search into a simple features dataframe
#'
#' @description Converts BOLD search results with coordinate data to a simple features dataframe with point geometry (sf object).
#'
#' @details This function transforms the search results from `bold_parquet_search` into an `sf` object. A data chunking  option is available to manage large sizes to avoid memory issues. The function creates point geometries in the WGS84 coordinate system (EPSG:4326). Records that don't have coordinate data are removed during processing.
#'
#' @param bold.search.res A tbl_sql object containing BOLD search results.
#' @param chunk.size Number of records to process in each chunk (default: 100000).
#'
#' @return An `sf` object with point geometry in the WGS84 coordinate reference system (EPSG:4326).
#'
#' @importFrom dplyr filter mutate select collect
#' @importFrom sf st_as_sf
#' @importFrom dbplyr sql
#' @importFrom dplyr compute
#'
#' @examples
#' \dontrun{
#'
#' # Search the BOLD data package
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   taxonomy = "Odonata",
#'   geography = "Malaysia"
#' )
#'
#' # Get the occurrence matrix#'
#' sf_data <- bcdm_to_sf(bold_search, chunk.size = 100000)
#' }
#'
#' @export
bcdm_to_sf <- function(bold.search.res, chunk.size = 100000) {
  # checking the data type
  check.tbl.sql(bold.search.res)
  # Check if 'coord' field is in the input data
  if (!"coord" %in% colnames(bold.search.res)) {
    stop("'coord' column not found. 'coord' column is required for generating 'sf' object")
  }
  # establish a temporary connection (for duckdb)
  con <- dbplyr::remote_con(bold.search.res)
  # duckDB progress bar disabled
  DBI::dbExecute(con, "PRAGMA disable_progress_bar;")
  # Chunk indexes
  indices <- get_chunk_indices(
    input_file = bold.search.res,
    chunksize = chunk.size
  )
  # duckDB progress bar enabled
  DBI::dbExecute(con, "PRAGMA enable_progress_bar;")
  #  sql query
  geo_data <- bold.search.res %>%
    dplyr::filter(!is.na(coord)) %>%
    dplyr::mutate(
      coord_clean = sql("replace(replace(trim(coord), '[', ''), ']', '')"),
      # instr here gives the position of the substring. Here it checks that the string contains a comma before splitting and casting as a double
      lat = sql("CASE WHEN instr(replace(replace(trim(coord), '[', ''), ']', ''), ',') > 0
          THEN CAST(split_part(replace(replace(trim(coord), '[', ''), ']', ''), ',', 1) AS DOUBLE)
          ELSE NULL END"),
      lon = sql("CASE WHEN instr(replace(replace(trim(coord), '[', ''), ']', ''), ',') > 0
          THEN CAST(split_part(replace(replace(trim(coord), '[', ''), ']', ''), ',', 2) AS DOUBLE)
          ELSE NULL END"),
      row_num = sql("row_number() over (order by coord)")
    ) %>%
    dplyr::compute()
  # Using the chunk index for converting the data into sf object
  sf_data_list <- lapply(indices$chunk_indices, function(i) {
    start <- (i - 1) * indices$chunk_size + 1
    end <- min(i * indices$chunk_size, indices$total_rows)
    geo_data %>%
      dplyr::filter(row_num >= start & row_num <= end) %>%
      dplyr::select(-row_num, coord_clean) %>%
      filter(!is.na(coord_clean)) %>%
      dplyr::collect() %>%
      dplyr::filter(!is.na(lat) | !is.na(lon)) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  })
  sf_data <- dplyr::bind_rows(sf_data_list)
  if (is.data.frame(sf_data) && nrow(sf_data) == 0) stop("No data retrieved.Please re-check the search criteria.")
  return(sf_data)
}
