#' Collect and export parquet search results
#'
#' @description Collects, outputs and exports the results of a `bold_parquet_search` query, processing large datasets in user-defined chunks to improve memory efficiency.
#' @details This function collects the results of a `bold_parquet_search` query into the current R session. To facilitate the handling of large datasets, records can be processed in user-defined chunks, with optional pauses between chunks to help manage memory usage and system resources. The function also supports exporting results in TSV or parquet format. When export = FALSE (default), the collected data are returned only within the R session. When export = TRUE, a complete file path, including a file name and extension, must be provided via output.path.
#' \emph{Important Note}: Some queries (for example, all records from the order Diptera) may produce very large result sets that exceed the available RAM on lower-specification systems (e.g., 8 GB RAM), regardless of the chunking and system sleep settings.
#' @param bold.search.res A `tbl_sql` object obtained from `bold_parquet_search`.
#' @param chunk.size Maximum number of rows to process in each chunk (default: 1e6).
#' @param sys.sleep Time to sleep between chunks in seconds (default: 0).
#' @param export Logical value that allows user to export the output locally (default: FALSE).
#' @param export.type Character string specifying the data type of the exported file (tsv or parquet). Required when export=TRUE.
#' @param output.path Character string specifying the local path for data export along with the file name and extension. Required when export=TRUE.
#'
#' @return A data frame containing all collected results. If export = TRUE, the results are also exported locally as either a TSV or Parquet file.
#' @importFrom dplyr summarise collect pull bind_rows %>%
#' @importFrom DBI dbExecute
#' @importFrom dbplyr remote_con sql_render
#' @importFrom progressr with_progress progressor
#' @examples
#' \dontrun{
#'
#'
#' # Search the BOLD data package
#' bold_search <- bold_parquet_search(
#'   input.parquet = parquet_file,
#'   taxonomy = "Coleoptera",
#'   geography = "Canada",
#'   marker = "COI-5P",
#'   basecount = c(500, 660)
#' )
#'
#' # Collect the data  (no export)
#' bold_search_collect(
#'   bold_search,
#'   chunk.size = 50000,
#'   export = FALSE
#' )
#'
#' # Collect and export
#' bold_search_collect(
#'   bold_search,
#'   chunk.size = 50000,
#'   export = TRUE,
#'   export.type = "parquet",
#'   output.path = "path/to/output_file.parquet"
#' )
#'
#' }
#' @export

bold_search_collect <- function(
  bold.search.res,
  chunk.size = 1000000,
  sys.sleep = 0,
  export = FALSE,
  export.type = c("tsv", "parquet"),
  output.path = NULL
) {
  # Match the type of file for export
  export.type <- match.arg(export.type)
  # checking the data type
  check.tbl.sql(bold.search.res)
  # establish a temporary connection (for duckdb)
  con <- dbplyr::remote_con(bold.search.res)
  # parquet export (direct and doesnt need chunks)
  if (export && export.type == "parquet") {
    if (is.null(output.path)) {
      stop("Please provide output.path for parquet export")
    }
    # Render the tbl_sql into a SQL like object for using to export the collected result back as a parquet file
    query <- dbplyr::sql_render(bold.search.res)
    tryCatch(
      DBI::dbExecute(
        con,
        paste0("COPY (", query, ") TO '", output.path, "'(FORMAT PARQUET, COMPRESSION ZSTD)")
      ),
      error = function(e) {
        stop(
          "Error: Parquet export failed for '", output.path, "'. ",
          "Details: ", conditionMessage(e)
        )
      }
    )
    message("Parquet export complete.")
    return(invisible(NULL))
  }
  # TSV export needs chunking of data for collection before export
  # 1 Getting chunks
  # Disable progress bar
  DBI::dbExecute(con, "PRAGMA disable_progress_bar;")
  # chunking the data
  chunk_info <- get_chunk_indices(
    input_file = bold.search.res,
    chunksize = chunk.size
  )
  total_rows <- chunk_info$total_rows
  chunk_size <- chunk_info$chunk_size
  chunk_indices <- chunk_info$chunk_indices
  total_chunks <- length(chunk_indices)
  # Enable the progress bar
  DBI::dbExecute(con, "PRAGMA enable_progress_bar;")
  # Collecting based on the number of chunks
  if (total_chunks == 1) {
    message(sprintf("Collecting all %d rows in a single chunk...", total_rows))
    res <- bold.search.res %>% dplyr::collect()
  } else {
    tbl_sql <- dbplyr::sql_render(bold.search.res)
    res_chunks <- progressr::with_progress({
      p <- progressr::progressor(steps = total_chunks)
      lapply(chunk_indices, function(i) {
        offset <- (i - 1) * chunk_size
        size <- min(chunk_size, total_rows - offset)
        p(sprintf("Chunk %d/%d (%d rows)", i, total_chunks, size))
        sql_query <- paste0(
          "SELECT * FROM (",
          tbl_sql,
          ") AS sub_tbl ",
          "LIMIT ",
          size,
          " OFFSET ",
          offset
        )
        out <- tryCatch(
          DBI::dbGetQuery(con, sql_query),
          error = function(e) {
            stop(
              "Error: Failed to collect chunk ", i, "/", total_chunks, ". ",
              "Details: ", conditionMessage(e)
            )
          }
        )
        # If sys.sleep is provided
        if (sys.sleep > 0) {
          Sys.sleep(sys.sleep)
        }
        out
      })
    })
    # Combine all the chunks
    res <- dplyr::bind_rows(res_chunks)
  }
  # Exporting as a TSV
  if (export && export.type == "tsv") {
    if (is.null(output.path)) {
      stop("Please provide output.path for TSV export")
    }
    utils::write.table(
      res,
      file = output.path,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    message("TSV export complete.")
  }

  return(invisible(res))
}
