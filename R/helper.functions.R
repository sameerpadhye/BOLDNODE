#' @keywords internal
#' @importFrom dplyr filter select between tbl summarize
#' @importFrom DBI dbConnect
#' @importFrom tools file_ext
#' @importFrom rlang !!
#' @importFrom utils head

# 1 Import Parquet File
import_parquet_data <- function(path) {
  # 1 Establish a temporary connection
  temp_connection <- DBI::dbConnect(duckdb::duckdb())
  # 2 Create a query
  get_parquet_data <- sprintf("Select * from parquet_scan('%s')", path)
  parquet_data <- tbl(
    temp_connection,
    sql(get_parquet_data)
  ) %>%
    # change the country and province names for ease of use ('/' against '.')
    dplyr::rename(
      "country.ocean" = "country/ocean",
      "province.state" = "province/state"
    ) %>%
    mutate(
      coord = sql("replace(replace(trim(coord), '[', ''), ']', '')"),
      bold_recordset_code_arr = sql("replace(replace(replace(replace(trim(bold_recordset_code_arr),
                                    '[', ''),
                                    ']', ''),
                                    '''', ''),
                                    ', ', ',')")
    )
  return(parquet_data)
  # Ensure disconnection at the end, even if an error occurs
  on.exit(
    DBI::dbDisconnect(temp_connection,
      shutdown = TRUE
    ),
    add = TRUE
  )
}

# 2 Check tbl_sql object
check.tbl.sql <- function(tbl) {
  if (!inherits(tbl, c("tbl_sql", "tbl_dbi"))) {
    stop("Error: Input must be a bold_parquet_search output (tbl_sql / dbplyr table)")
  }
  TRUE
}

# 3 BOLD search filters
bold.search.filters <- function(bold.df,
                                ids = NULL,
                                bins = NULL,
                                taxonomy.kingdom = NULL,
                                taxonomy.names = NULL,
                                geography.level = "any",
                                geography.names = NULL,
                                institutes = NULL,
                                identified.by = NULL,
                                seq.source = NULL,
                                marker = NULL,
                                basecount = NULL,
                                biogeo_cat = NULL,
                                dataset.projects = NULL,
                                bounding.box = NULL,
                                ambi.base.cutoff = NULL) {
  # 1 ids
  # ids is your vector of IDs to filter
  if (!is.null(ids)) {
    bold.df <- bold.df %>%
      filter(processid %in% !!ids |
        sampleid %in% !!ids)
  }
  # 2 BIN
  if (!is.null(bins)) {
    bold.df <- bold.df %>%
      filter(bin_uri %in% !!bins)
  }
  # 3 taxon name
  # condition to check if the taxon name is of the correct data type

  if (!is.null(taxonomy.names)) {
    # Restrict to a kingdom only if requested
    if (!is.null(taxonomy.kingdom)) {
      bold.df <- bold.df %>%
        dplyr::filter(kingdom %in% taxonomy.kingdom)
    } else {
      bold.df
    }

    # Apply taxonomy name filter
    bold.df <- bold.df %>%
      dplyr::filter(
        phylum %in% taxonomy.names |
          class %in% taxonomy.names |
          `order` %in% taxonomy.names |
          family %in% taxonomy.names |
          subfamily %in% taxonomy.names |
          genus %in% taxonomy.names |
          species %in% taxonomy.names
      )
  }

  # 4 specific country/region/site/sector

  if (!is.null(geography.names)) {
    # Restrict to a geography level only if requested
    if (geography.level != "any") {
      if (geography.level == "country.ocean") {
        bold.df <- bold.df %>%
          dplyr::filter(country.ocean %in% geography.names)
      } else if (geography.level == "province.state") {
        bold.df <- bold.df %>%
          dplyr::filter(province.state %in% geography.names)
      } else if (geography.level == "region") {
        bold.df <- bold.df %>%
          dplyr::filter(region %in% geography.names)
      } else if (geography.level == "sector") {
        bold.df <- bold.df %>%
          dplyr::filter(sector %in% geography.names)
      } else if (geography.level == "site") {
        bold.df <- bold.df %>%
          dplyr::filter(site %in% geography.names)
      }
    } else {
      bold.df <- bold.df %>%
        dplyr::filter(
          country.ocean %in% geography.names |
            province.state %in% geography.names |
            region %in% geography.names |
            sector %in% geography.names |
            site %in% geography.names
        )
    }
  }
  # 5 Latitude/Longitude bounding box
  if (!is.null(bounding.box)) {
    if (!is.numeric(bounding.box) || length(bounding.box) != 4) {
      stop(
        "Error: 'bounding.box' must be a numeric vector of length 4: ",
        "c(min_lon, max_lon, min_lat, max_lat)."
      )
    }
    min_lon <- bounding.box[1]
    max_lon <- bounding.box[2]
    min_lat <- bounding.box[3]
    max_lat <- bounding.box[4]
    bold.df <- bold.df %>%
      mutate(
        coord_clean = regexp_replace(coord, "\\[|\\]|\\s", ""),
        lat = sql("split_part(coord_clean, ',', 1)"),
        # Extract longitude = second part
        lon = sql("replace(split_part(coord_clean, ',', 2), ']', '')"),
        lat = as.numeric(lat),
        lon = as.numeric(lon)
      ) %>%
      filter(
        between(lon, min_lon, max_lon),
        between(lat, min_lat, max_lat)
      ) %>%
      select(-c(coord_clean, lat, lon))
  }
  # 6 Institutes storing the specimen
  if (!is.null(institutes)) {
    bold.df <- bold.df %>%
      dplyr::filter(inst %in% !!institutes)
  }
  # 7 Identified by
  if (!is.null(identified.by)) {
    bold.df <- bold.df %>%
      dplyr::filter(identified_by %in% !!identified.by)
  }
  # 8 sequence source
  if (!is.null(seq.source)) {
    bold.df <- bold.df %>%
      dplyr::filter(sequence_run_site %in% !!seq.source)
  }
  # 9 Type of marker
  if (!is.null(marker)) {
    bold.df <- bold.df %>%
      dplyr::filter(marker_code %in% !!marker)
  }
  # 10 basecount
  if (!is.null(basecount)) {
    # if user has selected multiple markers, lengths should be provided as a named list
    if (is.list(basecount)) {
      for (i in seq_along(basecount)) {
        marker <- names(basecount[i])
        vals <- unname(unlist(basecount[i]))
        if (length(vals) == 1) {
          bold.df <- bold.df %>%
            dplyr::filter(((marker_code == marker) & (nuc_basecount %in% vals)) |
              (marker_code != marker))
        } else if (length(vals) == 2) {
          first_val <- vals[1]
          last_val <- vals[2]
          bold.df <- bold.df %>%
            dplyr::filter(((marker_code == marker) &
              (dplyr::between(nuc_basecount, first_val, last_val))) |
              (marker_code != marker))
        }
      }
    } else {
      if (length(basecount) == 1) {
        bold.df <- bold.df %>%
          dplyr::filter(nuc_basecount %in% basecount)
      } else if (length(basecount) == 2) {
        first_val <- basecount[1]
        last_val <- basecount[2]
        bold.df <- bold.df %>%
          dplyr::filter(dplyr::between(
            nuc_basecount,
            first_val,
            last_val
          ))
      } else {
        stop("Incorrect value input")
      }
    }
  }
  # 11 Biogeo/ecological categories
  if (!is.null(biogeo_cat)) {
    # filter condition to select the specific taxon name.
    bold.df <- bold.df %>%
      dplyr::filter(biome %in% !!biogeo_cat |
        realm %in% !!biogeo_cat |
        ecoregion %in% !!biogeo_cat)
  }
  # 12 Dataset or project code
  if (!is.null(dataset.projects)) {
    codes_regex <- paste(dataset.projects, collapse = "|")
    bold.df <- bold.df %>%
      mutate(
        bold_recordset_code_arr = regexp_replace(bold_recordset_code_arr, "\\[|\\]", "")
      ) %>%
      filter(regexp_matches(
        bold_recordset_code_arr,
        codes_regex
      ))
  }
  # 13 Ambiguous bases cutoff
  # Return unchanged if no cutoff provided
  if (!is.null(ambi.base.cutoff)) {
    bold.df <- bold.df %>%
      mutate(
        count_ambi = sql("
        ARRAY_LENGTH(
          ARRAY_FILTER(
            REGEXP_SPLIT_TO_ARRAY(UPPER(nuc), ''),
            LAMBDA x: REGEXP_MATCHES(x, '[NRYSWKMBDHV]')
          ))"),
        percent_ambi = (count_ambi * 100.0 / nuc_basecount)
      ) %>%
      filter(case_when(
        percent_ambi < 1.00 ~ "<1%",
        percent_ambi >= 1.00 & percent_ambi <= 5.00 ~ "1-5%",
        percent_ambi > 5.00 ~ ">5%"
      ) == !!ambi.base.cutoff)
  }
  return(bold.df)
}

# 4 Chunked data indices

get_chunk_indices <- function(input_file,
                              chunksize = 100000) {
  chunk_size <- chunksize
  # Determine total rows
  total_rows <- input_file %>%
    summarize(n = n()) %>%
    collect() %>%
    pull(n)
  # Calculating the number of chunks required
  n_chunks <- ceiling(total_rows / chunk_size)
  # Generate chunk indices
  chunk_indices <- seq_len(n_chunks)
  chunk_res <- list(
    total_rows = total_rows,
    chunk_size = chunk_size,
    chunk_indices = chunk_indices
  )
  return(chunk_res)
}

# 5 Regular expression for interim names

re_int <- "\\.\\Z|\\S{2,}\\.\\S|[0-9]|\\s[A-Z]|[A-Z]\\Z|[A-Z]{2}|[a-z][A-Z]|_(?!(hn|sl|ss)\\Z)|%|\\?|!|\\[|\\]|\\{|\\}|\\(|\\)|,|\\s(?:aff|agg|cf|complex|group|grp|gr|gp|cmplx|pr|ms|cfr|nr|nsp|near|nomen|hybrid|voucher|form|from|ss|ssl|see|spp?|sample)\\.?(?:\\s|\\Z)"

# 6 Select BIN reps from `tbl_sql` using DBI & dbplyr

get_sql_sort <- function(df, criteria) {
  make_sort_key <- function(col, levels = NULL, grepl_pat = NULL) {
    col_sym <- sym(col)

    key <- if (!is.null(levels)) {
      cases <- lapply(seq_along(levels), function(i) {
        dplyr::expr(!!col_sym == !!levels[[i]] ~ !!i)
      })
      dplyr::expr(case_when(!!!cases))
    } else if (!is.null(grepl_pat)) {
      dplyr::expr(grepl(!!grepl_pat, !!col_sym))
    } else {
      col_sym
    }

    key
  }

  sort_sequence <- lapply(setNames(names(criteria), names(criteria)), function(step) {
    if (step == "vouchered") {
      list(key = make_sort_key(grepl_pat = "(?i)GenBank", col = "inst"))
    } else if (step == "seq_length") {
      if (criteria$seq_length == "COI_auto") {
        list(
          prep = function(df) {
            df %>%
              dplyr::mutate(
                .rounded = (round((nuc_basecount - 1) / 3) * 3) + 1
              ) %>%
              dplyr::group_by(bin_uri, .rounded) %>%
              dplyr::mutate(
                .freq = n()
              ) %>%
              dplyr::group_by(bin_uri) %>%
              dplyr::mutate(
                .tiebreak = abs(.rounded - 658)
              ) %>%
              dbplyr::window_order(desc(.freq), .tiebreak) %>%
              dplyr::mutate(
                .mode_rank = dplyr::row_number(),
                .mode = dplyr::if_else(.mode_rank == 1L, .rounded, NA_real_),
                .bin_mode = max(.mode, na.rm = TRUE)
              ) %>%
              dplyr::ungroup() %>%
              dplyr::select(-.rounded, -.freq, -.tiebreak, -.mode_rank, -.mode)
          },
          key = dplyr::expr(abs(nuc_basecount - .bin_mode))
        )
      } else if (criteria$seq_length == "longest") {
        list(key = make_sort_key(col = "nuc_basecount"), desc = TRUE)
      } else if (criteria$seq_length == "shortest") {
        list(key = make_sort_key(col = "nuc_basecount"), desc = FALSE)
      } else {
        list(key = dplyr::expr(abs(nuc_basecount - round(criteria$seq_length, 0))))
      }
    } else if (step == "id_method") {
      list(key = make_sort_key(
        col = "identification_method",
        levels = criteria$id_method
      ))
    } else if (step == "inst") {
      list(key = make_sort_key(
        col = "inst",
        levels = criteria$inst
      ))
    } else if (step == "coll_date") {
      list(
        key = if (criteria$coll_date == "latest") {
          key <- make_sort_key(col = "collection_date_start")
          dplyr::expr(dplyr::coalesce(as.Date(!!key), as.Date("0001-01-01")))
        } else {
          key <- make_sort_key(col = "collection_date_start")
        },
        desc = criteria$coll_date == "latest"
      )
    } else if (step == "seq_date") {
      list(
        key = if (criteria$coll_date == "latest") {
          key <- make_sort_key(col = "sequence_upload_date")
          dplyr::expr(dplyr::coalesce(as.Date(!!key), as.Date("0001-01-01")))
        } else {
          key <- make_sort_key(col = "sequence_upload_date")
        },
        desc = criteria$seq_date == "latest"
      )
    }
  })

  sort_sequence
}

# 7 Select BIN reps from data in memory using data.table

get_dt_sort <- function(dt, criteria) {
  lapply(names(criteria), function(step) {
    if (step == "vouchered") {
      (grepl("(?i)GenBank", dt$inst, perl = T) == criteria$vouchered)
    } else if (step == "seq_length") {
      seq_sort <- if (criteria$seq_length == "COI_auto") {
        # find the modal sequence length for each BIN
        mode_by_bin <- dt[, .(mode = {
          t <- table(sapply(nuc_basecount, function(bp) (round((bp - 1) / 3) * 3) + 1))
          modal <- names(t)[t == max(t)]

          # in case of multiple modes, select the one nearest to 658bp
          as.numeric(unname(modal)[which.min(sapply(modal, function(m) abs(as.numeric(m) - 658)))])
        }), by = "bin_uri"]

        # output absolute difference between sequence length and modal length for BIN
        dt[, .SD, .SDcols = c("bin_uri", "nuc_basecount")][mode_by_bin, on = "bin_uri"][, ({
          abs(nuc_basecount - mode)
        })]
      } else if (criteria$seq_length == "longest") {
        dt[["nuc_basecount"]] * -1
      } else if (criteria$seq_length == "shortest") {
        dt[["nuc_basecount"]]
      } else {
        abs(dt[["nuc_basecount"]] - round(criteria$seq_length, 0))
      }
    } else if (step == "id_method") {
      factor(as.character(dt$identification_method), levels = criteria$id_method)
    } else if (step == "inst") {
      factor(as.character(dt$inst), levels = criteria$inst)
    } else if (step == "coll_date") {
      dt$collection_date_start
    } else if (step == "seq_date") {
      dt$sequence_upload_date
    }
  })
}
