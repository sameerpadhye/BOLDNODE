# BOLDNODE 1.0.2
## Bug Fixes
- `bold_parquet_search()`, `bcdm_field_values()`: Data package paths containing an apostrophe no longer break the search.
- `bold_parquet_search()`: `dataset.projects` now matches whole codes (`DS-ABC` no longer matches `DS-ABCD`) and codes with special characters no longer cause errors.
- `get_bin_consensus()`, `get_bin_reps()`: Now stop with a clear message when the data has no records with BINs.
## Improvements
- `bcdm_field_names()`, `bcdm_to_dwc()`: BCDM reference tables are now shipped with the package and used when offline (a cached copy is refreshed from GitHub weekly). Set `options(BOLDNODE.offline = TRUE)` to skip the online check.
## Documentation
- `Using_BOLDNODE` vignette: Skip the `muscle` alignment and NJ tree chunks unless the required packages are available.
- `Using_BOLDNODE` vignette: The `muscle` alignment and NJ tree chunks now also require the `MultipleAlignment` package (needed by `Biostrings::DNAMultipleAlignment()` in Bioconductor 3.24 or later), fixing the vignette re-build error on r-devel.

# BOLDNODE 1.0.1
## Bug Fixes
- `bcdm_to_occmatrix()`: Fixed `pre.abs = TRUE` by correcting the object used for presence/absence conversion.
- `get_bin_consensus()`: Added the missing `data.table::setorderv` import.
- `get_bin_consensus()`: The function no longer throws an error when all identifications in the first group are empty.
- `get_bin_consensus()`: Previously, a discordant identification at the highest rank would sometimes also report a concordant identification for a lower one. This has been corrected.
## Compatibility
- `bold_parquet_search()`: Replaced deprecated DuckDB lambda syntax in the `ambi.base.cutoff` filter.
- `DESCRIPTION`: Set the minimum DuckDB version to `1.5.5`.
## Documentation
- `Using_BOLDNODE` vignette: Skip code chunks requiring `Biostrings` and `muscle` when these packages are unavailable, preventing CRAN build failures.
