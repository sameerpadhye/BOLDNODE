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
