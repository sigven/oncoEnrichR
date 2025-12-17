# Load oncoEnrichR data repository

Function that fetches version-tagged oncoEnrichR data repository from
Google Drive

## Usage

``` r
load_db(cache_dir = NA, force_download = FALSE, googledrive = FALSE)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should force download (i.e. set to
  TRUE to re-download even if data exists in cache)

- googledrive:

  logical indicating to use googledrive or UiO server as host for
  downloading

## Value

A `list` object with oncoEnrichR datasets, to be used as the *oeDB*
argument for
[`onco_enrich()`](https://sigven.github.io/oncoEnrichR/reference/onco_enrich.md)
and [`write()`](https://sigven.github.io/oncoEnrichR/reference/write.md)
