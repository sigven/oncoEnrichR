# Write oncoEnrichR report object to output file

Function that writes an oncoEnrichR report object to file, either as an
interactive HTML report or as an Excel workbook.

## Usage

``` r
write(
  report,
  oeDB,
  file = "testReport.html",
  ignore_file_extension = F,
  overwrite = F,
  render_quarto_quiet = T,
  format = "html",
  ...
)
```

## Arguments

- report:

  object with oncoEnrichR report data (returned by oeDB\$onco_enrich)

- oeDB:

  oncoEnrichR data repository object - as returned from
  [`load_db()`](https://sigven.github.io/oncoEnrichR/reference/load_db.md)

- file:

  full filename for report output (e.g. "oe_report.html" or
  "oe_report.xlsx")

- ignore_file_extension:

  logical to accept any type of filaname extensions (for Galaxy
  integration)

- overwrite:

  logical indicating if existing output files may be overwritten

- render_quarto_quiet:

  logical indicating if Quarto rendering should be done quietly

- format:

  file format of output (html/excel)

- ...:

  options for Galaxy/non self-contained HTML. Only applicable for use in
  Galaxy
