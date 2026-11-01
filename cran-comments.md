## Submission

This is a new submission.

## Test environments

* macOS Tahoe 26.7 (aarch64), R 4.6.2: `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: "New submission". The maintainer is Andrew Elmore.

## Notes for the reviewer

* Nothing in the examples, tests or vignette needs the network. `ld_proxy()` and
  `region_plot(ld = "ldlink")` query the LDlink web service (https://ldlink.nih.gov),
  which needs a personal token; their tests use saved responses, and the test
  and the examples that call the service are skipped on CRAN or unless
  `LDLINK_TOKEN` is set.
* The bundled gene tables (`genes_grch37`, `genes_grch38`; Ensembl, 0.4 MB) are
  built by `data-raw/genes.R`.
