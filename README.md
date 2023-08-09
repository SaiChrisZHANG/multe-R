# multe-R

The R version of Multiple Treatment Effects regression with saturated group control, as in [Goldsmith-Pinkham, P., Hull, P., &amp; Kolesár, M. (2022)](https://www.nber.org/papers/w30108).

See vignette [multeR](multeR/doc/MulteR.pdf) for description of the package (available through `browseVignettes("multeR")` once package is installed), and the package manual for documentation of the package functions.

The author of this package acknowledge the support.

## Installation

You can install the current development version from Github:

```r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("SaiChrisZHANG/multe-R", subdir="multeR", build_vignettes=TRUE)
```
