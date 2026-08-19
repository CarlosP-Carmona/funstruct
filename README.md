# funstruct

<!-- badges: start -->
[![R-CMD-check](https://github.com/CarlosP-Carmona/funstruct/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CarlosP-Carmona/funstruct/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**funstruct** provides a unified workflow for trait-based ecology: build
trait spaces, evaluate their quality and dimensionality, estimate
functional structure using trait probability densities (TPD) and classical
indices, and partition functional diversity across nested scales.

funstruct supersedes and extends the
[TPD](https://cran.r-project.org/package=TPD) and
[funspace](https://cran.r-project.org/package=funspace) packages.

## The pipeline

```r
library(funstruct)

space <- fs_space(traits, method = "pca")   # 1. build the FULL space
# (or: fs_space(fs_dist(traits), "pcoa")    #    dissimilarities are always
#  -- fs_dist gives 0-1 bounded Gower)      #    computed explicitly, upfront
fs_dimensionality(space)                    # 2. evaluate before choosing dims
fs_quality(space)
space <- fs_reduce(space, dims = 3)         # 3. reduce (and optionally rotate)
space <- fs_rotate(space)                   #    varimax, PCA spaces only
space <- fs_tpd(space, obs = individuals,   # 4. trait probability densities
                ids = "species")
str_   <- fs_structure(space, comm)         # 5. functional structure
part   <- fs_partition(space, comm,         # 6. partition across scales
                       fs_hierarchy(plot, site))
```

Steps 4-6 are under active development; the current version implements
space building, import (`as_fspace()`), reduction and rotation.

## Design principles

- **Evaluate before you reduce.** `fs_space()` builds the full space; the
  workflow makes dimensionality an explicit, justified choice.
- **Two engines.** Indices are available from trait probability densities
  (`engine = "prob"`) and from point-based methods (`engine = "points"`),
  behind the same interface.
- **Dissimilarities are explicit and bounded.** `fs_space()` never computes
  dissimilarities internally; `fs_dist()` builds them from any mix of trait
  types (including intraspecific variability via distribution overlap), and
  every trait contributes a distance bounded 0-1, so Rao/MPD keep an
  interpretable scale.
- **Bandwidths are never silent.** Smoothing parameters are attached to the
  ecological entity they describe, stored in every object, printed in every
  summary, and stress-testable with `fs_sensitivity()`.
- **Base R.** No tidyverse dependencies; base graphics throughout.

## Installation

```r
# development version
# install.packages("remotes")
remotes::install_github("CarlosP-Carmona/funstruct")
```

## For developers

```r
# after cloning
devtools::document()   # regenerate man/ and NAMESPACE
devtools::test()       # run the test suite
devtools::check()      # full R CMD check
```
