# Package datasets -------------------------------------------------------------

#' Mediterranean grassland: individual-level traits
#'
#' Individual-level trait measurements from a Mediterranean annual
#' grassland (Valdeloshielos, central Spain). In each of 40 plots, up to
#' ten individuals of the most abundant species were measured. Only
#' individuals with both traits measured are included.
#'
#' @format A data frame with 2,444 rows and 5 columns:
#' \describe{
#'   \item{plot}{Plot identifier (integer, 1-40).}
#'   \item{species}{Species name.}
#'   \item{individual}{Individual number within plot and species (1-10).}
#'   \item{height}{Vegetative height (cm).}
#'   \item{sla}{Specific leaf area (mm2 mg-1).}
#' }
#' @source Carmona, C.P., Rota, C., Azcarate, F.M. & Peco, B. (2015)
#'   More for less: sampling strategies of plant functional traits across
#'   local environmental gradients. *Functional Ecology*, 29, 579-588.
#'   \doi{10.1111/1365-2435.12366}. See also [grassland_abund] for the
#'   corresponding relative covers.
"grassland"

#' Mediterranean grassland: relative covers
#'
#' Relative cover of each species in each of the 40 plots of the
#' Valdeloshielos grassland (rows sum to 1). Companion to [grassland].
#'
#' @format A numeric matrix with 40 rows (plots) and 51 columns (species).
#' @source Carmona, C.P., Rota, C., Azcarate, F.M. & Peco, B. (2015)
#'   *Functional Ecology*, 29, 579-588. \doi{10.1111/1365-2435.12366}.
"grassland_abund"

#' Global spectrum of plant form and function: complete trait data
#'
#' Six aboveground traits for 2,630 plant species with complete trait
#' information, from the global spectrum of plant form and function
#' (Diaz et al. 2016), as compiled from the TRY database and distributed
#' with the \pkg{funspace} package. All traits are log10-transformed and
#' scaled.
#'
#' @format A data frame with 2,630 rows and 6 columns: `la` (leaf area),
#'   `ln` (leaf nitrogen content), `ph` (plant height), `sla` (specific
#'   leaf area), `ssd` (specific stem density), `sm` (seed mass).
#' @source \doi{10.6084/m9.figshare.13140146}; TRY database
#'   (https://www.try-db.org). Carmona et al. (2021) *Nature*, 597,
#'   683-687.
"gspff"

#' Global spectrum of plant form and function: incomplete trait data
#'
#' The same six aboveground traits for 10,746 species with incomplete
#' trait information (species with at least three traits measured);
#' log10-transformed and scaled. Useful together with [gspff_phylo] to
#' demonstrate imputation with [fs_impute()].
#'
#' @format A data frame with 10,746 rows and 6 columns; see [gspff].
#' @source \doi{10.6084/m9.figshare.13140146}.
"gspff_missing"

#' Taxonomy for the complete global spectrum species
#'
#' Genus, family and order for the 2,630 species in [gspff].
#'
#' @format A data frame with 2,630 rows and 3 columns: `genus`, `family`,
#'   `order`.
#' @source \doi{10.6084/m9.figshare.13140146}.
"gspff_tax"

#' Phylogeny for the incomplete global spectrum species
#'
#' Phylogenetic tree (class `phylo`, \pkg{ape}) for the 10,746 species in
#' [gspff_missing].
#'
#' @format An object of class `phylo` with 10,746 tips.
#' @source Distributed with the \pkg{funspace} package.
"gspff_phylo"
