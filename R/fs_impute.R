# fs_impute(): trait imputation ------------------------------------------------

#' Impute missing trait values
#'
#' Random-forest imputation of missing trait values via
#' \pkg{missForest} (Stekhoven & Buhlmann 2012), optionally informed by
#' phylogeny through phylogenetic eigenvectors (the funspace approach).
#' Which values were imputed is recorded in the result and propagated by
#' [fs_space()] into the space's unit metadata, so imputation is never
#' invisible downstream.
#'
#' This is a convenience wrapper, not an imputation research tool: for
#' serious gap-filling exercises, evaluate imputation quality for your
#' data and traits.
#'
#' @param traits Data.frame (units x traits) with unique row names;
#'   numeric and factor columns are supported. Missing values (`NA`) are
#'   imputed.
#' @param phylo Optional phylogeny (class `phylo`, \pkg{ape}) whose tip
#'   labels include all row names of `traits`. Its cophenetic distances
#'   are converted to eigenvectors used as additional predictors (and
#'   dropped from the output).
#' @param n_eigen Number of phylogenetic eigenvectors (default 10, capped
#'   at units - 1).
#' @param seed Random seed.
#' @param ... Passed to [missForest::missForest()].
#'
#' @return The completed trait data.frame, with attribute `imputed`: a
#'   list with `cells` (a units x traits logical matrix marking imputed
#'   values), `units` (row names with at least one imputed value),
#'   `complete` (row names that were fully empirical, useful to build a
#'   space from empirical data only and project the rest with
#'   [fs_project()]) and `oob` (missForest out-of-bag error estimate).
#' @references Stekhoven, D.J. & Buhlmann, P. (2012) MissForest:
#'   non-parametric missing value imputation for mixed-type data.
#'   *Bioinformatics*, 28, 112-118.
#' @seealso [fs_space()]
#' @examples
#' \donttest{
#' data(gspff_missing)
#' data(gspff_phylo)
#' set.seed(1)
#' sub <- gspff_missing[sample(nrow(gspff_missing), 150), ]
#'
#' out <- fs_impute(sub, seed = 1)
#' attr(out, "imputed")$oob            # out-of-bag error estimate
#'
#' # informed by the phylogeny:
#' outp <- fs_impute(sub, phylo = gspff_phylo, n_eigen = 5, seed = 1)
#' head(attr(outp, "imputed")$units)
#'
#' # recommended workflow: build the space from EMPIRICAL data only and
#' # project the imputed species into it (same centring and scaling):
#' data(gspff)
#' sp <- fs_reduce(fs_space(gspff, method = "pca"), 2)
#' proj <- fs_project(sp, out)
#' head(proj)
#' }
#' @export
fs_impute <- function(traits, phylo = NULL, n_eigen = 10L, seed = NULL,
                      ...) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("fs_impute() needs the 'missForest' package: ",
         "install.packages(\"missForest\").", call. = FALSE)
  }
  traits <- as.data.frame(traits)
  if (is.null(rownames(traits)) ||
      anyDuplicated(rownames(traits)) > 0L ||
      identical(rownames(traits), as.character(seq_len(nrow(traits))))) {
    stop("`traits` must have unique, informative row names identifying ",
         "the units.", call. = FALSE)
  }
  na_mask <- is.na(traits)
  if (!any(na_mask)) {
    message("No missing values found; returning `traits` unchanged.")
    attr(traits, "imputed") <- list(cells = na_mask,
                                    units = character(0),
                                    complete = rownames(traits),
                                    oob = NULL)
    return(traits)
  }
  ximp_input <- traits
  eig_cols <- character(0)
  if (!is.null(phylo)) {
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("Phylogenetic imputation needs the 'ape' package.",
           call. = FALSE)
    }
    if (!inherits(phylo, "phylo")) {
      stop("`phylo` must be a phylo object.", call. = FALSE)
    }
    missing_tips <- setdiff(rownames(traits), phylo$tip.label)
    if (length(missing_tips)) {
      stop("Units absent from the phylogeny: ",
           paste(utils::head(missing_tips, 5L), collapse = ", "),
           if (length(missing_tips) > 5L) ", ...", call. = FALSE)
    }
    ph <- ape::keep.tip(phylo, rownames(traits))
    pd <- ape::cophenetic.phylo(ph)
    k <- min(as.integer(n_eigen), nrow(traits) - 1L)
    pe <- suppressWarnings(
      stats::cmdscale(stats::as.dist(pd), k = k)
    )
    pe <- pe[rownames(traits), , drop = FALSE]
    eig_cols <- paste0(".phylo", seq_len(ncol(pe)))
    colnames(pe) <- eig_cols
    ximp_input <- cbind(traits, as.data.frame(pe))
  }
  if (!is.null(seed)) set.seed(seed)
  fit <- missForest::missForest(ximp_input, ...)
  out <- fit$ximp[, setdiff(colnames(fit$ximp), eig_cols), drop = FALSE]
  rownames(out) <- rownames(traits)
  attr(out, "imputed") <- list(
    cells = na_mask,
    units = rownames(traits)[rowSums(na_mask) > 0L],
    complete = rownames(traits)[rowSums(na_mask) == 0L],
    oob = fit$OOBerror
  )
  out
}
