# fs_space(): build the full trait space -------------------------------------

#' Build a trait space
#'
#' Builds the full trait space, either from a units-by-traits table
#' (`method = "pca"` or `"raw"`) or from a precomputed dissimilarity
#' matrix (`method = "pcoa"` or `"nmds"`). Dimensionality is deliberately
#' not chosen here: evaluate the space first (see [fs_dimensionality()]
#' and [fs_quality()]) and then trim it with [fs_reduce()] and, for PCA
#' spaces, optionally rotate it with [fs_rotate()].
#'
#' `fs_space()` never computes dissimilarities itself. For `"pcoa"` and
#' `"nmds"`, pass a [stats::dist] object as `traits` -- typically the
#' output of [fs_dist()], which combines multiple trait types into a
#' Gower-style dissimilarity bounded between 0 and 1 (any `dist` object,
#' e.g. from [stats::dist()] or the gawdis package, is accepted).
#' Keeping the dissimilarity step explicit makes the choice of
#' coefficient a visible part of the analysis instead of a hidden
#' default.
#'
#' @section Negative eigenvalues and their corrections (PCoA):
#'
#' Metric scaling of a non-Euclidean dissimilarity matrix (Gower
#' dissimilarities usually are non-Euclidean, especially with categorical
#' traits or missing values) produces negative eigenvalues: the units
#' cannot be embedded exactly in any Euclidean space, and the axes
#' associated with negative eigenvalues represent "imaginary" variation.
#' Left uncorrected (`correction = "none"`), distances in the resulting
#' space underestimate the original dissimilarities, and quantities
#' computed from the axes (functional richness, dispersion) absorb the
#' distortion. Two classical fixes are offered:
#'
#' * `"cailliez"` (default): adds the smallest constant to all
#'   off-diagonal *dissimilarities* that makes the matrix Euclidean
#'   (Cailliez 1983). This is the correction used by
#'   `ape::pcoa(..., correction = "cailliez")` and is a common default in
#'   trait-based work.
#' * `"lingoes"`: adds the smallest constant to all off-diagonal
#'   *squared* dissimilarities that makes the matrix Euclidean (Lingoes
#'   1971). It distorts large dissimilarities relatively less than
#'   Cailliez; the two usually give very similar spaces.
#'
#' Either correction alters all dissimilarities slightly, so the
#' corrected space no longer reproduces the input distances exactly;
#' [fs_quality()] measures how much is lost. If the input is already
#' Euclidean (e.g. `stats::dist()` on scaled numeric traits), no
#' correction is applied in practice because the added constant is zero
#' (Cailliez) or the eigenvalues are already non-negative (Lingoes).
#' Axes with eigenvalues that remain non-positive after correction are
#' dropped.
#'
#' @param traits For `method = "pca"` and `"raw"`: a matrix or data.frame
#'   (units x traits) with unique row names and all-numeric columns
#'   (impute missing values first, see [fs_impute()]). For
#'   `method = "pcoa"` and `"nmds"`: a [stats::dist] object, typically
#'   from [fs_dist()]. Passing a trait table for `"pcoa"`/`"nmds"` is an
#'   error -- compute the dissimilarities first with [fs_dist()].
#' @param method Ordination method: `"pca"` (default; via
#'   [stats::princomp()], which requires more units than traits),
#'   `"pcoa"`, `"nmds"`, or `"raw"` (scaled trait axes used directly).
#'   For PCA spaces the returned object stores both `$loadings` (the
#'   eigenvectors scaled by the component standard deviations, i.e. the
#'   trait-axis correlations when `scale = TRUE` -- loadings in the
#'   factor-analytic sense) and `$eigenvectors` (the unit-norm
#'   eigenvectors used internally for projection).
#' @param scale Logical; scale traits to unit variance (default `TRUE`).
#'   Used by `"pca"` (correlation vs covariance matrix) and `"raw"` only.
#'   With `scale = FALSE` the PCA is computed on the covariance matrix,
#'   so traits with larger variances dominate the leading axes; a warning
#'   is issued when the trait variances differ noticeably.
#' @param correction Correction for negative eigenvalues in PCoA:
#'   `"cailliez"` (default), `"lingoes"`, or `"none"`. See the section
#'   below.
#' @param k Number of dimensions for NMDS only (NMDS has no full space; it
#'   is fitted at a chosen dimensionality). Default `2`, the usual NMDS
#'   convention; [fs_reduce()] refits NMDS at any other dimensionality.
#' @param seed Random seed used for methods with a stochastic component.
#'
#' @return An object of class `fspace`.
#' @references Cailliez, F. (1983) The analytical solution of the
#'   additive constant problem. *Psychometrika*, 48, 305-308.
#'
#'   Lingoes, J.C. (1971) Some boundary conditions for a monotone
#'   analysis of symmetric matrices. *Psychometrika*, 36, 195-203.
#' @seealso [fs_dist()] to compute bounded dissimilarities,
#'   [as_fspace()] to import existing ordinations, [fs_reduce()],
#'   [fs_rotate()].
#' @examples
#' data(gspff)
#' sp <- fs_space(gspff, method = "pca")
#' sp
#' plot(sp)
#'
#' # PCoA: compute the dissimilarities explicitly, then ordinate
#' d <- fs_dist(gspff[1:150, ])
#' sp2 <- fs_space(d, method = "pcoa")
#' @export
fs_space <- function(traits,
                     method = c("pca", "pcoa", "nmds", "raw"),
                     scale = TRUE,
                     correction = c("cailliez", "lingoes", "none"),
                     k = NULL,
                     seed = NULL) {
  method <- match.arg(method)
  correction <- match.arg(correction)
  is_dist <- inherits(traits, "dist")

  if (method %in% c("pcoa", "nmds")) {
    if (!is_dist) {
      stop("`", method, "` requires a precomputed dissimilarity matrix ",
           "(a `dist` object) as the `traits` input. fs_space() does not ",
           "compute dissimilarities; use fs_dist() first.", call. = FALSE)
    }
    d <- traits
    if (is.null(labels(d))) {
      stop("The `dist` object must carry unit labels ",
           "(see ?dist, argument `Labels`).", call. = FALSE)
    }
    if (anyNA(d)) {
      stop("The dissimilarity matrix contains NAs.", call. = FALSE)
    }
    traits <- NULL
  } else {
    if (is_dist) {
      stop("A `dist` object can only be used with method = \"pcoa\" or ",
           "\"nmds\"; `", method, "` needs the trait values themselves.",
           call. = FALSE)
    }
    traits <- .check_traits(traits)
  }
  if (!is.null(seed)) set.seed(seed)

  unit_ids <- if (is.null(traits)) labels(d) else rownames(traits)
  units <- data.frame(
    id = unit_ids,
    n_obs = 1L,
    has_own_obs = FALSE,
    imputed_traits = if (is.null(traits)) FALSE else .imputed_flags(traits),
    stringsAsFactors = FALSE
  )

  out <- switch(method,
    pca  = .space_pca(traits, scale),
    pcoa = .space_pcoa(d, correction),
    nmds = .space_nmds(d, k),
    raw  = .space_raw(traits, scale)
  )

  new_fspace(
    coords = out$coords, method = method, traits = traits, units = units,
    loadings = out$loadings, eigenvectors = out$eigenvectors,
    eig = out$eig, center = out$center,
    scale_values = out$scale_values, proj = out$proj,
    stress = out$stress, dist = out$dist, scale = scale,
    call = match.call()
  )
}

# internal builders -----------------------------------------------------------

.space_pca <- function(x, scale) {
  m <- as.matrix(x)
  if (nrow(m) <= ncol(m)) {
    stop("PCA (princomp) requires more units than traits (got ",
         nrow(m), " units x ", ncol(m), " traits).", call. = FALSE)
  }
  if (!scale) {
    v <- apply(m, 2L, stats::var)
    if (max(v) / min(v) > 1.5) {
      warning("PCA on the covariance matrix (scale = FALSE) with unequal ",
              "trait variances: axes will be dominated by the ",
              "high-variance traits (variances range from ",
              signif(min(v), 3), " to ", signif(max(v), 3),
              "). Use scale = TRUE unless the traits share a scale ",
              "deliberately.", call. = FALSE)
    }
  }
  p <- stats::princomp(m, cor = scale)
  V <- unclass(p$loadings)          # unit eigenvectors
  coords <- p$scores
  colnames(coords) <- paste0("PC", seq_len(ncol(coords)))
  colnames(V) <- colnames(coords)
  # loadings sensu factor analysis: eigenvectors scaled by the component
  # standard deviations (= trait-axis correlations when cor = TRUE);
  # same convention as funspace: t(sdev * t(loadings))
  loadings <- t(p$sdev * t(V))
  colnames(loadings) <- colnames(coords)
  list(coords = coords, loadings = loadings, eigenvectors = V,
       eig = unname(p$sdev^2), center = p$center,
       scale_values = p$scale, proj = V, stress = NULL, dist = NULL)
}

.space_pcoa <- function(d, correction) {
  n <- attr(d, "Size")
  k <- n - 1L
  # cmdscale warns when k exceeds the positive-eigenvalue rank; we handle
  # rank explicitly below, so the warning is suppressed.
  res <- suppressWarnings(switch(correction,
    cailliez = stats::cmdscale(d, k = k, eig = TRUE, add = TRUE),
    lingoes  = stats::cmdscale(.lingoes(d), k = k, eig = TRUE, add = FALSE),
    none     = stats::cmdscale(d, k = k, eig = TRUE, add = FALSE)
  ))
  keep <- which(res$eig > sqrt(.Machine$double.eps))
  keep <- keep[keep <= ncol(res$points)]
  coords <- res$points[, keep, drop = FALSE]
  colnames(coords) <- paste0("PCo", seq_len(ncol(coords)))
  rownames(coords) <- labels(d)
  list(coords = coords, loadings = NULL, eig = res$eig[keep],
       stress = NULL, dist = d)
}

.space_nmds <- function(d, k) {
  if (is.null(k)) k <- 2L
  fit <- MASS::isoMDS(d, k = k, trace = FALSE)
  coords <- fit$points
  colnames(coords) <- paste0("NMDS", seq_len(ncol(coords)))
  rownames(coords) <- labels(d)
  list(coords = coords, loadings = NULL, eig = NULL,
       stress = fit$stress / 100, dist = d)
}

.space_raw <- function(x, scale) {
  m <- as.matrix(x)
  coords <- base::scale(m, center = TRUE, scale = scale)
  attr(coords, "scaled:center") <- NULL
  attr(coords, "scaled:scale") <- NULL
  list(coords = coords, loadings = NULL, eig = NULL,
       stress = NULL, dist = NULL)
}

# helpers ---------------------------------------------------------------------

.check_traits <- function(traits) {
  if (is.null(dim(traits))) {
    stop("`traits` must be a matrix or data.frame (units x traits).",
         call. = FALSE)
  }
  traits <- as.data.frame(traits)
  if (is.null(rownames(traits)) ||
      anyDuplicated(rownames(traits)) > 0L ||
      identical(rownames(traits), as.character(seq_len(nrow(traits))))) {
    stop("`traits` must have unique, informative row names ",
         "identifying the units.", call. = FALSE)
  }
  if (ncol(traits) < 2L) {
    stop("At least two traits are required to build a space.", call. = FALSE)
  }
  numeric_ok <- vapply(traits, is.numeric, logical(1L))
  if (!all(numeric_ok)) {
    stop("Non-numeric traits found: ",
         paste(names(traits)[!numeric_ok], collapse = ", "),
         ". For mixed trait types, compute Gower dissimilarities with ",
         "fs_dist() and use method = \"pcoa\".", call. = FALSE)
  }
  if (anyNA(traits)) {
    stop("Missing trait values found. Impute first (see fs_impute()), or ",
         "compute dissimilarities with fs_dist() (which tolerates missing ",
         "values) and use method = \"pcoa\".", call. = FALSE)
  }
  traits
}

.lingoes <- function(d) {
  # Lingoes correction: add constant c1 to squared dissimilarities
  n <- attr(d, "Size")
  D <- as.matrix(d)
  G <- -0.5 * .center_mat(n) %*% (D^2) %*% .center_mat(n)
  e <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  c1 <- max(0, -min(e))
  if (c1 <= sqrt(.Machine$double.eps)) return(d)
  Dc <- sqrt(D^2 + 2 * c1)
  diag(Dc) <- 0
  stats::as.dist(Dc)
}

.center_mat <- function(n) diag(n) - matrix(1 / n, n, n)

.imputed_flags <- function(traits) {
  imp <- attr(traits, "imputed")
  if (is.null(imp)) return(rep(FALSE, nrow(traits)))
  rownames(traits) %in% imp$imputed
}
