# fs_space(): build the full trait space -------------------------------------

#' Build a trait space
#'
#' Builds the full trait space from a units-by-traits table. Dimensionality
#' is deliberately not chosen here: evaluate the space first (see
#' [fs_dimensionality()] and [fs_quality()]) and then trim it with
#' [fs_reduce()] and, for PCA spaces, optionally rotate it with
#' [fs_rotate()].
#'
#' @param traits A matrix or data.frame (units x traits) with unique row
#'   names identifying the units (species, populations, individuals). All
#'   columns must be numeric for `method = "pca"`, `"nmds"` and `"raw"`;
#'   mixed types are allowed for `method = "pcoa"` (Gower dissimilarity is
#'   used internally via [cluster::daisy()]).
#' @param method Ordination method: `"pca"` (default), `"pcoa"`, `"nmds"`,
#'   or `"raw"` (scaled trait axes used directly).
#' @param scale Logical; scale traits to unit variance (default `TRUE`).
#'   Ignored for `method = "pcoa"` with mixed data (Gower scales internally).
#' @param dist Optional distance object (`dist`, e.g. from [stats::dist()] or
#'   gawdis) used for `"pcoa"`/`"nmds"`; overrides internal computation.
#' @param correction Correction for negative eigenvalues in PCoA:
#'   `"cailliez"` (default), `"lingoes"`, or `"none"`.
#' @param k Number of dimensions for NMDS only (NMDS has no full space; it is
#'   fitted at a chosen dimensionality). Default `min(4, ncol(traits) - 1)`.
#'   [fs_reduce()] refits NMDS at the requested dimensionality.
#' @param seed Random seed used for methods with a stochastic component.
#'
#' @return An object of class `fspace`.
#' @seealso [as_fspace()] to import existing ordinations, [fs_reduce()],
#'   [fs_rotate()].
#' @export
fs_space <- function(traits,
                     method = c("pca", "pcoa", "nmds", "raw"),
                     scale = TRUE,
                     dist = NULL,
                     correction = c("cailliez", "lingoes", "none"),
                     k = NULL,
                     seed = NULL) {
  method <- match.arg(method)
  correction <- match.arg(correction)
  traits <- .check_traits(traits, method)
  if (!is.null(seed)) set.seed(seed)

  units <- data.frame(
    id = rownames(traits),
    n_obs = 1L,
    has_own_obs = FALSE,
    imputed_traits = .imputed_flags(traits),
    stringsAsFactors = FALSE
  )

  out <- switch(method,
    pca  = .space_pca(traits, scale),
    pcoa = .space_pcoa(traits, dist, correction),
    nmds = .space_nmds(traits, dist, k),
    raw  = .space_raw(traits, scale)
  )

  new_fspace(
    coords = out$coords, method = method, traits = traits, units = units,
    loadings = out$loadings, eig = out$eig, stress = out$stress,
    dist = out$dist, scale = scale, call = match.call()
  )
}

# internal builders -----------------------------------------------------------

.space_pca <- function(x, scale) {
  m <- as.matrix(x)
  p <- stats::prcomp(m, center = TRUE, scale. = scale)
  coords <- p$x
  colnames(coords) <- paste0("PC", seq_len(ncol(coords)))
  loadings <- p$rotation
  colnames(loadings) <- colnames(coords)
  list(coords = coords, loadings = loadings, eig = p$sdev^2,
       stress = NULL, dist = NULL)
}

.space_pcoa <- function(x, d, correction) {
  if (is.null(d)) d <- .default_dist(x)
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

.space_nmds <- function(x, d, k) {
  if (is.null(d)) d <- .default_dist(x)
  if (is.null(k)) k <- max(1L, min(4L, ncol(x) - 1L))
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

.check_traits <- function(traits, method) {
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
  if (method != "pcoa" && !all(numeric_ok)) {
    stop("Non-numeric traits found: ",
         paste(names(traits)[!numeric_ok], collapse = ", "),
         ". Use method = \"pcoa\" (Gower) for mixed trait types.",
         call. = FALSE)
  }
  if (method != "pcoa" && anyNA(traits)) {
    stop("Missing trait values found. Impute first (see fs_impute()) or ",
         "use method = \"pcoa\", which tolerates missing values via Gower.",
         call. = FALSE)
  }
  traits
}

.default_dist <- function(x) {
  numeric_ok <- vapply(as.data.frame(x), is.numeric, logical(1L))
  if (all(numeric_ok) && !anyNA(x)) {
    stats::dist(base::scale(as.matrix(x)))
  } else {
    cluster::daisy(as.data.frame(x), metric = "gower")
  }
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
  rownames(traits) %in% imp$units
}
