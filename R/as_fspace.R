# as_fspace(): import existing ordinations ------------------------------------

#' Convert existing ordinations to an fspace
#'
#' Entry point for users with an existing ordination: objects from
#' [stats::prcomp()], [stats::princomp()], `ade4::dudi.pca()`,
#' `vegan::metaMDS()`, a `dist` object, or a plain coordinates matrix.
#'
#' @param x The object to convert.
#' @param traits Optional units-by-traits table stored alongside the space
#'   (required by some downstream evaluation functions).
#' @param correction Correction for negative eigenvalues when converting a
#'   `dist` object via principal coordinates analysis: `"cailliez"`
#'   (default), `"lingoes"`, or `"none"`.
#' @param ... Passed to methods.
#' @return An object of class `fspace`.
#' @seealso [fs_space()]
#' @examples
#' data(gspff)
#' p <- prcomp(gspff, center = TRUE, scale. = TRUE)
#' sp <- as_fspace(p, traits = gspff)
#' sp
#'
#' # from a distance matrix (PCoA)
#' d <- dist(scale(gspff[1:100, ]))
#' as_fspace(d)
#' @export
as_fspace <- function(x, ...) UseMethod("as_fspace")

#' @rdname as_fspace
#' @export
as_fspace.prcomp <- function(x, traits = NULL, ...) {
  coords <- x$x
  if (is.null(rownames(coords))) {
    stop("The prcomp object has no row names; refit with named data.",
         call. = FALSE)
  }
  V <- x$rotation
  scl <- if (identical(x$scale, FALSE)) {
    stats::setNames(rep(1, length(x$center)), names(x$center))
  } else {
    x$scale
  }
  new_fspace(coords = coords, method = "pca",
             traits = traits, units = .units_from(coords),
             loadings = t(x$sdev * t(V)), eigenvectors = V,
             eig = unname(x$sdev^2), center = x$center,
             scale_values = scl, proj = V,
             scale = !identical(x$scale, FALSE), call = match.call())
}

#' @rdname as_fspace
#' @export
as_fspace.princomp <- function(x, traits = NULL, ...) {
  coords <- x$scores
  colnames(coords) <- paste0("PC", seq_len(ncol(coords)))
  V <- unclass(x$loadings)
  colnames(V) <- colnames(coords)
  loadings <- t(x$sdev * t(V))
  colnames(loadings) <- colnames(coords)
  new_fspace(coords = coords, method = "pca",
             traits = traits, units = .units_from(coords),
             loadings = loadings, eigenvectors = V,
             eig = unname(x$sdev^2), center = x$center,
             scale_values = x$scale, proj = V,
             scale = NA, call = match.call())
}

#' @rdname as_fspace
#' @export
as_fspace.dist <- function(x, correction = c("cailliez", "lingoes", "none"),
                           traits = NULL, ...) {
  correction <- match.arg(correction)
  out <- .space_pcoa(NULL, x, correction)
  new_fspace(coords = out$coords, method = "pcoa",
             traits = traits, units = .units_from(out$coords),
             eig = out$eig, dist = x, call = match.call())
}

#' @rdname as_fspace
#' @export
as_fspace.matrix <- function(x, traits = NULL, ...) {
  if (is.null(rownames(x))) {
    stop("Coordinate matrix must have row names identifying units.",
         call. = FALSE)
  }
  new_fspace(coords = x, method = "raw",
             traits = traits, units = .units_from(x), call = match.call())
}

#' @rdname as_fspace
#' @export
as_fspace.default <- function(x, traits = NULL, ...) {
  cls <- class(x)[1L]
  # vegan::metaMDS
  if (inherits(x, "metaMDS")) {
    coords <- x$points
    colnames(coords) <- paste0("NMDS", seq_len(ncol(coords)))
    sp <- new_fspace(coords = coords, method = "nmds",
                     traits = traits, units = .units_from(coords),
                     stress = x$stress, call = match.call())
    return(sp)
  }
  # ade4::dudi.pca (class 'pca' 'dudi')
  if (inherits(x, "dudi")) {
    coords <- as.matrix(x$li)
    colnames(coords) <- paste0("PC", seq_len(ncol(coords)))
    sp <- new_fspace(coords = coords, method = "pca",
                     traits = traits, units = .units_from(coords),
                     loadings = as.matrix(x$c1), eig = x$eig,
                     call = match.call())
    return(sp)
  }
  stop("Don't know how to convert an object of class '", cls,
       "' to an fspace.", call. = FALSE)
}

.units_from <- function(coords) {
  data.frame(id = rownames(coords), n_obs = 1L, has_own_obs = FALSE,
             imputed_traits = FALSE, stringsAsFactors = FALSE)
}
