# fs_reduce() and fs_rotate(): dimensionality choice after evaluation ---------

#' Reduce a trait space to a chosen dimensionality
#'
#' Trims a full space (PCA, PCoA, raw) to the first `dims` axes, or refits
#' NMDS at `dims` dimensions (NMDS solutions cannot be truncated). Run the
#' evaluation functions (`fs_dimensionality()`, `fs_quality()`,
#' `fs_adequacy()`, forthcoming) before choosing `dims`.
#'
#' @param space An `fspace` object.
#' @param dims Number of dimensions to keep (positive integer).
#' @return The reduced `fspace`.
#' @export
fs_reduce <- function(space, dims) {
  stopifnot(is_fspace(space))
  dims <- as.integer(dims)
  d_now <- ncol(space$coords)
  if (length(dims) != 1L || is.na(dims) || dims < 1L) {
    stop("`dims` must be a single positive integer.", call. = FALSE)
  }
  if (space$method != "nmds" && dims > d_now) {
    stop("`dims` (", dims, ") exceeds the available dimensions (",
         d_now, ").", call. = FALSE)
  }
  if (!is.null(space$tpds)) {
    warning("Reducing the space invalidates previously estimated TPDs; ",
            "they have been removed. Re-run fs_tpd().", call. = FALSE)
    space$tpds <- NULL
    space$bw <- NULL
  }
  if (space$method == "nmds") {
    if (is.null(space$dist)) {
      stop("This NMDS space stores no distance matrix, so it cannot be ",
           "refitted. Rebuild it with fs_space() or import the refit via ",
           "as_fspace().", call. = FALSE)
    }
    fit <- MASS::isoMDS(space$dist, k = dims, trace = FALSE)
    coords <- fit$points
    colnames(coords) <- paste0("NMDS", seq_len(ncol(coords)))
    rownames(coords) <- labels(space$dist)
    space$coords <- coords
    space$stress <- fit$stress / 100
  } else {
    space$dims_full <- d_now
    space$coords <- space$coords[, seq_len(dims), drop = FALSE]
    if (!is.null(space$loadings)) {
      space$loadings <- space$loadings[, seq_len(dims), drop = FALSE]
    }
    if (!is.null(space$eig)) {
      space$eig_full <- space$eig
      space$eig <- space$eig[seq_len(dims)]
    }
  }
  space$reduced <- TRUE
  space
}

#' Rotate a reduced PCA space
#'
#' Applies varimax rotation to the retained axes of a reduced PCA space.
#' Rotation is only meaningful after reduction (rotating the full space and
#' then truncating gives different, incorrect results), and only for PCA
#' (PCoA/NMDS axes have no loadings to rotate).
#'
#' @param space A reduced PCA `fspace`.
#' @param method Rotation method; only `"varimax"` in this version.
#' @return The rotated `fspace`.
#' @export
fs_rotate <- function(space, method = "varimax") {
  stopifnot(is_fspace(space))
  method <- match.arg(method)
  if (space$method != "pca") {
    stop("Rotation applies to PCA spaces only (", space$method,
         " axes have no trait loadings to rotate).", call. = FALSE)
  }
  if (!isTRUE(space$reduced)) {
    stop("Rotate after reduction: run fs_reduce() first. ",
         "(Rotating the full space and then truncating is not equivalent ",
         "to truncating and then rotating, and only the latter is valid.)",
         call. = FALSE)
  }
  if (ncol(space$coords) < 2L) {
    stop("Rotation requires at least two retained dimensions.", call. = FALSE)
  }
  if (!is.null(space$tpds)) {
    warning("Rotating the space invalidates previously estimated TPDs; ",
            "they have been removed. Re-run fs_tpd().", call. = FALSE)
    space$tpds <- NULL
    space$bw <- NULL
  }
  rot <- stats::varimax(space$loadings, normalize = TRUE)
  space$loadings <- unclass(rot$loadings)
  space$coords <- space$coords %*% rot$rotmat
  colnames(space$coords) <- paste0("RC", seq_len(ncol(space$coords)))
  colnames(space$loadings) <- colnames(space$coords)
  # eigenvalues no longer apply axis-wise after rotation
  space$eig_rotated <- colSums(space$loadings^2)
  space$rotation <- "varimax"
  space
}
