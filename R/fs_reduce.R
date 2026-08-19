# fs_reduce() and fs_rotate(): dimensionality choice after evaluation ---------

#' Reduce a trait space to a chosen dimensionality
#'
#' Trims a full space (PCA, PCoA, raw) to the first `dims` axes, or refits
#' NMDS at `dims` dimensions (NMDS solutions cannot be truncated). Run the
#' evaluation functions ([fs_dimensionality()], [fs_quality()],
#' [fs_adequacy()]) before choosing `dims`.
#'
#' @param space An `fspace` object.
#' @param dims Number of dimensions to keep (positive integer).
#' @return The reduced `fspace`.
#' @examples
#' data(gspff)
#' sp <- fs_space(gspff[1:300, ], method = "pca")
#' fs_quality(sp)             # evaluate first...
#' sp2 <- fs_reduce(sp, 2)    # ...then choose the dimensionality
#' sp2
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
    for (slot in c("loadings", "eigenvectors", "proj")) {
      if (!is.null(space[[slot]])) {
        space[[slot]] <- space[[slot]][, seq_len(dims), drop = FALSE]
      }
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
#' (PCoA/NMDS axes have no loadings to rotate). Two conventions are
#' offered:
#'
#' * `type = "rigid"` (default): the varimax criterion is applied to the
#'   eigenvectors and the scores are rotated by the same orthogonal
#'   matrix. The configuration of units is untouched -- all pairwise
#'   distances are exactly preserved -- and the reported loadings are the
#'   correlations between traits and the rotated scores. Note that
#'   rigidly rotated components are in general no longer uncorrelated.
#' * `type = "rescaled"`: the factor-analytic convention (as in
#'   \pkg{psych}): varimax on the scaled loadings, with scores rescaled
#'   so each rotated component's variance equals its sum of squared
#'   loadings. Components remain closer to the factor-analysis reading,
#'   but pairwise distances between units are not preserved.
#'
#' @param space A reduced PCA `fspace`.
#' @param method Rotation method; only `"varimax"` in this version.
#' @param type Rotation convention; see Description.
#' @return The rotated `fspace` (fields `loadings`, `eig_rotated`,
#'   `rotmat` and `rotation_type` describe the rotation).
#' @examples
#' data(gspff)
#' sp2 <- fs_reduce(fs_space(gspff, method = "pca"), 2)
#'
#' rig <- fs_rotate(sp2)                      # rigid (default)
#' rig$loadings                               # trait-axis correlations
#' # geometry untouched:
#' all.equal(dist(rig$coords), dist(sp2$coords), check.attributes = FALSE)
#'
#' res <- fs_rotate(sp2, type = "rescaled")   # factor-analytic
#' colSums(res$loadings^2)                    # variance per rotated axis
#' @export
fs_rotate <- function(space, method = "varimax",
                      type = c("rigid", "rescaled")) {
  stopifnot(is_fspace(space))
  method <- match.arg(method)
  type <- match.arg(type)
  if (space$rotation != "none") {
    stop("This space is already rotated; rebuild and reduce it to rotate ",
         "with a different convention.", call. = FALSE)
  }
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
  k <- ncol(space$coords)
  if (type == "rigid") {
    # varimax criterion on the eigenvectors; scores rotated by the same
    # orthogonal matrix -> the configuration (all pairwise distances) is
    # exactly preserved.
    R <- stats::varimax(space$eigenvectors, normalize = TRUE)$rotmat
    coords <- space$coords %*% R
    colnames(coords) <- paste0("RC", seq_len(k))
    ss <- apply(coords, 2L, function(z) mean((z - mean(z))^2))
    Vr <- space$eigenvectors %*% R
    colnames(Vr) <- colnames(coords)
    # loadings = correlations between traits and the rotated scores
    # (reduces to eigenvectors x sdev when R = identity)
    L_rot <- (space$eigenvectors %*% diag(space$eig, k) %*% R) %*%
      diag(1 / sqrt(ss), k)
    colnames(L_rot) <- colnames(coords)
    space$coords <- coords
    space$eigenvectors <- Vr
    space$loadings <- L_rot
    if (!is.null(space$proj)) {
      space$proj <- Vr
    }
  } else {
    # factor-analytic convention: varimax on the SCALED loadings, scores
    # rescaled so each rotated component's variance equals its sum of
    # squared loadings. Distances between units are NOT preserved.
    L <- space$loadings
    rot <- stats::varimax(L, normalize = TRUE)
    R <- rot$rotmat
    L_rot <- unclass(rot$loadings)
    ss <- colSums(L_rot^2)
    sdev_ret <- sqrt(space$eig)
    S_std <- sweep(space$coords, 2L, sdev_ret, "/")
    coords <- S_std %*% R %*% diag(sqrt(ss), k)
    colnames(coords) <- paste0("RC", seq_len(k))
    colnames(L_rot) <- colnames(coords)
    space$coords <- coords
    space$loadings <- L_rot
    if (!is.null(space$proj)) {
      space$proj <- space$eigenvectors %*% diag(1 / sdev_ret, k) %*%
        R %*% diag(sqrt(ss), k)
      colnames(space$proj) <- colnames(coords)
    }
  }
  space$eig_rotated <- ss
  space$rotmat <- R
  space$rotation <- "varimax"
  space$rotation_type <- type
  space
}
