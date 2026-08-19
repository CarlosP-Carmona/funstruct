# fs_angles(): pairwise angles between traits in the retained space ----------

#' Angles between traits in the retained space
#'
#' Computes the pairwise angles that traits form in the retained
#' dimensions of a space, from their loading vectors across the axes.
#' Bivariate trait plots show relationships two traits at a time; the
#' angles summarize how the whole trait battery is organized in the
#' multidimensional space (Bueno et al. 2023; Beccari & Carmona 2024):
#' 0 degrees means two traits point in the same direction (perfectly
#' positively associated as represented in the space), 90 degrees means
#' they are orthogonal (unrelated), and 180 degrees means they are
#' opposed.
#'
#' The trait vectors are the rows of the space's loadings over the
#' retained axes. For a PCA space these are the construction loadings
#' (eigenvectors scaled by the component standard deviations); for PCoA
#' and NMDS spaces run [fs_loadings()] first, and the post-hoc
#' correlation loadings are used (numeric traits only). In the *full*
#' PCA space the angle between two traits is exactly
#' `acos(cor(trait_i, trait_j))` -- the classical geometric reading of a
#' correlation; in a *reduced* space the angles describe the trait
#' relationships as the retained axes represent them, which is the
#' interesting object: run [fs_reduce()] first.
#'
#' By default the angle is the direct vector angle,
#' `acos(sum(l_i * l_j) / (|l_i| |l_j|))` (the uncentered cosine, which
#' gives the exact correspondence with trait correlations above). Set
#' `center = TRUE` to use the Pearson correlation between the two
#' loading vectors instead (centering them across axes first), as done
#' in some published implementations.
#'
#' With `profile = TRUE` (on the full space, before reducing) the
#' angles double as a dimensionality diagnostic: the angle matrix is
#' recomputed at every dimensionality from 2 up to the full space, and
#' each is compared with the full-space angles -- which are exactly the
#' raw trait correlations -- through their correlation over trait pairs
#' and the mean absolute angular deviation (mad, in the units of the
#' angles). A high correlation at `k` dimensions means the reduced
#' space represents the correlation structure among the original traits
#' faithfully. This is the trait-side complement of [fs_quality()],
#' which asks the same question about the distances among *units*; the
#' two need not agree, and a good reduced space should pass both. Note
#' the correlation runs over the `p (p - 1) / 2` trait pairs, so it is
#' a coarse statistic for small trait batteries -- read it together
#' with mad.
#'
#' @param space An `fspace` with loadings: a PCA space, or a PCoA/NMDS
#'   space after [fs_loadings()]. Typically reduced with [fs_reduce()].
#' @param dims Axes to use (integer indices); default all axes present
#'   in the loadings.
#' @param degrees Logical; return degrees (default) or radians.
#' @param center Logical; `FALSE` (default) for the direct vector angle,
#'   `TRUE` for the angle of the centered (Pearson) correlation between
#'   loading vectors.
#' @param profile Logical; if `TRUE`, also compute the angle matrices at
#'   every dimensionality from 2 to the full space and their agreement
#'   with the full-space angles (see Details). Cannot be combined with
#'   `dims`.
#'
#' @return A list of class `fs_angles` with `angles` (traits x traits
#'   symmetric matrix), `cosines`, `dims`, `degrees` and `center`. With
#'   `profile = TRUE`, also `by_dims` (list of angle matrices, one per
#'   dimensionality) and `profile` (data.frame with `dims`, `cor`,
#'   `mad`). `print()` and `plot()` methods are provided
#'   (`plot(x, which = "profile")` draws the preservation curve).
#' @references Bueno, C.G., Toussaint, A., Träger, S., et al. (2023)
#'   Reply to: The importance of trait selection in ecology. *Nature*,
#'   618, E31-E34. \doi{10.1038/s41586-023-06149-7}
#'
#'   Beccari, E. & Carmona, C.P. (2024) Aboveground and belowground
#'   sizes are aligned in the unified spectrum of plant form and
#'   function. *Nature Communications*, 15.
#'   \doi{10.1038/s41467-024-53180-x}
#' @seealso [fs_loadings()], [fs_reduce()], [fs_trait_dim()]
#' @examples
#' data(gspff)
#' sp2 <- fs_reduce(fs_space(gspff, method = "pca"), 2)
#' an <- fs_angles(sp2)
#' an
#' plot(an)
#'
#' # angle preservation across dimensionalities (run on the FULL space):
#' ap <- fs_angles(fs_space(gspff, method = "pca"), profile = TRUE)
#' ap$profile
#' plot(ap, which = "profile")
#'
#' # PCoA route: post-hoc loadings first
#' spo <- fs_loadings(fs_space(fs_dist(gspff[1:150, ]), "pcoa"),
#'                    gspff[1:150, ])
#' fs_angles(fs_reduce(spo, 2))
#' @export
fs_angles <- function(space, dims = NULL, degrees = TRUE,
                      center = FALSE, profile = FALSE) {
  stopifnot(is_fspace(space))
  L <- space$loadings
  if (is.null(L) || nrow(L) < 2L) {
    stop("This space carries no trait loadings. PCA spaces have them ",
         "from construction; for PCoA/NMDS run fs_loadings() first.",
         call. = FALSE)
  }
  if (profile && !is.null(dims)) {
    stop("`dims` and `profile = TRUE` cannot be combined: the profile ",
         "scans all dimensionalities.", call. = FALSE)
  }
  if (is.null(dims)) dims <- seq_len(ncol(L))
  dims <- as.integer(dims)
  if (any(is.na(dims)) || any(dims < 1L) || any(dims > ncol(L))) {
    stop("`dims` must be axis indices between 1 and ", ncol(L), ".",
         call. = FALSE)
  }
  if (length(dims) < 2L) {
    stop("At least two axes are needed to compute angles.", call. = FALSE)
  }
  M <- unclass(L)
  res <- .angle_mat(M[, dims, drop = FALSE], center, degrees)
  out <- list(angles = res$angles, cosines = res$cosines, dims = dims,
              degrees = degrees, center = center,
              space_method = space$method, call = match.call())
  if (profile) {
    K <- ncol(M)
    ut <- upper.tri(res$angles)
    by_dims <- vector("list", K - 1L)
    names(by_dims) <- paste0("d", 2:K)
    prof <- data.frame(dims = 2:K, cor = NA_real_, mad = NA_real_)
    for (k in 2:K) {
      Ak <- .angle_mat(M[, seq_len(k), drop = FALSE], center,
                       degrees)$angles
      by_dims[[paste0("d", k)]] <- Ak
      prof$cor[prof$dims == k] <- stats::cor(Ak[ut], res$angles[ut])
      prof$mad[prof$dims == k] <- mean(abs(Ak[ut] - res$angles[ut]))
    }
    out$by_dims <- by_dims
    out$profile <- prof
  }
  class(out) <- "fs_angles"
  out
}

# angle matrix from trait-loading rows over a set of axes ---------------------

.angle_mat <- function(M, center, degrees) {
  if (center) {
    cosines <- stats::cor(t(M))
  } else {
    nrm <- sqrt(rowSums(M^2))
    if (any(nrm <= 0)) {
      stop("Trait(s) with zero loading vector on the selected axes: ",
           paste(rownames(M)[nrm <= 0], collapse = ", "), call. = FALSE)
    }
    cosines <- (M %*% t(M)) / outer(nrm, nrm)
  }
  cosines <- pmin(pmax(cosines, -1), 1)
  ang <- acos(cosines)
  if (degrees) ang <- ang * 180 / pi
  diag(ang) <- 0
  dimnames(ang) <- dimnames(cosines) <- list(rownames(M), rownames(M))
  list(angles = ang, cosines = cosines)
}

#' @export
print.fs_angles <- function(x, ...) {
  cat("<fs_angles> pairwise trait angles over axes ",
      paste(x$dims, collapse = ", "), " (", x$space_method, " space)\n",
      sep = "")
  print(round(x$angles, 1))
  cat(if (x$degrees) "0 = aligned, 90 = unrelated, 180 = opposed ",
      "(degrees).\n", sep = "")
  if (!is.null(x$profile)) {
    cat("\nAngle preservation across dimensionalities ",
        "(reference: all ", max(x$profile$dims), " axes):\n", sep = "")
    pr <- x$profile
    pr$cor <- round(pr$cor, 3)
    pr$mad <- round(pr$mad, 2)
    print.data.frame(pr, row.names = FALSE)
    cat("cor: correlation with the full-space angles; mad: mean ",
        "absolute angular deviation",
        if (x$degrees) " (degrees)", ".\n", sep = "")
  }
  invisible(x)
}

#' @export
plot.fs_angles <- function(x, which = c("angles", "profile"), ...) {
  which <- match.arg(which)
  if (which == "profile") {
    if (is.null(x$profile)) {
      stop("No profile stored; re-run fs_angles() with profile = TRUE.",
           call. = FALSE)
    }
    pr <- x$profile
    graphics::plot(pr$dims, pr$cor, type = "b", pch = 16,
                   ylim = range(pr$cor, 1), las = 1,
                   xlab = "Dimensions retained",
                   ylab = "Correlation with full-space angles", ...)
    graphics::abline(h = 1, lty = 3, col = "grey60")
    return(invisible(pr))
  }
  A <- x$angles
  k <- nrow(A)
  full <- if (x$degrees) 180 else pi
  pal <- grDevices::hcl.colors(101, "Blue-Red 2")
  z <- t(A[k:1, , drop = FALSE])   # z[x, y] drawn at (x, y)
  graphics::image(seq_len(k), seq_len(k), z,
                  zlim = c(0, full), col = pal, axes = FALSE,
                  xlab = "", ylab = "", ...)
  graphics::axis(1L, at = seq_len(k), labels = colnames(A), las = 2,
                 tick = FALSE)
  graphics::axis(2L, at = seq_len(k), labels = rev(rownames(A)), las = 1,
                 tick = FALSE)
  graphics::text(rep(seq_len(k), times = k), rep(seq_len(k), each = k),
                 labels = round(as.vector(z), 0), cex = 0.8)
  graphics::box()
  invisible(A)
}
