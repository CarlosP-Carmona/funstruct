# fs_dimensionality(): how many dimensions does the space need? ---------------

#' Evaluate the dimensionality of a trait space
#'
#' Suggests how many dimensions are needed to describe a trait space,
#' before reducing it with [fs_reduce()]. Five criteria are available:
#'
#' * `"auc"` (default): elbow of the quality curve from [fs_quality()]
#'   (mSD against dimensionality), in the spirit of the quality-based
#'   assessment of Mouillot et al. (2021). Works for any space with
#'   recoverable initial dissimilarities.
#' * `"elbow"`: elbow of the eigenvalue scree (PCA/PCoA).
#' * `"parallel"`: Horn's parallel analysis (PCA only): keep axes whose
#'   eigenvalue exceeds the 95th percentile of eigenvalues obtained after
#'   permuting each trait independently.
#' * `"end"`: effective number of dimensions (PCA only): the Hill number
#'   of order `q` of the eigenvalue shares (Beccari & Carmona 2024). With
#'   the default `q = 2` this is the inverse Simpson index of the
#'   eigenvalues: a continuous measure of dimensionality, equal to the
#'   number of traits when all axes carry equal variance and approaching
#'   1 as one axis dominates. The suggested (integer) dimensionality is
#'   the rounded END; the continuous value is in `$end`. See
#'   [fs_trait_dim()] for how individual traits contribute to it.
#' * `"stress"`: NMDS stress against dimensionality (refits the NMDS at
#'   each candidate `k`); suggests the smallest `k` with stress <= 0.10
#'   (a common rule of thumb).
#'
#' @param space An `fspace` object.
#' @param method Criterion; see Details.
#' @param n_perm Number of permutations for `method = "parallel"`.
#' @param k_max Largest dimensionality scanned by `method = "stress"`
#'   (default `min(4, n - 2)`).
#' @param q Hill order for `method = "end"` (default 2, the inverse
#'   Simpson convention of Beccari & Carmona 2024; `q = 1` gives the
#'   exponential of the Shannon entropy of the eigenvalue shares).
#' @param seed Random seed for the stochastic criteria.
#'
#' @return A list of class `fs_dimensionality` with elements `method`,
#'   `suggested` and `curve` (a data.frame; its columns depend on the
#'   method); `method = "end"` adds `end` (the continuous value) and `q`.
#'   `print()` and `plot()` methods are provided.
#' @references Mouillot, D. et al. (2021) The dimensionality and structure
#'   of species trait spaces. *Ecology Letters*, 24, 1988-2009.
#'
#'   Beccari, E. & Carmona, C.P. (2024) Aboveground and belowground sizes
#'   are aligned in the unified spectrum of plant form and function.
#'   *Nature Communications*, 15. \doi{10.1038/s41467-024-53180-x}
#' @seealso [fs_quality()], [fs_adequacy()], [fs_reduce()],
#'   [fs_trait_dim()]
#' @examples
#' data(gspff)
#' sp <- fs_space(gspff[1:300, ], method = "pca")
#' fs_dimensionality(sp, method = "elbow")
#' dpar <- fs_dimensionality(sp, method = "parallel", n_perm = 49,
#'                           seed = 1)
#' dpar
#' plot(dpar)
#' @export
fs_dimensionality <- function(space,
                              method = c("auc", "elbow", "parallel",
                                         "end", "stress"),
                              n_perm = 199L, k_max = NULL, q = 2,
                              seed = NULL) {
  stopifnot(is_fspace(space))
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  out <- switch(method,
    auc      = .dim_auc(space),
    elbow    = .dim_elbow(space),
    parallel = .dim_parallel(space, n_perm),
    end      = .dim_end(space, q),
    stress   = .dim_stress(space, k_max)
  )
  out$method <- method
  out$call <- match.call()
  class(out) <- "fs_dimensionality"
  out
}

#' @export
print.fs_dimensionality <- function(x, ...) {
  cat("<fs_dimensionality> method:", x$method, "\n")
  if (!is.null(x$end)) {
    cat("Effective number of dimensions (q = ", x$q, "): ",
        round(x$end, 2), "\n", sep = "")
  }
  if (is.na(x$suggested)) {
    cat("No dimensionality satisfied the criterion; inspect the curve.\n")
  } else {
    cat("Suggested dimensionality:", x$suggested, "\n")
  }
  print.data.frame(x$curve, row.names = FALSE, digits = 4)
  invisible(x)
}

#' @export
plot.fs_dimensionality <- function(x, ...) {
  cv <- x$curve
  ylab <- setdiff(colnames(cv), c("dims", "axis"))[1L]
  xcol <- intersect(c("dims", "axis"), colnames(cv))[1L]
  graphics::plot(cv[[xcol]], cv[[ylab]], type = "b", pch = 16,
                 xlab = "Dimensions", ylab = ylab, las = 1, ...)
  if (!is.na(x$suggested)) {
    graphics::abline(v = x$suggested, lty = 2, col = "firebrick")
  }
  invisible(cv)
}

# criteria --------------------------------------------------------------------

.dim_auc <- function(space) {
  q <- fs_quality(space)
  sug <- if (nrow(q) < 3L) attr(q, "best") else .elbow_point(q$dims, q$mSD)
  list(suggested = sug, curve = data.frame(dims = q$dims, mSD = q$mSD))
}

.dim_elbow <- function(space) {
  if (is.null(space$eig)) {
    stop("This space has no eigenvalues (method '", space$method,
         "'); use method = \"auc\" or \"stress\" instead.", call. = FALSE)
  }
  ev <- space$eig[space$eig > 0]
  ax <- seq_along(ev)
  sug <- if (length(ev) < 3L) 1L else .elbow_point(ax, ev)
  list(suggested = sug, curve = data.frame(axis = ax, eigenvalue = ev))
}

.dim_parallel <- function(space, n_perm) {
  if (space$method != "pca" || is.null(space$traits)) {
    stop("Parallel analysis requires a PCA space with stored traits.",
         call. = FALSE)
  }
  scl <- !identical(space$scale, FALSE)
  M <- as.matrix(space$traits)
  cmat <- if (scl) stats::cor(M) else stats::cov(M)
  obs <- eigen(cmat, symmetric = TRUE, only.values = TRUE)$values
  p <- ncol(M)
  null_eig <- matrix(NA_real_, n_perm, p)
  for (b in seq_len(n_perm)) {
    Mp <- apply(M, 2L, sample)
    cp <- if (scl) stats::cor(Mp) else stats::cov(Mp)
    null_eig[b, ] <- eigen(cp, symmetric = TRUE, only.values = TRUE)$values
  }
  thr <- apply(null_eig, 2L, stats::quantile, probs = 0.95)
  keep <- obs > thr
  sug <- if (!keep[1L]) NA_integer_ else
    which(!c(keep, FALSE))[1L] - 1L  # leading run of TRUE
  list(suggested = sug,
       curve = data.frame(axis = seq_len(p),
                          eigenvalue = obs,
                          null95 = thr))
}

.dim_end <- function(space, q) {
  if (space$method != "pca" || is.null(space$eig)) {
    stop("The effective number of dimensions is computed from PCA ",
         "eigenvalues; build the space with method = \"pca\".",
         call. = FALSE)
  }
  ev <- space$eig[space$eig > 0]
  end <- .end_hill(ev, q)
  list(suggested = max(1L, as.integer(round(end))),
       end = end, q = q,
       curve = data.frame(axis = seq_along(ev), eigenvalue = ev,
                          share = ev / sum(ev)))
}

# Hill number of order q of the eigenvalue shares (END; Beccari &
# Carmona 2024, Nat Commun). q = 2: inverse Simpson.
.end_hill <- function(ev, q) {
  ev <- ev[ev > sqrt(.Machine$double.eps)]
  p <- ev / sum(ev)
  if (abs(q - 1) < 1e-10) return(exp(-sum(p * log(p))))
  sum(p^q)^(1 / (1 - q))
}

.dim_stress <- function(space, k_max) {
  d <- space$dist
  if (is.null(d)) {
    stop("Stress scanning requires a space with a stored distance matrix ",
         "(PCoA or NMDS built by fs_space(), or as_fspace() on a dist).",
         call. = FALSE)
  }
  n <- attr(d, "Size")
  if (is.null(k_max)) k_max <- max(1L, min(4L, n - 2L))
  stress <- vapply(seq_len(k_max), function(k) {
    MASS::isoMDS(d, k = k, trace = FALSE)$stress / 100
  }, numeric(1L))
  ok <- which(stress <= 0.10)
  sug <- if (length(ok)) ok[1L] else NA_integer_
  list(suggested = sug,
       curve = data.frame(dims = seq_len(k_max), stress = stress))
}

# elbow as the point with maximum perpendicular distance to the chord ---------

.elbow_point <- function(x, y) {
  n <- length(x)
  x1 <- x[1L]; y1 <- y[1L]; x2 <- x[n]; y2 <- y[n]
  num <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1)
  den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  x[which.max(num / den)]
}
