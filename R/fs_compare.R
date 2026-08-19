# fs_compare(): compare two trait spaces --------------------------------------

#' Compare two trait spaces
#'
#' Quantifies the agreement between two spaces containing (at least
#' partially) the same units, e.g. spaces built with different methods,
#' different trait sets, or with/without imputed values.
#'
#' * `"procrustes"`: symmetric Procrustes analysis on the common units
#'   (both configurations centred and scaled to unit sum of squares).
#'   Reports the Procrustes statistic `m2` (0 = identical configurations)
#'   and the correlation `sqrt(1 - m2)`.
#' * `"correlation"`: Pearson correlation between the pairwise distance
#'   vectors of the two spaces (a Mantel-type statistic, without the test).
#'
#' @param space1,space2 `fspace` objects with overlapping unit names.
#' @param method Comparison method.
#' @param dims Number of dimensions used from each space (default: the
#'   smaller of the two).
#'
#' @return A list of class `fs_compare` with the statistic(s), the method,
#'   and the number of shared units.
#' @seealso [fs_quality()]
#' @examples
#' data(gspff)
#' s1 <- fs_space(gspff[1:200, ], method = "pca")
#' s2 <- fs_space(fs_dist(gspff[1:200, ]), method = "pcoa")
#' fs_compare(s1, s2)
#' fs_compare(s1, s2, method = "correlation")
#' @export
fs_compare <- function(space1, space2,
                       method = c("procrustes", "correlation"),
                       dims = NULL) {
  stopifnot(is_fspace(space1), is_fspace(space2))
  method <- match.arg(method)
  shared <- intersect(rownames(space1$coords), rownames(space2$coords))
  if (length(shared) < 3L) {
    stop("The two spaces share fewer than 3 unit names.", call. = FALSE)
  }
  k <- min(ncol(space1$coords), ncol(space2$coords))
  if (!is.null(dims)) k <- min(k, as.integer(dims))
  X <- space1$coords[shared, seq_len(k), drop = FALSE]
  Y <- space2$coords[shared, seq_len(k), drop = FALSE]

  out <- switch(method,
    procrustes = {
      Xc <- base::scale(X, center = TRUE, scale = FALSE)
      Yc <- base::scale(Y, center = TRUE, scale = FALSE)
      Xn <- Xc / sqrt(sum(Xc^2))
      Yn <- Yc / sqrt(sum(Yc^2))
      m2 <- 1 - sum(svd(crossprod(Xn, Yn))$d)^2
      m2 <- max(0, min(1, m2))
      list(m2 = m2, correlation = sqrt(1 - m2))
    },
    correlation = {
      r <- stats::cor(as.vector(stats::dist(X)), as.vector(stats::dist(Y)))
      list(correlation = r)
    }
  )
  out <- c(out, list(method = method, n_shared = length(shared), dims = k))
  class(out) <- "fs_compare"
  out
}

#' @export
print.fs_compare <- function(x, ...) {
  cat("<fs_compare> method: ", x$method,
      " | shared units: ", x$n_shared,
      " | dimensions: ", x$dims, "\n", sep = "")
  if (x$method == "procrustes") {
    cat("Procrustes m2 = ", round(x$m2, 4),
        " | correlation = ", round(x$correlation, 4), "\n", sep = "")
  } else {
    cat("Distance correlation = ", round(x$correlation, 4), "\n", sep = "")
  }
  invisible(x)
}
