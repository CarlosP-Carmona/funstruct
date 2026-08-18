# fs_adequacy(): can the data support density estimation at this D? -----------

# Sample sizes required to estimate a standard multivariate normal density
# at the origin with relative mean squared error < 0.1 (Silverman 1986,
# Table 4.2), for dimensions 1 to 10.
.SILVERMAN_N <- c(4, 19, 67, 223, 768, 2790, 10700, 43700, 187000, 842000)

#' Sample-size adequacy for density estimation
#'
#' The number of observations needed for kernel density estimation grows
#' exponentially with dimensionality: roughly 4 observations suffice in one
#' dimension, ~220 are needed in four, and ~2800 in six (Silverman 1986).
#' This diagnostic compares the available observations per unit with those
#' requirements, to inform the choice of dimensionality *before* estimating
#' TPDs. Point-based methods degrade over the same gradient, only silently
#' (convex hulls degenerate and pairwise distances concentrate), so an
#' inadequate `n` is a problem for every engine, not only the probabilistic
#' one.
#'
#' @param space Optional `fspace`; if supplied, `d` defaults to its current
#'   number of dimensions.
#' @param n Observations available per unit (a single value or a vector,
#'   e.g. individuals per species).
#' @param d Dimensionality to assess (1 to 10).
#'
#' @return A list of class `fs_adequacy` with elements `d`, `n_required`,
#'   `n`, `adequate` (logical, per element of `n`), and `table` (the full
#'   requirement table). `print()` gives a plain-language verdict.
#' @references Silverman, B.W. (1986) *Density Estimation for Statistics
#'   and Data Analysis*. Chapman & Hall. (Table 4.2.)
#' @seealso [fs_dimensionality()], [fs_quality()]
#' @examples
#' fs_adequacy(n = 50, d = 2)   # feasible
#' fs_adequacy(n = 50, d = 4)   # not feasible
#'
#' # per-species sample sizes from the grassland data
#' data(grassland)
#' fs_adequacy(n = table(grassland$species), d = 2)
#' @export
fs_adequacy <- function(space = NULL, n = NULL, d = NULL) {
  if (!is.null(space)) {
    stopifnot(is_fspace(space))
    if (is.null(d)) d <- ncol(space$coords)
  }
  if (is.null(d)) {
    stop("Supply `d` (or a space to take it from).", call. = FALSE)
  }
  d <- as.integer(d)
  if (length(d) != 1L || is.na(d) || d < 1L) {
    stop("`d` must be a single positive integer.", call. = FALSE)
  }
  tab <- data.frame(dims = seq_along(.SILVERMAN_N),
                    n_required = .SILVERMAN_N)
  if (d > length(.SILVERMAN_N)) {
    warning("Requirements are tabulated up to 10 dimensions; at d = ", d,
            " density estimation is not realistic for ecological sample ",
            "sizes.", call. = FALSE)
    n_req <- Inf
  } else {
    n_req <- .SILVERMAN_N[d]
  }
  adequate <- if (is.null(n)) NULL else n >= n_req
  out <- list(d = d, n_required = n_req, n = n, adequate = adequate,
              table = tab)
  class(out) <- "fs_adequacy"
  out
}

#' @export
print.fs_adequacy <- function(x, ...) {
  cat("<fs_adequacy>\n")
  cat("Estimating a density in ", x$d, " dimension(s) requires roughly ",
      format(x$n_required, big.mark = ","),
      " observations per unit (Silverman 1986).\n", sep = "")
  if (!is.null(x$n)) {
    n_ok <- sum(x$adequate)
    n_tot <- length(x$adequate)
    if (n_tot == 1L) {
      cat("With n = ", x$n, ": ",
          if (x$adequate) "adequate." else "NOT adequate.", "\n", sep = "")
    } else {
      cat(n_ok, " of ", n_tot, " units meet the requirement ",
          "(n ranges ", min(x$n), "-", max(x$n), ").\n", sep = "")
    }
    if (n_ok < n_tot) {
      cat("Consider fewer dimensions, pooling observations, or ",
          "interpreting results with caution.\n", sep = "")
    }
  }
  invisible(x)
}
