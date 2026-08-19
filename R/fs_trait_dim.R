# fs_trait_dim(): trait contributions to the effective dimensionality --------

#' How much does each trait add to the dimensionality of the space?
#'
#' Quantifies the information a trait contributes to a PCA trait space
#' through its effect on the effective number of dimensions (END; the
#' Hill number of order `q` of the eigenvalue shares -- with the default
#' `q = 2`, the inverse Simpson index; Beccari & Carmona 2024). For a
#' target trait, the observed contribution is
#'
#' \deqn{\Delta = END(all traits) - END(all traits except the target)}
#'
#' and its yardstick is a null model that permutes the target trait's
#' values across units: a permuted trait is uncorrelated with the rest,
#' so it adds close to one full effective dimension. Comparing the
#' observed \eqn{\Delta} with the null distribution tells you where the
#' trait sits between two poles:
#'
#' * \eqn{\Delta} close to the null (ratio near 1, `p_random` not
#'   small): the trait behaves like an independent axis of variation --
#'   the other traits do not capture it, and dropping it would discard
#'   information.
#' * \eqn{\Delta} far below the null (ratio near 0, small `p_random`):
#'   the trait is largely redundant with the existing set -- the space
#'   without it is nearly the same space.
#'
#' With `trait = NULL` the analysis runs for every trait, which screens
#' the whole battery at once. Used iteratively, this supports backward
#' selection of the simplest sufficient space: drop the trait with the
#' smallest contribution, rebuild, and repeat while the smallest
#' contribution stays below your retention criterion. The stopping
#' threshold is a scientific choice, not a statistical one, so the
#' package deliberately does not automate it. Note that contributions
#' are conditional on the current trait set: after dropping a trait, the
#' remaining contributions change, so re-run at each step rather than
#' dropping several traits in one pass.
#'
#' @param traits Units x traits matrix or data.frame (all numeric, no
#'   missing values -- impute first, see [fs_impute()]), or a PCA
#'   `fspace` whose stored traits are used.
#' @param trait Optional name of a single target trait; default `NULL`
#'   analyzes every trait in turn.
#' @param q Hill order for the END (default 2; see
#'   [fs_dimensionality()]).
#' @param n_null Number of permutations of the target trait (default
#'   199).
#' @param scale Logical; compute the PCA on correlations (default
#'   `TRUE`) or covariances.
#' @param seed Random seed for the permutations.
#'
#' @return A list of class `fs_trait_dim` with `end_full` (END of the
#'   complete space), `q`, `n_null`, and `table`: one row per target
#'   trait with `end_without` (END of the space without it), `delta`
#'   (observed contribution), `null_mean`, `null_lo`, `null_hi` (mean
#'   and 2.5/97.5 percentiles of the permutation contributions),
#'   `ratio` (`delta / null_mean`; 0 = fully redundant, 1 = adds as much
#'   as a random independent trait) and `p_random` (proportion of null
#'   contributions <= the observed one; small values mean the trait adds
#'   significantly less than a random trait). `print()` and `plot()`
#'   methods are provided.
#' @references Beccari, E. & Carmona, C.P. (2024) Aboveground and
#'   belowground sizes are aligned in the unified spectrum of plant form
#'   and function. *Nature Communications*, 15.
#'   \doi{10.1038/s41467-024-53180-x}
#' @seealso [fs_dimensionality()] (`method = "end"`), [fs_space()]
#' @examples
#' data(gspff)
#' x <- gspff[1:300, ]
#' td <- fs_trait_dim(x, n_null = 99, seed = 1)
#' td
#' plot(td)
#'
#' # a single target trait:
#' fs_trait_dim(x, trait = "ph", n_null = 99, seed = 1)
#' @export
fs_trait_dim <- function(traits, trait = NULL, q = 2, n_null = 199L,
                         scale = TRUE, seed = NULL) {
  if (is_fspace(traits)) {
    if (traits$method != "pca" || is.null(traits$traits)) {
      stop("When `traits` is an fspace it must be a PCA space with ",
           "stored traits.", call. = FALSE)
    }
    scale <- !identical(traits$scale, FALSE)
    traits <- traits$traits
  }
  M <- as.matrix(as.data.frame(traits))
  if (!is.numeric(M)) {
    stop("All traits must be numeric.", call. = FALSE)
  }
  if (anyNA(M)) {
    stop("Missing trait values found; impute first (see fs_impute()).",
         call. = FALSE)
  }
  if (ncol(M) < 3L) {
    stop("At least three traits are needed (the space without the ",
         "target must still have two or more).", call. = FALSE)
  }
  n_null <- as.integer(n_null)
  if (is.null(colnames(M))) colnames(M) <- paste0("tr", seq_len(ncol(M)))
  targets <- if (is.null(trait)) colnames(M) else trait
  bad <- setdiff(targets, colnames(M))
  if (length(bad)) {
    stop("Trait(s) not found: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  end_of <- function(X) {
    cm <- if (scale) stats::cor(X) else stats::cov(X)
    .end_hill(eigen(cm, symmetric = TRUE, only.values = TRUE)$values, q)
  }
  end_full <- end_of(M)

  rows <- lapply(targets, function(tn) {
    j <- match(tn, colnames(M))
    end_wo <- end_of(M[, -j, drop = FALSE])
    delta <- end_full - end_wo
    dnull <- vapply(seq_len(n_null), function(b) {
      Mp <- M
      Mp[, j] <- sample(Mp[, j])
      end_of(Mp) - end_wo
    }, numeric(1L))
    data.frame(
      trait = tn,
      end_without = end_wo,
      delta = delta,
      null_mean = mean(dnull),
      null_lo = unname(stats::quantile(dnull, 0.025)),
      null_hi = unname(stats::quantile(dnull, 0.975)),
      ratio = delta / mean(dnull),
      p_random = (1 + sum(dnull <= delta)) / (n_null + 1),
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  out <- list(end_full = end_full, q = q, n_null = n_null,
              scale = scale, table = tab, call = match.call())
  class(out) <- "fs_trait_dim"
  out
}

#' @export
print.fs_trait_dim <- function(x, ...) {
  cat("<fs_trait_dim> trait contributions to the effective number of ",
      "dimensions\n", sep = "")
  cat("END of the full space (q = ", x$q, "): ", round(x$end_full, 2),
      " | null: ", x$n_null, " permutations per trait\n", sep = "")
  tb <- x$table
  tb[, -1L] <- round(tb[, -1L], 3)
  print.data.frame(tb, row.names = FALSE)
  cat("ratio: 0 = redundant with the other traits, ",
      "1 = adds as much as a random trait.\n", sep = "")
  invisible(x)
}

#' @export
plot.fs_trait_dim <- function(x, ...) {
  tb <- x$table[order(x$table$delta), , drop = FALSE]
  k <- nrow(tb)
  xlim <- range(0, tb$delta, tb$null_lo, tb$null_hi)
  graphics::plot(NA, xlim = xlim, ylim = c(0.5, k + 0.5), yaxt = "n",
                 xlab = "Added effective dimensions", ylab = "", las = 1,
                 ...)
  graphics::axis(2L, at = seq_len(k), labels = tb$trait, las = 1)
  graphics::abline(v = 0, lty = 3, col = "grey60")
  graphics::segments(tb$null_lo, seq_len(k), tb$null_hi, seq_len(k),
                     col = "grey70", lwd = 3)
  graphics::points(tb$null_mean, seq_len(k), pch = 3, col = "grey40")
  graphics::points(tb$delta, seq_len(k), pch = 16, col = "firebrick")
  graphics::legend("bottomright",
                   legend = c("observed", "random-trait null (95%)"),
                   pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 3),
                   col = c("firebrick", "grey70"), bty = "n")
  invisible(tb)
}
