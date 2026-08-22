# fs_sensitivity(): bandwidth sensitivity analysis ----------------------------

#' Bandwidth sensitivity of functional structure
#'
#' Recomputes functional structure with all unit bandwidths scaled by a
#' range of factors (default 0.5 to 1.5, i.e. +/- 50%), packaging the
#' sensitivity protocol of Tordoni et al.: if conclusions survive halving
#' and inflating the smoothing, they do not hinge on the bandwidth.
#' Bandwidths are scaled, never re-selected.
#'
#' @param space An `fspace` with TPDs estimated by [fs_tpd()] (which
#'   stores the observations needed to re-evaluate the densities).
#' @param comm Community matrix or long data.frame; see [fs_structure()].
#' @param range Multiplicative range for the kernel standard deviations
#'   (default `c(0.5, 1.5)`).
#' @param steps Number of factors evaluated across `range` (default 5).
#' @param indices Passed to [fs_structure()].
#'
#' @return A list of class `fs_sensitivity`: `values` (long data.frame:
#'   factor, assemblage, one column per index) and `stability` (per
#'   index: correlation of assemblage values between the two extreme
#'   factors, and the mean relative range across factors). `print()` and
#'   `plot()` methods are provided.
#' @seealso [fs_tpd()], [fs_structure()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 2, 4), shape = c(0, 1, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)
#'
#' sen <- fs_sensitivity(tp, steps = 3, indices = "richness")
#' sen
#' plot(sen, index = "richness")
#' @export
fs_sensitivity <- function(space, comm = NULL, range = c(0.5, 1.5),
                           steps = 5L, indices = NULL) {
  stopifnot(is_fspace(space))
  if (is.null(space$tpds) || is.null(space$tpds$X)) {
    stop("fs_sensitivity() needs TPDs estimated with fs_tpd() (which ",
         "stores the observations).", call. = FALSE)
  }
  factors <- seq(range[1L], range[2L], length.out = as.integer(steps))
  grid <- space$tpds$grid
  alpha <- space$tpds$alpha

  res <- lapply(factors, function(f) {
    sp2 <- space
    units2 <- .unit_tpds(space)
    for (u in names(units2)) {
      H <- space$bw$values[[u]] * f^2
      dens <- .kde_at_cells(grid$cells, space$tpds$X[[u]], H)
      p <- dens * grid$cell_volume
      p <- p / sum(p)
      units2[[u]] <- .alpha_trim(p, alpha)
    }
    # swapped unit TPDs invalidate any aggregated levels: keep only the
    # bottom level in the internal copy
    sp2$tpds$levels <- space$tpds$levels["unit"]
    sp2$tpds$levels$unit$tpds <- units2
    st <- fs_structure(sp2, comm, engine = "prob", indices = indices)
    cbind(data.frame(factor = f, assemblage = rownames(st),
                     stringsAsFactors = FALSE),
          as.data.frame(st), row.names = NULL)
  })
  values <- do.call(rbind, res)

  idx_cols <- setdiff(colnames(values), c("factor", "assemblage"))
  f_lo <- factors[1L]
  f_hi <- factors[length(factors)]
  stability <- do.call(rbind, lapply(idx_cols, function(j) {
    v_lo <- values[values$factor == f_lo, j]
    v_hi <- values[values$factor == f_hi, j]
    r <- if (length(v_lo) > 2L && stats::sd(v_lo) > 0 &&
             stats::sd(v_hi) > 0) stats::cor(v_lo, v_hi) else NA_real_
    per_a <- split(values[[j]], values$assemblage)
    rel_range <- mean(vapply(per_a, function(v) {
      m <- mean(abs(v))
      if (m == 0) return(NA_real_)
      (max(v) - min(v)) / m
    }, numeric(1L)), na.rm = TRUE)
    data.frame(index = j, cor_extremes = r, mean_rel_range = rel_range)
  }))
  out <- list(values = values, stability = stability, factors = factors)
  class(out) <- "fs_sensitivity"
  out
}

#' @export
print.fs_sensitivity <- function(x, ...) {
  cat("<fs_sensitivity> bandwidth factors: ",
      paste(round(x$factors, 2), collapse = ", "), "\n", sep = "")
  cat("Stability across the extreme factors:\n")
  print.data.frame(cbind(x$stability["index"],
                         round(x$stability[, -1L], 3)),
                   row.names = FALSE)
  invisible(x)
}

#' @export
plot.fs_sensitivity <- function(x, index = NULL, ...) {
  idx_cols <- setdiff(colnames(x$values), c("factor", "assemblage"))
  if (is.null(index)) index <- idx_cols[1L]
  index <- match.arg(index, idx_cols)
  ag <- tapply(x$values[[index]], x$values$factor, mean, na.rm = TRUE)
  graphics::plot(as.numeric(names(ag)), ag, type = "b", pch = 16,
                 xlab = "Bandwidth factor", ylab = paste("Mean", index),
                 las = 1, ...)
  graphics::abline(v = 1, lty = 3, col = "grey60")
  invisible(ag)
}
