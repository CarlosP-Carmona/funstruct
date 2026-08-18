# fspace class: internal constructor, validator, and methods ------------------

new_fspace <- function(coords, method, traits, units,
                       loadings = NULL, eig = NULL, stress = NULL,
                       dist = NULL, scale = NA, rotation = "none",
                       reduced = FALSE, dims_full = ncol(coords),
                       tpds = NULL, bw = NULL, call = NULL) {
  stopifnot(is.matrix(coords))
  out <- list(
    coords    = coords,
    method    = method,
    loadings  = loadings,
    eig       = eig,
    stress    = stress,
    rotation  = rotation,
    reduced   = reduced,
    dims_full = dims_full,
    traits    = traits,
    units     = units,
    dist      = dist,
    scale     = scale,
    tpds      = tpds,
    bw        = bw,
    call      = call,
    version   = as.character(utils::packageVersion("funstruct"))
  )
  class(out) <- "fspace"
  out
}

#' Test whether an object is an fspace
#'
#' @param x An object.
#' @return `TRUE` if `x` inherits from class `fspace`.
#' @examples
#' data(gspff)
#' is_fspace(fs_space(gspff[1:50, ], method = "pca"))
#' is_fspace(gspff)
#' @export
is_fspace <- function(x) inherits(x, "fspace")

#' @export
print.fspace <- function(x, ...) {
  d <- ncol(x$coords)
  cat("<fspace> ", toupper(x$method),
      if (x$rotation != "none") paste0(" (", x$rotation, "-rotated)"),
      "\n", sep = "")
  cat("  units: ", nrow(x$coords),
      " | dimensions: ", d,
      if (x$reduced) paste0(" (reduced from ", x$dims_full, ")"),
      "\n", sep = "")
  if (!is.null(x$eig)) {
    pv <- round(100 * x$eig[seq_len(min(d, 4L))] / sum(x$eig[x$eig > 0]), 1)
    cat("  variance explained (first axes): ",
        paste0(pv, "%", collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$stress)) {
    cat("  stress: ", round(x$stress, 3), "\n", sep = "")
  }
  cat("  TPDs: ", if (is.null(x$tpds)) "not estimated" else
    paste0("estimated (alpha = ", x$tpds$alpha, ")"), "\n", sep = "")
  invisible(x)
}

#' @export
summary.fspace <- function(object, ...) {
  print(object)
  if (!is.null(object$bw)) {
    cat("\nBandwidth information\n")
    cat("  attachment: ", object$bw$attachment,
        " | selector: ", object$bw$selector, "\n", sep = "")
    if (any(object$bw$imputed)) {
      imp_names <- names(object$bw$imputed)[object$bw$imputed]
      cat("  imputed widths for ", sum(object$bw$imputed), " unit(s): ",
          paste(utils::head(imp_names, 5L), collapse = ", "),
          if (sum(object$bw$imputed) > 5L) ", ...", "\n", sep = "")
    }
  }
  n_imp <- sum(object$units$imputed_traits, na.rm = TRUE)
  if (isTRUE(n_imp > 0)) {
    cat("\n", n_imp, " unit(s) contain imputed trait values.\n", sep = "")
  }
  invisible(object)
}

#' @export
plot.fspace <- function(x, dims = c(1L, 2L), color_by = NULL,
                        loadings = TRUE, cex = 0.8, ...) {
  d <- ncol(x$coords)
  if (d < 2L) stop("Plotting requires at least two dimensions.", call. = FALSE)
  dims <- as.integer(dims[1:2])
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  xy <- x$coords[, dims, drop = FALSE]
  col <- "grey30"
  if (!is.null(color_by)) {
    f <- as.factor(color_by)
    col <- grDevices::hcl.colors(nlevels(f), "Dark 3")[as.integer(f)]
  }
  lab <- .axis_labels(x, dims)
  graphics::plot(xy, col = col, pch = 16, cex = cex,
                 xlab = lab[1L], ylab = lab[2L], las = 1, ...)
  graphics::abline(h = 0, v = 0, lty = 3, col = "grey70")

  if (isTRUE(loadings) && !is.null(x$loadings)) {
    ld <- x$loadings[, dims, drop = FALSE]
    sc <- 0.8 * min(apply(abs(xy), 2L, max)) / max(abs(ld))
    graphics::arrows(0, 0, ld[, 1L] * sc, ld[, 2L] * sc,
                     length = 0.08, col = "firebrick")
    graphics::text(ld[, 1L] * sc * 1.08, ld[, 2L] * sc * 1.08,
                   labels = rownames(ld), col = "firebrick", cex = 0.8)
  }
  invisible(list(coords = xy,
                 loadings = if (!is.null(x$loadings))
                   x$loadings[, dims, drop = FALSE]))
}

.axis_labels <- function(x, dims) {
  if (!is.null(x$eig)) {
    pv <- round(100 * x$eig / sum(x$eig[x$eig > 0]), 1)
    paste0(toupper(x$method), " ", dims, " (", pv[dims], "%)")
  } else {
    paste0("Axis ", dims)
  }
}
