# Novel indices: redundancy profiles and ITV contribution ---------------------

#' Continuous functional redundancy profiles
#'
#' *Novel index.* Extends the TPD redundancy of [fs_structure()]
#' (abundance-weighted mean number of units co-occupying a cell, minus 1)
#' into a profile over the dominance parameter `q`: abundances are raised
#' to `q` and renormalized before mixing, so `q = 0` weighs all present
#' units equally and larger `q` concentrates the assemblage TPD on
#' dominant units. The profile shows whether redundancy is carried by the
#' whole assemblage or only by dominant units, linking redundancy to
#' functional insurance.
#'
#' @param space An `fspace` with TPDs ([fs_tpd()]).
#' @param comm Community matrix or long data.frame; see [fs_structure()].
#' @param q Numeric vector of dominance parameters (default
#'   `seq(0, 2, 0.1)`).
#'
#' @return A list of class `fs_redundancy`: `values` (long data.frame:
#'   assemblage, q, redundancy). `print()` and `plot()` methods provided.
#' @seealso [fs_structure()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 0.3, 4), shape = c(0, 0.2, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)
#' comm <- rbind(A = c(sp1 = 0.6, sp2 = 0.3, sp3 = 0.1))
#'
#' red <- fs_redundancy(tp, comm, q = c(0, 1, 2))
#' red
#' plot(red)
#' @export
fs_redundancy <- function(space, comm = NULL, q = seq(0, 2, 0.1)) {
  stopifnot(is_fspace(space))
  if (is.null(space$tpds)) {
    stop("fs_redundancy() needs TPDs: run fs_tpd() first.", call. = FALSE)
  }
  units <- .unit_tpds(space)
  W <- .as_comm(comm, names(units), relative = TRUE)
  rows <- list()
  for (a in seq_len(nrow(W))) {
    w <- W[a, ]
    present <- which(w > 0)
    for (qq in q) {
      wq <- w[present]^qq
      wq <- wq / sum(wq)
      wfull <- numeric(length(units))
      wfull[present] <- wq
      mixed <- .mix_tpds(units, wfull, counts = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        assemblage = rownames(W)[a], q = qq,
        redundancy = sum(mixed$probs * mixed$counts) - 1)
    }
  }
  out <- list(values = do.call(rbind, rows), q = q)
  class(out) <- "fs_redundancy"
  out
}

#' @export
print.fs_redundancy <- function(x, ...) {
  cat("<fs_redundancy> profiles over q = ",
      min(x$q), "-", max(x$q), "\n", sep = "")
  at1 <- x$values[abs(x$values$q - 1) < 1e-9, ]
  if (nrow(at1)) {
    cat("Redundancy at q = 1:\n")
    print.data.frame(cbind(at1["assemblage"],
                           redundancy = round(at1$redundancy, 3)),
                     row.names = FALSE)
  }
  invisible(x)
}

#' @export
plot.fs_redundancy <- function(x, ...) {
  v <- x$values
  assemblages <- unique(v$assemblage)
  cols <- grDevices::hcl.colors(max(2L, length(assemblages)), "Dark 3")
  graphics::plot(range(v$q), range(v$redundancy), type = "n",
                 xlab = "q", ylab = "Functional redundancy", las = 1, ...)
  for (i in seq_along(assemblages)) {
    z <- v[v$assemblage == assemblages[i], ]
    graphics::lines(z$q, z$redundancy, col = cols[i], lwd = 2)
  }
  if (length(assemblages) <= 8L) {
    graphics::legend("topright", legend = assemblages,
                     col = cols[seq_along(assemblages)], lwd = 2,
                     bty = "n")
  }
  invisible(v)
}

#' Contribution of intraspecific variability to functional structure
#'
#' *Novel index.* Quantifies how much of each assemblage's functional
#' structure is contributed by within-unit (intraspecific) trait
#' variability. Structure is computed twice with identical bandwidths:
#' from the full unit TPDs (individual observations) and from TPDs
#' collapsed to each unit's mean position. The contribution is
#' `(full - mean_based) / full` per index and assemblage: 0 means the
#' index is fully explained by mean trait positions, values toward 1 mean
#' within-unit variability dominates. Requires TPDs estimated from
#' observations ([fs_tpd()] with `ids`).
#'
#' @param space An `fspace` with observation-based TPDs.
#' @param comm Community matrix or long data.frame; see [fs_structure()].
#' @param indices Indices to evaluate (default `"richness"` and
#'   `"dispersion"`, where the interpretation is cleanest).
#'
#' @return A list of class `fs_itv_contribution`: `full`, `mean_based`
#'   (both `fstructure`) and `contribution` (data.frame, same shape).
#' @seealso [fs_itv()] for the effect of ITV on the space itself.
#' @examples
#' set.seed(1)
#' tr <- data.frame(t1 = rnorm(60, rep(c(0, 3, 6), each = 20), 0.5),
#'                  t2 = rnorm(60, rep(c(0, 3, 0), each = 20), 0.5),
#'                  row.names = paste0("i", 1:60))
#' ids <- rep(c("a", "b", "c"), each = 20)
#' sp <- fs_space(tr, method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, ids = ids)
#'
#' fs_itv_contribution(tp)
#' @export
fs_itv_contribution <- function(space, comm = NULL,
                                indices = c("richness", "dispersion")) {
  stopifnot(is_fspace(space))
  if (is.null(space$tpds) || !identical(space$tpds$route, "obs")) {
    stop("fs_itv_contribution() needs TPDs estimated from individual ",
         "observations (fs_tpd() with `ids`).", call. = FALSE)
  }
  full <- fs_structure(space, comm, engine = "prob", indices = indices)

  grid <- space$tpds$grid
  alpha <- space$tpds$alpha
  spm <- space
  units_m <- .unit_tpds(space)
  for (u in names(units_m)) {
    Xm <- matrix(colMeans(space$tpds$X[[u]]), 1L, ncol(space$coords))
    dens <- .kde_at_cells(grid$cells, Xm, space$bw$values[[u]])
    p <- dens * grid$cell_volume
    p <- p / sum(p)
    units_m[[u]] <- .alpha_trim(p, alpha)
  }
  # swapped unit TPDs invalidate any aggregated levels: keep only the
  # bottom level in the internal copy
  spm$tpds$levels <- space$tpds$levels["unit"]
  spm$tpds$levels$unit$tpds <- units_m
  mean_based <- fs_structure(spm, comm, engine = "prob",
                             indices = indices)

  fm <- as.matrix(full)
  mm <- as.matrix(mean_based)
  contribution <- (fm - mm) / fm
  contribution[!is.finite(contribution)] <- NA_real_
  out <- list(full = full, mean_based = mean_based,
              contribution = as.data.frame(contribution))
  class(out) <- "fs_itv_contribution"
  out
}

#' @export
print.fs_itv_contribution <- function(x, ...) {
  cat("<fs_itv_contribution> share of structure from within-unit ",
      "variability\n", sep = "")
  print.data.frame(round(x$contribution, 3))
  invisible(x)
}
