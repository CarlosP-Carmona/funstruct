# fs_pool(): one-step pooled assemblage TPDs (use with caution) ---------------

#' One-step pooled assemblage TPDs
#'
#' Builds one TPD per assemblage directly from the pooled individual
#' observations of its constituent units, ignoring unit identities and
#' abundances (every individual contributes equally).
#'
#' **This is usually the wrong tool.** Pooling observations without
#' accounting for unit identity introduces pseudoreplication (individuals
#' of the same species cluster in trait space), makes results depend on
#' how many units happen to be sampled, and, if bandwidths were selected
#' per pooled sample, lets the smoothing covary with sample composition
#' (Tordoni et al.; see also Gross et al. 2024 and its reanalysis). The
#' default workflow ([fs_tpd()] then [fs_structure()]) avoids all three
#' problems. `fs_pool()` exists for comparisons and for reproducing
#' pooled-observation analyses; it uses one single bandwidth selected
#' once from all observations, so at least the smoothing cannot covary
#' with composition.
#'
#' @param space An `fspace` whose rows are individual observations.
#' @param ids Vector assigning each row to a unit.
#' @param comm Community matrix (assemblages x units) or long data.frame.
#' @param bw Optional bandwidth (numeric per-axis SDs or a d x d matrix);
#'   default: one plug-in selection on all observations pooled.
#' @param grid An `fs_grid`; default `fs_grid(space)`.
#' @param alpha Probability mass retained (default 0.99).
#'
#' @return The `fspace` with `tpds` whose units are the *assemblages*
#'   (one pooled TPD each); downstream functions such as [fs_beta()] can
#'   then compare assemblages directly (use an identity community
#'   matrix).
#' @seealso [fs_tpd()] for the recommended two-step construction.
#' @export
fs_pool <- function(space, ids, comm, bw = NULL, grid = NULL,
                    alpha = 0.99) {
  stopifnot(is_fspace(space))
  d <- ncol(space$coords)
  .check_tpd_dims(d)
  ids <- as.character(ids)
  if (length(ids) != nrow(space$coords)) {
    stop("`ids` must have one entry per row of the space.", call. = FALSE)
  }
  unit_names <- sort(unique(ids))
  W <- .as_comm(comm, unit_names, relative = TRUE)
  if (is.null(grid)) grid <- fs_grid(space)

  warning("fs_pool() pools observations ignoring unit identities and ",
          "abundances; this induces pseudoreplication and richness-",
          "dependent estimates. Prefer fs_tpd() + fs_structure(). ",
          "See ?fs_pool.", call. = FALSE)

  H <- if (is.null(bw)) {
    message("Using one plug-in bandwidth selected from all observations ",
            "pooled (constant across assemblages).")
    .select_bw(space$coords, "plugin")
  } else {
    .as_H(bw, d)
  }
  if (is.null(H)) {
    stop("Bandwidth selection failed; supply `bw`.", call. = FALSE)
  }

  units_tpd <- vector("list", nrow(W))
  names(units_tpd) <- rownames(W)
  X_list <- vector("list", nrow(W))
  names(X_list) <- rownames(W)
  Hs <- vector("list", nrow(W))
  names(Hs) <- rownames(W)
  for (a in seq_len(nrow(W))) {
    present <- colnames(W)[W[a, ] > 0]
    Xa <- space$coords[ids %in% present, , drop = FALSE]
    dens <- .kde_at_cells(grid$cells, Xa, H)
    p <- dens * grid$cell_volume
    p <- p / sum(p)
    units_tpd[[a]] <- .alpha_trim(p, alpha)
    X_list[[a]] <- Xa
    Hs[[a]] <- H
  }
  space$tpds <- list(grid = grid, alpha = alpha, units = units_tpd,
                     n_obs = vapply(X_list, nrow, integer(1L)),
                     route = "pooled", X = X_list)
  space$bw <- list(attachment = "common", selector =
                     if (is.null(bw)) "plugin" else "user",
                   values = Hs,
                   imputed = stats::setNames(rep(FALSE, nrow(W)),
                                             rownames(W)))
  space
}
