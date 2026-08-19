# fs_structure(): functional structure of assemblages -------------------------

#' Functional structure of assemblages
#'
#' Computes functional structure indices per assemblage, with two engines
#' behind the same interface:
#'
#' * `engine = "prob"` (default when TPDs are present): assemblage TPDs are
#'   built as abundance-weighted mixtures of unit TPDs (the TPDc logic;
#'   Carmona et al. 2016, 2019) and indices are computed from them.
#'   Between-unit dissimilarities are overlap-based (1 minus the shared
#'   probability mass).
#' * `engine = "points"`: classical point-based indices on the unit
#'   coordinates (convex-hull richness, FDis, Rao, MPD, FEve, FDiv, CWM).
#'
#' Available indices (same names in both engines; `NA` where an index has
#' no meaning in an engine): `richness` (occupied volume / hull volume),
#' `evenness`, `divergence`, `dispersion`, `rao`, `mpd`, `originality`,
#' `redundancy` (prob engine only), and `cwm` (one column per dimension).
#'
#' @param space An `fspace`; for the probabilistic engine it must carry
#'   TPDs (see [fs_tpd()]).
#' @param comm Assemblage composition: a matrix (assemblages x units,
#'   units as columns, the canonical format) or a long data.frame with
#'   columns assemblage, unit, abundance. `NULL` treats all units as one
#'   assemblage with equal abundances.
#' @param engine `"prob"` or `"points"`.
#' @param indices Character vector choosing which indices to compute
#'   (default: all available for the engine).
#' @param relative Logical; normalize abundances within assemblages to sum
#'   to 1 (default `TRUE`).
#'
#' @return A data.frame of class `fstructure` (assemblages x indices) with
#'   attributes `engine` and `settings`.
#' @references Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J.
#'   (2016) *Trends in Ecology & Evolution*, 31, 382-394. Carmona, C.P.
#'   et al. (2019) *Ecology*, 100, e02876. Villeger, S., Mason, N.W.H. &
#'   Mouillot, D. (2008) *Ecology*, 89, 2290-2301. Laliberte, E. &
#'   Legendre, P. (2010) *Ecology*, 91, 299-305.
#' @seealso [fs_tpd()], [fs_space()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 2, 4), shape = c(0, 1, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)
#' comm <- rbind(A = c(sp1 = 1, sp2 = 1, sp3 = 0),
#'               B = c(sp1 = 1, sp2 = 1, sp3 = 1))
#'
#' fs_structure(tp, comm)                      # probabilistic engine
#' fs_structure(sp, comm, engine = "points")   # point-based engine
#' @export
fs_structure <- function(space, comm = NULL,
                         engine = c("prob", "points"),
                         indices = NULL, relative = TRUE) {
  stopifnot(is_fspace(space))
  engine <- match.arg(engine)
  all_idx <- c("richness", "evenness", "divergence", "dispersion",
               "rao", "mpd", "originality", "redundancy", "cwm")
  if (is.null(indices)) indices <- all_idx
  indices <- match.arg(indices, all_idx, several.ok = TRUE)

  if (engine == "prob") {
    if (is.null(space$tpds)) {
      stop("The probabilistic engine needs TPDs: run fs_tpd() first, or ",
           "use engine = \"points\".", call. = FALSE)
    }
    unit_names <- names(space$tpds$units)
  } else {
    unit_names <- rownames(space$coords)
  }

  W <- .as_comm(comm, unit_names, relative)
  out <- if (engine == "prob") {
    .structure_prob(space, W, indices)
  } else {
    .structure_points(space, W, indices)
  }
  attr(out, "engine") <- engine
  attr(out, "settings") <- list(relative = relative, indices = indices)
  class(out) <- c("fstructure", "data.frame")
  out
}

#' @export
print.fstructure <- function(x, ...) {
  cat("<fstructure> engine: ", attr(x, "engine"),
      " | assemblages: ", nrow(x), "\n", sep = "")
  print.data.frame(utils::head(cbind(assemblage = rownames(x),
                                     round(x, 4)), 10L),
                   row.names = FALSE)
  if (nrow(x) > 10L) cat("... and", nrow(x) - 10L, "more assemblages\n")
  invisible(x)
}

# community matrix handling ----------------------------------------------------

.as_comm <- function(comm, unit_names, relative) {
  if (is.null(comm)) {
    W <- matrix(1, 1L, length(unit_names),
                dimnames = list("all", unit_names))
  } else if (is.data.frame(comm) && ncol(comm) == 3L &&
             !is.numeric(comm[[1L]])) {
    # long format: assemblage, unit, abundance
    a <- as.character(comm[[1L]])
    u <- as.character(comm[[2L]])
    ab <- comm[[3L]]
    W <- matrix(0, length(unique(a)), length(unit_names),
                dimnames = list(unique(a), unit_names))
    unknown <- setdiff(unique(u), unit_names)
    if (length(unknown)) {
      stop("Unknown unit(s) in `comm`: ",
           paste(utils::head(unknown, 5L), collapse = ", "), call. = FALSE)
    }
    for (i in seq_along(a)) W[a[i], u[i]] <- W[a[i], u[i]] + ab[i]
  } else {
    W <- as.matrix(comm)
    if (is.null(colnames(W))) {
      stop("`comm` must have unit names as column names.", call. = FALSE)
    }
    unknown <- setdiff(colnames(W), unit_names)
    if (length(unknown)) {
      stop("Unknown unit(s) in `comm`: ",
           paste(utils::head(unknown, 5L), collapse = ", "), call. = FALSE)
    }
    if (is.null(rownames(W))) {
      rownames(W) <- paste0("assemblage", seq_len(nrow(W)))
    }
    # complete missing units with zeros, in canonical order
    full <- matrix(0, nrow(W), length(unit_names),
                   dimnames = list(rownames(W), unit_names))
    full[, colnames(W)] <- W
    W <- full
  }
  if (any(W < 0) || anyNA(W)) {
    stop("Abundances must be non-negative and complete.", call. = FALSE)
  }
  rs <- rowSums(W)
  if (any(rs == 0)) {
    stop("Assemblage(s) with zero total abundance: ",
         paste(utils::head(rownames(W)[rs == 0], 5L), collapse = ", "),
         call. = FALSE)
  }
  if (relative) W <- W / rs
  W
}

# probabilistic engine ---------------------------------------------------------

.structure_prob <- function(space, W, indices) {
  tp <- space$tpds
  grid <- tp$grid
  d <- grid$d
  cellvol <- grid$cell_volume
  units <- names(tp$units)

  need_dissim <- any(c("rao", "mpd", "originality") %in% indices)
  D <- if (need_dissim) .overlap_dissim(tp$units) else NULL

  res <- vector("list", nrow(W))
  for (a in seq_len(nrow(W))) {
    w <- W[a, ]
    present <- which(w > 0)
    # assemblage TPD: abundance-weighted mixture on the union of cells
    mix <- list()
    for (u in present) {
      tu <- tp$units[[u]]
      mix[[length(mix) + 1L]] <- list(cells = tu$cells,
                                      p = tu$probs * w[u])
    }
    all_cells <- sort(unique(unlist(lapply(mix, `[[`, "cells"))))
    pc <- numeric(length(all_cells))
    n_units_cell <- integer(length(all_cells))
    pos <- seq_along(all_cells)
    names(pos) <- all_cells
    for (m in mix) {
      j <- pos[as.character(m$cells)]
      pc[j] <- pc[j] + m$p
      n_units_cell[j] <- n_units_cell[j] + 1L
    }
    pc <- pc / sum(pc)
    coords_c <- grid$cells[all_cells, , drop = FALSE]

    row <- list()
    if ("richness" %in% indices) {
      row$richness <- length(all_cells) * cellvol
    }
    if ("evenness" %in% indices) {
      M <- length(all_cells)
      row$evenness <- sum(pmin(pc, 1 / M))
    }
    if ("dispersion" %in% indices || "cwm" %in% indices) {
      # density-weighted centre: the community-weighted mean position
      cog_w <- colSums(coords_c * pc)
      if ("dispersion" %in% indices) {
        dev_w <- sweep(coords_c, 2L, cog_w)
        row$dispersion <- sum(pc * sqrt(rowSums(dev_w^2)))
      }
      if ("cwm" %in% indices) {
        for (j in seq_len(d)) row[[paste0("cwm_", j)]] <- cog_w[j]
      }
    }
    if ("divergence" %in% indices) {
      # TPD::REND convention: centre of gravity is the UNWEIGHTED mean of
      # the occupied cells; densities weight only the deviances
      cog_u <- colMeans(coords_c)
      dcog <- sqrt(rowSums(sweep(coords_c, 2L, cog_u)^2))
      dbar <- mean(dcog)
      deltaD <- sum(pc * (dcog - dbar))
      deltaAbs <- sum(pc * abs(dcog - dbar))
      row$divergence <- (deltaD + dbar) / (deltaAbs + dbar)
    }
    if (need_dissim) {
      wp <- w[present]
      Dp <- D[present, present, drop = FALSE]
      ww <- outer(wp, wp)
      diag(ww) <- 0
      if ("rao" %in% indices) row$rao <- sum(outer(wp, wp) * Dp)
      if ("mpd" %in% indices) {
        row$mpd <- if (length(present) > 1L) sum(ww * Dp) / sum(ww) else NA_real_
      }
      if ("originality" %in% indices) {
        if (length(present) > 1L) {
          orig_u <- vapply(seq_along(present), function(i) {
            sum(wp[-i] * Dp[i, -i]) / sum(wp[-i])
          }, numeric(1L))
          row$originality <- sum(wp * orig_u)
        } else {
          row$originality <- NA_real_
        }
      }
    }
    if ("redundancy" %in% indices) {
      # abundance-weighted mean number of units co-occupying a cell, - 1
      row$redundancy <- sum(pc * n_units_cell) - 1
    }
    res[[a]] <- as.data.frame(row)
  }
  out <- do.call(rbind, res)
  rownames(out) <- rownames(W)
  out
}

# overlap-based dissimilarity between unit TPDs (1 - shared mass) -------------

.overlap_dissim <- function(units) {
  n <- length(units)
  D <- matrix(0, n, n, dimnames = list(names(units), names(units)))
  for (i in seq_len(n - 1L)) {
    ci <- units[[i]]$cells
    pi_ <- units[[i]]$probs
    for (j in seq(i + 1L, n)) {
      cj <- units[[j]]$cells
      shared <- intersect(ci, cj)
      ov <- if (length(shared)) {
        sum(pmin(pi_[match(shared, ci)],
                 units[[j]]$probs[match(shared, cj)]))
      } else 0
      D[i, j] <- D[j, i] <- 1 - ov
    }
  }
  D
}

# point-based engine -----------------------------------------------------------

.structure_points <- function(space, W, indices) {
  co <- space$coords
  d <- ncol(co)
  Dfull <- if (any(c("rao", "mpd", "originality", "evenness") %in% indices)) {
    as.matrix(stats::dist(co))
  } else NULL

  res <- vector("list", nrow(W))
  for (a in seq_len(nrow(W))) {
    w <- W[a, ]
    present <- which(w > 0)
    S <- length(present)
    X <- co[present, , drop = FALSE]
    wp <- w[present]

    row <- list()
    if ("richness" %in% indices) {
      row$richness <- .hull_volume(X, d)
    }
    if ("dispersion" %in% indices || "cwm" %in% indices ||
        "divergence" %in% indices) {
      centroid <- colSums(X * wp)
      dev <- sweep(X, 2L, centroid)
      dcent <- sqrt(rowSums(dev^2))
      if ("dispersion" %in% indices) row$dispersion <- sum(wp * dcent)
      if ("cwm" %in% indices) {
        for (j in seq_len(d)) row[[paste0("cwm_", j)]] <- centroid[j]
      }
      if ("divergence" %in% indices) {
        row$divergence <- .fdiv_points(X, wp, d)
      }
    }
    if (!is.null(Dfull)) {
      Dp <- Dfull[present, present, drop = FALSE]
      ww <- outer(wp, wp)
      diag(ww) <- 0
      if ("rao" %in% indices) row$rao <- sum(outer(wp, wp) * Dp)
      if ("mpd" %in% indices) {
        row$mpd <- if (S > 1L) sum(ww * Dp) / sum(ww) else NA_real_
      }
      if ("originality" %in% indices) {
        if (S > 1L) {
          orig_u <- vapply(seq_len(S), function(i) {
            sum(wp[-i] * Dp[i, -i]) / sum(wp[-i])
          }, numeric(1L))
          row$originality <- sum(wp * orig_u)
        } else {
          row$originality <- NA_real_
        }
      }
      if ("evenness" %in% indices) {
        row$evenness <- .feve_points(Dp, wp, S)
      }
    }
    if ("redundancy" %in% indices) row$redundancy <- NA_real_
    res[[a]] <- as.data.frame(row)
  }
  out <- do.call(rbind, res)
  rownames(out) <- rownames(W)
  out
}

.hull_volume <- function(X, d) {
  if (nrow(X) <= d) return(NA_real_)
  vol <- tryCatch(geometry::convhulln(X, options = "FA")$vol,
                  error = function(e) NA_real_)
  vol
}

.fdiv_points <- function(X, wp, d) {
  if (nrow(X) <= d) return(NA_real_)
  hull <- tryCatch(geometry::convhulln(X), error = function(e) NULL)
  if (is.null(hull)) return(NA_real_)
  verts <- unique(as.vector(hull))
  cog <- colMeans(X[verts, , drop = FALSE])
  dcog <- sqrt(rowSums(sweep(X, 2L, cog)^2))
  dbar <- mean(dcog)
  deltaD <- sum(wp * (dcog - dbar))
  deltaAbs <- sum(wp * abs(dcog - dbar))
  (deltaD + dbar) / (deltaAbs + dbar)
}

.feve_points <- function(Dp, wp, S) {
  if (S < 3L) return(NA_real_)
  edges <- .mst_edges(Dp)
  ew <- apply(edges, 1L, function(e) {
    Dp[e[1L], e[2L]] / (wp[e[1L]] + wp[e[2L]])
  })
  pew <- ew / sum(ew)
  thr <- 1 / (S - 1L)
  (sum(pmin(pew, thr)) - thr) / (1 - thr)
}

# Prim's minimum spanning tree on a distance matrix ---------------------------

.mst_edges <- function(D) {
  n <- nrow(D)
  in_tree <- c(1L)
  out_tree <- setdiff(seq_len(n), in_tree)
  edges <- matrix(0L, n - 1L, 2L)
  for (k in seq_len(n - 1L)) {
    sub <- D[in_tree, out_tree, drop = FALSE]
    m <- which(sub == min(sub), arr.ind = TRUE)[1L, ]
    from <- in_tree[m[1L]]
    to <- out_tree[m[2L]]
    edges[k, ] <- c(from, to)
    in_tree <- c(in_tree, to)
    out_tree <- setdiff(out_tree, to)
  }
  edges
}
