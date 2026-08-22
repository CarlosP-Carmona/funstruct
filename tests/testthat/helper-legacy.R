# Frozen legacy implementations (pre-fs_aggregate refactor) -------------------
#
# These are verbatim copies of the internal mixture/dissimilarity code as of
# commit "fixed colors" (2026-08-22), with one mechanical change: they take
# `grid` and `units` (the sparse per-unit TPD list) as explicit arguments
# instead of reading them from an fspace, so they stay valid whatever the
# storage layout of the package becomes. They pin the CURRENT numerical
# behaviour of fs_structure() and fs_beta() (probabilistic engine) and
# fs_partition(); the refactored code must reproduce them exactly.
#
# DO NOT EDIT the function bodies: they are the reference.

legacy_as_comm <- function(comm, unit_names, relative) {
  if (is.null(comm)) {
    W <- matrix(1, 1L, length(unit_names),
                dimnames = list("all", unit_names))
  } else if (is.data.frame(comm) && ncol(comm) == 3L &&
             !is.numeric(comm[[1L]])) {
    a <- as.character(comm[[1L]])
    u <- as.character(comm[[2L]])
    ab <- comm[[3L]]
    W <- matrix(0, length(unique(a)), length(unit_names),
                dimnames = list(unique(a), unit_names))
    for (i in seq_along(a)) W[a[i], u[i]] <- W[a[i], u[i]] + ab[i]
  } else {
    W <- as.matrix(comm)
    if (is.null(rownames(W))) {
      rownames(W) <- paste0("assemblage", seq_len(nrow(W)))
    }
    full <- matrix(0, nrow(W), length(unit_names),
                   dimnames = list(rownames(W), unit_names))
    full[, colnames(W)] <- W
    W <- full
  }
  rs <- rowSums(W)
  if (relative) W <- W / rs
  W
}

legacy_overlap_dissim <- function(units) {
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

legacy_tpdc_sparse <- function(units, w) {
  present <- which(w > 0)
  parts <- lapply(present, function(u) {
    list(cells = units[[u]]$cells, p = units[[u]]$probs * w[u])
  })
  all_cells <- sort(unique(unlist(lapply(parts, `[[`, "cells"))))
  pc <- numeric(length(all_cells))
  pos <- seq_along(all_cells)
  names(pos) <- all_cells
  for (m in parts) {
    j <- pos[as.character(m$cells)]
    pc[j] <- pc[j] + m$p
  }
  list(cells = all_cells, probs = pc / sum(pc))
}

legacy_structure_prob <- function(grid, units_list, W, indices) {
  tp <- list(grid = grid, units = units_list)
  d <- grid$d
  cellvol <- grid$cell_volume

  need_dissim <- any(c("rao", "mpd", "originality") %in% indices)
  D <- if (need_dissim) legacy_overlap_dissim(tp$units) else NULL

  res <- vector("list", nrow(W))
  for (a in seq_len(nrow(W))) {
    w <- W[a, ]
    present <- which(w > 0)
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
      row$redundancy <- sum(pc * n_units_cell) - 1
    }
    res[[a]] <- as.data.frame(row)
  }
  out <- do.call(rbind, res)
  rownames(out) <- rownames(W)
  out
}

legacy_beta_prob <- function(units_list, W, decompose = TRUE) {
  n <- nrow(W)
  nm <- rownames(W)
  Dm <- matrix(0, n, n, dimnames = list(nm, nm))
  Psh <- Pns <- if (decompose) {
    matrix(NA_real_, n, n, dimnames = list(nm, nm))
  } else NULL

  tpdc <- lapply(seq_len(n), function(a) {
    legacy_tpdc_sparse(units_list, W[a, ])
  })
  for (i in seq_len(n - 1L)) {
    ci <- tpdc[[i]]$cells; pi_ <- tpdc[[i]]$probs
    for (j in seq(i + 1L, n)) {
      cj <- tpdc[[j]]$cells; pj_ <- tpdc[[j]]$probs
      common <- intersect(ci, cj)
      ii <- match(common, ci); jj <- match(common, cj)
      O <- if (length(common)) sum(pmin(pi_[ii], pj_[jj])) else 0
      dis <- 1 - O
      Dm[i, j] <- Dm[j, i] <- dis
      if (!is.null(Psh)) {
        if (dis > 1e-12) {
          A <- if (length(common)) sum(pmax(pi_[ii], pj_[jj])) - O else 0
          B <- 1 - if (length(common)) sum(pi_[ii]) else 0
          C <- 1 - if (length(common)) sum(pj_[jj]) else 0
          pns <- (2 * min(B, C)) / (A + 2 * min(B, C))
          Pns[i, j] <- Pns[j, i] <- pns
          Psh[i, j] <- Psh[j, i] <- 1 - pns
        }
      }
    }
  }
  list(dissimilarity = Dm, P_shared = Psh, P_non_shared = Pns)
}

legacy_sparse_max <- function(parts) {
  all_cells <- sort(unique(unlist(lapply(parts, `[[`, "cells"))))
  vmax <- numeric(length(all_cells))
  pos <- seq_along(all_cells)
  names(pos) <- all_cells
  for (p in parts) {
    j <- pos[as.character(p$cells)]
    vmax[j] <- pmax(vmax[j], p$vals)
  }
  list(cells = all_cells, vals = vmax)
}

legacy_eqv_levels <- function(units, W, H, q) {
  maxfuns <- vector("list", nrow(W))
  eqv_a <- numeric(nrow(W))
  for (a in seq_len(nrow(W))) {
    w <- W[a, ]
    present <- which(w > 0)
    wmax <- max(w[present])
    parts <- lapply(present, function(u) {
      list(cells = units[[u]]$cells,
           vals = units[[u]]$probs * (w[u] / wmax)^q)
    })
    maxfuns[[a]] <- legacy_sparse_max(parts)
    eqv_a[a] <- sum(maxfuns[[a]]$vals)
  }
  names(maxfuns) <- names(eqv_a) <- rownames(W)

  level_names <- c("assemblage", if (!is.null(H)) colnames(H), "total")
  alpha <- numeric(length(level_names))
  names(alpha) <- level_names
  alpha["assemblage"] <- mean(eqv_a)
  values <- list(assemblage = eqv_a)

  if (!is.null(H)) {
    for (l in seq_len(ncol(H))) {
      gr <- split(seq_len(nrow(W)), H[[l]])
      g_eqv <- vapply(gr, function(idx) {
        sum(legacy_sparse_max(maxfuns[idx])$vals)
      }, numeric(1L))
      alpha[colnames(H)[l]] <- mean(g_eqv)
      values[[colnames(H)[l]]] <- g_eqv
    }
  }
  gamma <- sum(legacy_sparse_max(maxfuns)$vals)
  alpha["total"] <- gamma
  values$total <- gamma

  comp <- data.frame(
    component = c(paste0("alpha_", level_names[1L]),
                  paste0("beta_", level_names[-length(level_names)],
                         ":", level_names[-1L])),
    value = c(alpha[1L], diff(alpha)))
  list(alpha = alpha, values = values, table = comp)
}

legacy_partition_rao <- function(units_list, W, H) {
  D <- legacy_overlap_dissim(units_list)
  rao_of <- function(w) as.numeric(t(w) %*% D %*% w)
  eqv_of <- function(Q) 1 / (1 - pmin(Q, 1 - 1e-12))

  rao_a <- apply(W, 1L, rao_of)
  level_names <- c("assemblage", if (!is.null(H)) colnames(H), "total")
  alpha <- numeric(length(level_names))
  names(alpha) <- level_names
  alpha["assemblage"] <- mean(rao_a)
  values <- list(assemblage = rao_a)

  if (!is.null(H)) {
    for (l in seq_len(ncol(H))) {
      gr <- split(seq_len(nrow(W)), H[[l]])
      g_rao <- vapply(gr, function(idx) {
        rao_of(colMeans(W[idx, , drop = FALSE]))
      }, numeric(1L))
      alpha[colnames(H)[l]] <- mean(g_rao)
      values[[colnames(H)[l]]] <- g_rao
    }
  }
  gamma <- rao_of(colMeans(W))
  alpha["total"] <- gamma
  values$total <- gamma

  E <- eqv_of(alpha)
  comp <- data.frame(
    component = c(paste0("alpha_", level_names[1L]),
                  paste0("beta_", level_names[-length(level_names)],
                         ":", level_names[-1L])),
    value = c(alpha[1L], diff(alpha)),
    value_eqv = c(E[1L], diff(E)),
    prop_eqv = c(E[1L], diff(E)) / E[length(E)])
  list(alpha = alpha, table = comp)
}

# storage-layout shim (the ONLY part of this file that may be edited): return
# the sparse per-unit TPD list wherever the package currently stores it.
legacy_get_units <- function(space) {
  tp <- space$tpds
  if (!is.null(tp$levels)) tp$levels$unit$tpds else tp$units
}
