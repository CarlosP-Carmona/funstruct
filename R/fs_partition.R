# fs_partition(): diversity across nested scales ------------------------------

#' Declare a hierarchy of nested scales
#'
#' Groups assemblages (the rows of the community matrix) into nested,
#' increasingly coarse levels, e.g. plots within sites within regions.
#'
#' @param ... Named vectors (one per level, ordered fine to coarse), each
#'   with one entry per assemblage.
#' @return A data.frame of class `fs_hierarchy`.
#' @seealso [fs_partition()]
#' @export
fs_hierarchy <- function(...) {
  args <- list(...)
  if (!length(args)) {
    stop("Supply at least one grouping vector.", call. = FALSE)
  }
  if (length(unique(lengths(args))) != 1L) {
    stop("All grouping vectors must have the same length (one entry per ",
         "assemblage).", call. = FALSE)
  }
  nm <- names(args)
  if (is.null(nm) || any(nm == "")) {
    nm <- paste0("level", seq_along(args))
    names(args) <- nm
  }
  out <- as.data.frame(lapply(args, as.character),
                       stringsAsFactors = FALSE)
  class(out) <- c("fs_hierarchy", "data.frame")
  out
}

#' Partition functional diversity across nested scales
#'
#' Decomposes diversity additively across any number of nested levels
#' (assemblages, then the levels declared with [fs_hierarchy()], up to the
#' whole dataset).
#'
#' * `method = "rao"`: Rao quadratic entropy with the additive partition
#'   and the equivalent-number correction of de Bello et al. (2010).
#'   Dissimilarities between units are overlap-based (1 minus shared TPD
#'   mass) when TPDs are present, and max-scaled Euclidean distances in
#'   the space otherwise. `q` is not used by this method.
#' * `method = "tpd_eqv"`: TPD-based equivalent numbers (Castro
#'   Sanchez-Bermejo et al.): each unit's TPD is rescaled by
#'   `(p/p_max)^q` (relative abundances within the assemblage) and the
#'   equivalent number is the integral of the pointwise maximum across
#'   units. Aggregation to higher levels takes the pointwise maximum of
#'   the lower-level maximum functions, without re-weighting, giving a
#'   non-negative additive partition at every level.
#' * `method = "hill"`: not yet implemented (planned; Hill-Chao framework).
#'
#' @param space An `fspace`; `"tpd_eqv"` (and `"rao"` with overlap
#'   dissimilarities) require TPDs (see [fs_tpd()]).
#' @param comm Community matrix (assemblages x units) or long data.frame;
#'   see [fs_structure()].
#' @param hierarchy An [fs_hierarchy()] with one row per assemblage (rows
#'   in the order of `comm`), or `NULL` for a two-level partition
#'   (assemblages vs the whole dataset).
#' @param method Partition method; see Details.
#' @param q Dominance parameter for `"tpd_eqv"` (default 1; `q = 0`
#'   ignores abundances).
#' @param profile Logical; for `"tpd_eqv"`, also compute the full profile
#'   over `q = seq(0, 2, 0.1)`.
#'
#' @return A list of class `fpartition`: `method`, `levels`, `values`
#'   (per-assemblage and per-group values, and gamma), `table` (the
#'   additive components), and `profiles` (when requested).
#' @references de Bello, F. et al. (2010) The partitioning of diversity:
#'   showing Theseus a way out of the labyrinth. *Journal of Vegetation
#'   Science*, 21, 992-1000.
#' @seealso [fs_hierarchy()], [fs_structure()], [fs_tpd()]
#' @export
fs_partition <- function(space, comm, hierarchy = NULL,
                         method = c("rao", "hill", "tpd_eqv"),
                         q = 1, profile = FALSE) {
  stopifnot(is_fspace(space))
  method <- match.arg(method)
  if (method == "hill") {
    stop("method = \"hill\" is not yet implemented; use \"rao\" or ",
         "\"tpd_eqv\".", call. = FALSE)
  }
  if (method == "tpd_eqv" && is.null(space$tpds)) {
    stop("method = \"tpd_eqv\" needs TPDs: run fs_tpd() first.",
         call. = FALSE)
  }
  unit_names <- if (!is.null(space$tpds)) names(space$tpds$units) else
    rownames(space$coords)
  W <- .as_comm(comm, unit_names, relative = TRUE)
  H <- .check_hierarchy(hierarchy, nrow(W))

  out <- switch(method,
    rao     = .partition_rao(space, W, H),
    tpd_eqv = .partition_eqv(space, W, H, q,
                             q_grid = if (profile) seq(0, 2, 0.1) else NULL)
  )
  out$method <- method
  out$call <- match.call()
  class(out) <- "fpartition"
  out
}

#' @export
print.fpartition <- function(x, ...) {
  cat("<fpartition> method: ", x$method,
      if (!is.null(x$q)) paste0(" | q = ", x$q), "\n", sep = "")
  cat("Levels: ", paste(x$levels, collapse = " < "), "\n", sep = "")
  print.data.frame(cbind(x$table["component"],
                         round(x$table[, -1L, drop = FALSE], 4)),
                   row.names = FALSE)
  if (!is.null(x$profiles)) {
    cat("(q-profile computed for q = ",
        min(x$profiles$q), "-", max(x$profiles$q), ")\n", sep = "")
  }
  invisible(x)
}

#' @export
plot.fpartition <- function(x, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  if (!is.null(x$profiles)) {
    pr <- x$profiles
    scales <- unique(pr$scale)
    cols <- grDevices::hcl.colors(length(scales), "Dark 3")
    graphics::plot(range(pr$q), range(pr$value), type = "n",
                   xlab = "q", ylab = "Equivalent numbers", las = 1, ...)
    for (i in seq_along(scales)) {
      z <- pr[pr$scale == scales[i], ]
      graphics::lines(z$q, z$value, col = cols[i], lwd = 2)
    }
    graphics::legend("topright", legend = scales, col = cols, lwd = 2,
                     bty = "n")
    return(invisible(pr))
  }
  bp <- x$table$value
  names(bp) <- x$table$component
  graphics::barplot(bp, las = 2, ylab = "Diversity component", ...)
  invisible(x$table)
}

# hierarchy handling -----------------------------------------------------------

.check_hierarchy <- function(hierarchy, n_assemblages) {
  if (is.null(hierarchy)) return(NULL)
  H <- as.data.frame(hierarchy)
  if (nrow(H) != n_assemblages) {
    stop("`hierarchy` must have one row per assemblage (", n_assemblages,
         "), got ", nrow(H), ".", call. = FALSE)
  }
  if (ncol(H) > 1L) {
    for (l in seq(2L, ncol(H))) {
      spans <- tapply(H[[l]], H[[l - 1L]],
                      function(x) length(unique(x)))
      if (any(spans > 1L)) {
        stop("Levels must be strictly nested: group(s) at level '",
             colnames(H)[l - 1L], "' span several groups of level '",
             colnames(H)[l], "'.", call. = FALSE)
      }
    }
  }
  H
}

# TPD equivalent numbers -------------------------------------------------------

.partition_eqv <- function(space, W, H, q, q_grid = NULL) {
  units <- space$tpds$units
  main <- .eqv_levels(units, W, H, q)
  profiles <- NULL
  if (!is.null(q_grid)) {
    prof <- lapply(q_grid, function(qq) {
      lv <- .eqv_levels(units, W, H, qq)
      data.frame(q = qq, scale = names(lv$alpha), value = unname(lv$alpha))
    })
    profiles <- do.call(rbind, prof)
  }
  list(q = q, levels = names(main$alpha),
       values = main$values,
       table = main$table, profiles = profiles)
}

.eqv_levels <- function(units, W, H, q) {
  # per-assemblage maximum functions (sparse: cells + values)
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
    maxfuns[[a]] <- .sparse_max(parts)
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
        sum(.sparse_max(maxfuns[idx])$vals)
      }, numeric(1L))
      alpha[colnames(H)[l]] <- mean(g_eqv)
      values[[colnames(H)[l]]] <- g_eqv
    }
  }
  gamma <- sum(.sparse_max(maxfuns)$vals)
  alpha["total"] <- gamma
  values$total <- gamma

  comp <- data.frame(
    component = c(paste0("alpha_", level_names[1L]),
                  paste0("beta_", level_names[-length(level_names)],
                         ":", level_names[-1L])),
    value = c(alpha[1L], diff(alpha)))
  list(alpha = alpha, values = values, table = comp)
}

.sparse_max <- function(parts) {
  # parts: list of list(cells, vals); returns pointwise maximum
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

# Rao partition (de Bello et al. 2010) -----------------------------------------

.partition_rao <- function(space, W, H) {
  if (!is.null(space$tpds)) {
    D <- .overlap_dissim(space$tpds$units)
  } else {
    D <- as.matrix(stats::dist(space$coords))
    D <- D / max(D)
  }
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

  list(q = NULL, levels = level_names, values = values,
       table = comp, profiles = NULL)
}
