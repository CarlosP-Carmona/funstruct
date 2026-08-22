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
#' @examples
#' # four assemblages in two sites within one region
#' fs_hierarchy(site = c("s1", "s1", "s2", "s2"),
#'              region = c("r1", "r1", "r1", "r1"))
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
#' Decomposes diversity additively across any number of nested scales.
#' The scales can be declared in two mutually exclusive ways:
#'
#' * `comm` (+ optional `hierarchy`): assemblages are the rows of the
#'   community matrix, grouped by the [fs_hierarchy()] levels, up to the
#'   whole dataset.
#' * `levels`: a character vector of stored TPD levels (fine to coarse;
#'   see [fs_aggregate()]). The first element is the base level whose
#'   TPDs (and dissimilarities) underlie the calculation, the second is
#'   the assemblage scale, and any further elements are increasingly
#'   coarse grouping scales (they must be strictly nested and lie on the
#'   stored aggregation chain). The effective assemblage compositions are
#'   taken from [fs_level_weights()], so both routes feed the exact same
#'   engines.
#'
#' Methods:
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
#'   units. Aggregation to higher scales takes the pointwise maximum of
#'   the lower-scale maximum functions, without re-weighting, giving a
#'   non-negative additive partition at every scale.
#' * `method = "hill"`: not yet implemented (planned; Hill-Chao framework).
#'
#' **Unbalanced designs.** At every grouping scale, alpha is the mean of
#' the group values weighted by the number of assemblages each group
#' contains (the sample-size weighting of Crist et al. 2003). This
#' matches how gamma pools the data (all assemblages count equally), so
#' beta components cannot turn negative through imbalance alone. In
#' balanced designs it reduces to the plain mean.
#'
#' @param space An `fspace`; `"tpd_eqv"` (and `"rao"` with overlap
#'   dissimilarities) require TPDs (see [fs_tpd()]).
#' @param comm Community matrix (assemblages x units) or long data.frame;
#'   see [fs_structure()]. Mutually exclusive with `levels`.
#' @param hierarchy An [fs_hierarchy()] with one row per assemblage (rows
#'   in the order of `comm`), or `NULL` for a two-level partition
#'   (assemblages vs the whole dataset). Only with `comm`.
#' @param levels Character vector of stored TPD level names, fine to
#'   coarse (at least two: base TPDs, then the assemblage scale); see
#'   Details. Mutually exclusive with `comm`/`hierarchy`.
#' @param method Partition method; see Details.
#' @param q Dominance parameter for `"tpd_eqv"` (default 1; `q = 0`
#'   ignores abundances).
#' @param profile Logical; for `"tpd_eqv"`, also compute the full profile
#'   over `q = seq(0, 2, 0.1)`.
#'
#' @return A list of class `fpartition`: `method`, `levels` (the scale
#'   names), `values` (per-assemblage and per-group values, and gamma),
#'   `table` (the additive components), and `profiles` (when requested).
#' @references de Bello, F. et al. (2010) The partitioning of diversity:
#'   showing Theseus a way out of the labyrinth. *Journal of Vegetation
#'   Science*, 21, 992-1000. Crist, T.O., Veech, J.A., Gering, J.C. &
#'   Summerville, K.S. (2003) Partitioning species diversity across
#'   landscapes and regions: a hierarchical analysis of alpha, beta, and
#'   gamma diversity. *The American Naturalist*, 162, 734-743.
#' @seealso [fs_hierarchy()], [fs_aggregate()], [fs_structure()],
#'   [fs_tpd()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 2, 4), shape = c(0, 1, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)
#' comm <- rbind(A1 = c(sp1 = 1, sp2 = 0, sp3 = 0),
#'               A2 = c(sp1 = 0, sp2 = 1, sp3 = 0),
#'               A3 = c(sp1 = 0, sp2 = 0, sp3 = 1),
#'               A4 = c(sp1 = 1, sp2 = 1, sp3 = 1))
#' h <- fs_hierarchy(site = c("s1", "s1", "s2", "s2"))
#'
#' fs_partition(tp, comm, hierarchy = h, method = "tpd_eqv", q = 0)
#' fs_partition(tp, comm, hierarchy = h, method = "rao")
#'
#' # stored-level route: identical scales, built once with fs_aggregate()
#' tp <- fs_aggregate(tp, comm, name = "plot")
#' tp <- fs_aggregate(tp, c(A1 = "s1", A2 = "s1", A3 = "s2", A4 = "s2"),
#'                    name = "site")
#' fs_partition(tp, levels = c("unit", "plot", "site"), method = "rao")
#'
#' # full q-profile of the equivalent numbers
#' pt <- fs_partition(tp, comm, method = "tpd_eqv", profile = TRUE)
#' plot(pt)
#' @export
fs_partition <- function(space, comm = NULL, hierarchy = NULL,
                         levels = NULL,
                         method = c("rao", "hill", "tpd_eqv"),
                         q = 1, profile = FALSE) {
  stopifnot(is_fspace(space))
  method <- match.arg(method)
  if (method == "hill") {
    stop("method = \"hill\" is not yet implemented; use \"rao\" or ",
         "\"tpd_eqv\".", call. = FALSE)
  }
  if (!is.null(levels)) {
    if (!is.null(comm) || !is.null(hierarchy)) {
      stop("Supply either `comm` (+ optional `hierarchy`) or `levels`, ",
           "not both.", call. = FALSE)
    }
    sc <- .scales_from_levels(space, levels)
  } else {
    if (is.null(comm)) {
      stop("Supply `comm` (+ optional `hierarchy`) or `levels` (stored ",
           "TPD levels).", call. = FALSE)
    }
    if (method == "tpd_eqv" && is.null(space$tpds)) {
      stop("method = \"tpd_eqv\" needs TPDs: run fs_tpd() first.",
           call. = FALSE)
    }
    sc <- .scales_from_comm(space, comm, hierarchy)
  }

  out <- switch(method,
    rao     = .partition_rao(space, sc),
    tpd_eqv = .partition_eqv(space, sc, q,
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

# scale descriptors -------------------------------------------------------------
#
# Both input routes are reduced to one common descriptor `sc`, so the
# partition engines have a single code path:
#   $W           assemblages x base-units effective weight matrix (rows sum 1)
#   $units_level name of the base TPD level (NULL: no TPDs, points fallback)
#   $asm_name    display name of the assemblage scale
#   $groupings   named list of grouping vectors (one entry per assemblage),
#                ordered fine to coarse

.scales_from_comm <- function(space, comm, hierarchy) {
  has_tpds <- !is.null(space$tpds)
  unit_names <- if (has_tpds) names(.unit_tpds(space)) else
    rownames(space$coords)
  W <- .as_comm(comm, unit_names, relative = TRUE)
  H <- .check_hierarchy(hierarchy, nrow(W))
  groupings <- if (is.null(H)) {
    stats::setNames(list(), character(0))
  } else {
    as.list(H)
  }
  list(W = W, units_level = if (has_tpds) "unit" else NULL,
       asm_name = "assemblage", groupings = groupings)
}

# all levels on the aggregation chain of `level`, from itself down
.level_chain <- function(space, level) {
  lv <- space$tpds$levels
  chain <- level
  cur <- lv[[level]]$from
  while (!is.null(cur)) {
    chain <- c(chain, cur)
    cur <- lv[[cur]]$from
  }
  chain
}

.scales_from_levels <- function(space, levels) {
  if (is.null(space$tpds) || is.null(space$tpds$levels)) {
    stop("`levels` needs stored TPDs: run fs_tpd() (and fs_aggregate()) ",
         "first.", call. = FALSE)
  }
  if (!is.character(levels) || length(levels) < 2L) {
    stop("`levels` must be a character vector of at least two stored ",
         "level names, fine to coarse (base TPDs first, then the ",
         "assemblage scale).", call. = FALSE)
  }
  if (anyDuplicated(levels)) {
    stop("`levels` contains duplicated names.", call. = FALSE)
  }
  if ("total" %in% levels[-1L]) {
    stop("A scale cannot be named \"total\": that name is reserved for ",
         "the pooled gamma scale.", call. = FALSE)
  }
  for (l in levels) .check_level(space, l)
  for (k in seq(2L, length(levels))) {
    chain <- .level_chain(space, levels[k])
    if (!levels[k - 1L] %in% chain) {
      stop("Level '", levels[k], "' is not aggregated (directly or ",
           "indirectly) from '", levels[k - 1L], "'; `levels` must run ",
           "fine to coarse along the stored aggregation chain.",
           call. = FALSE)
    }
  }
  base <- levels[1L]
  asm <- levels[2L]
  W <- fs_level_weights(space, asm, to = base)

  groupings <- stats::setNames(list(), character(0))
  if (length(levels) > 2L) {
    for (k in seq(3L, length(levels))) {
      A <- fs_level_weights(space, levels[k], to = asm)
      n_parents <- colSums(A > 0)
      if (any(n_parents != 1L)) {
        stop("Level '", levels[k], "' is not a strict nesting of '", asm,
             "' (some groups share children). Scales above the ",
             "assemblage scale must be strictly nested to partition ",
             "diversity.", call. = FALSE)
      }
      groupings[[levels[k]]] <- rownames(A)[apply(A > 0, 2L, which.max)]
    }
  }
  list(W = W, units_level = base, asm_name = asm, groupings = groupings)
}

# TPD equivalent numbers -------------------------------------------------------

.partition_eqv <- function(space, sc, q, q_grid = NULL) {
  units <- space$tpds$levels[[sc$units_level]]$tpds
  main <- .eqv_levels(units, sc, q)
  profiles <- NULL
  if (!is.null(q_grid)) {
    prof <- lapply(q_grid, function(qq) {
      lv <- .eqv_levels(units, sc, qq)
      data.frame(q = qq, scale = names(lv$alpha), value = unname(lv$alpha))
    })
    profiles <- do.call(rbind, prof)
  }
  list(q = q, levels = names(main$alpha),
       values = main$values,
       table = main$table, profiles = profiles)
}

.eqv_levels <- function(units, sc, q) {
  W <- sc$W
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

  scale_names <- c(sc$asm_name, names(sc$groupings), "total")
  alpha <- numeric(length(scale_names))
  names(alpha) <- scale_names
  alpha[sc$asm_name] <- mean(eqv_a)
  values <- list()
  values[[sc$asm_name]] <- eqv_a

  for (g_name in names(sc$groupings)) {
    gr <- split(seq_len(nrow(W)), sc$groupings[[g_name]])
    g_eqv <- vapply(gr, function(idx) {
      sum(.sparse_max(maxfuns[idx])$vals)
    }, numeric(1L))
    # size-weighted alpha: groups count in proportion to the number of
    # assemblages they contain (Crist et al. 2003), matching gamma
    alpha[g_name] <- stats::weighted.mean(g_eqv, lengths(gr))
    values[[g_name]] <- g_eqv
  }
  gamma <- sum(.sparse_max(maxfuns)$vals)
  alpha["total"] <- gamma
  values$total <- gamma

  comp <- data.frame(
    component = c(paste0("alpha_", scale_names[1L]),
                  paste0("beta_", scale_names[-length(scale_names)],
                         ":", scale_names[-1L])),
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

.partition_rao <- function(space, sc) {
  W <- sc$W
  if (!is.null(sc$units_level)) {
    D <- .overlap_dissim(space$tpds$levels[[sc$units_level]]$tpds)
  } else {
    D <- as.matrix(stats::dist(space$coords))
    D <- D / max(D)
  }
  rao_of <- function(w) as.numeric(t(w) %*% D %*% w)
  eqv_of <- function(Q) 1 / (1 - pmin(Q, 1 - 1e-12))

  rao_a <- apply(W, 1L, rao_of)
  scale_names <- c(sc$asm_name, names(sc$groupings), "total")
  alpha <- numeric(length(scale_names))
  names(alpha) <- scale_names
  alpha[sc$asm_name] <- mean(rao_a)
  values <- list()
  values[[sc$asm_name]] <- rao_a

  for (g_name in names(sc$groupings)) {
    gr <- split(seq_len(nrow(W)), sc$groupings[[g_name]])
    g_rao <- vapply(gr, function(idx) {
      rao_of(colMeans(W[idx, , drop = FALSE]))
    }, numeric(1L))
    # size-weighted alpha (Crist et al. 2003): consistent with gamma's
    # colMeans pooling over all assemblages
    alpha[g_name] <- stats::weighted.mean(g_rao, lengths(gr))
    values[[g_name]] <- g_rao
  }
  gamma <- rao_of(colMeans(W))
  alpha["total"] <- gamma
  values$total <- gamma

  E <- eqv_of(alpha)
  comp <- data.frame(
    component = c(paste0("alpha_", scale_names[1L]),
                  paste0("beta_", scale_names[-length(scale_names)],
                         ":", scale_names[-1L])),
    value = c(alpha[1L], diff(alpha)),
    value_eqv = c(E[1L], diff(E)),
    prop_eqv = c(E[1L], diff(E)) / E[length(E)])

  list(q = NULL, levels = scale_names, values = values,
       table = comp, profiles = NULL)
}
