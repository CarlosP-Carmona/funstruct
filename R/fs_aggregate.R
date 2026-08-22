# fs_aggregate(): the level stack -- TPDs at every scale ----------------------
#
# Philosophy (Carmona et al. 2016): a TPD of a species and a TPD of a biome
# are the same kind of object on the same grid. Aggregation is ONE closed
# operation, applied repeatedly (individuals -> species -> communities ->
# habitats -> ... -> biomes). Aggregated TPDs are STORED at every level the
# user asks for, and every index function consumes the same stored stack.
#
# .mix_tpds() is THE single mixture operator of the package: fs_structure(),
# fs_beta(), fs_redundancy() and fs_aggregate() all call it, so there is
# exactly one place where "abundance-weighted mixture of TPDs" is defined.
# Its body reproduces the pre-refactor summation order verbatim; the
# regression tests in test-regression-levels.R pin this.

#' Aggregate stored TPDs into a coarser level
#'
#' Builds one TPD per group as the weighted mixture of the TPDs of a finer
#' level, and stores the result as a new named level in the space's TPD
#' stack. Because aggregation is closed (a mixture of TPDs is a TPD on the
#' same grid), it can be applied repeatedly to build any nesting:
#' individuals to species (done by [fs_tpd()]), species to communities,
#' communities to habitats, and so on up to biomes.
#'
#' **Grouping formats.** `groups` can be (a) a matrix (new groups x
#' children) with raw weight entries and 0 for absent children -- the
#' community-matrix format, needed when a child can belong to several
#' groups -- or (b) a named vector or factor assigning each child to a
#' single parent (`names` are the children). The vector form requires
#' *every* child to be assigned; the matrix form completes missing
#' children with zeros.
#'
#' **Weights.** `weights = "groups"` (default) uses the matrix entries,
#' normalized within each group (abundance behaviour); with the vector
#' form, which carries no weight information, it is equivalent to
#' `"equal"`. `weights = "equal"` gives every present child the same
#' weight. A named numeric vector (e.g. area or sampling effort per
#' child; vector form only) weights each child by its value, normalized
#' within each group. The rule used is recorded in the level object;
#' nothing is silent.
#'
#' **Weight composition.** Across levels, the effective weight of a
#' bottom-level unit in a group is the product of the within-level
#' relative weights along the chain (a mixture of mixtures is the mixture
#' with product weights). [fs_level_weights()] exposes these effective
#' weights for any stored level.
#'
#' **No `alpha` here.** Probability mass trimming happens once, at
#' estimation time in [fs_tpd()]. Aggregation only takes the union of the
#' children's kept cells and renormalizes floating-point drift, so upper
#' levels occupy more cells (reported by `summary()`).
#'
#' @param space An `fspace` with stored TPDs (see [fs_tpd()]).
#' @param groups Grouping of the children of level `from`: a matrix (new
#'   groups x children, children as column names) or a named vector /
#'   factor (child -> parent). See Details.
#' @param name Name of the new level (a single string, not yet in use).
#' @param from Name of the level whose TPDs are aggregated. Default: the
#'   most recently added level.
#' @param weights `"groups"`, `"equal"`, or a named numeric vector (vector
#'   form only). See Details.
#'
#' @return The `fspace` with the new level appended to
#'   `space$tpds$levels`: a list with the group TPDs (`tpds`), the parent
#'   level (`from`), the normalized membership matrix (`members`, rows
#'   summing to 1), the weight rule used (`weights_rule`) and the number
#'   of children per group (`n_children`).
#' @references Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J.
#'   (2016) Traits without borders: integrating functional diversity
#'   across scales. *Trends in Ecology & Evolution*, 31, 382-394.
#' @seealso [fs_level_weights()], [fs_get_tpd()], [fs_tpd()],
#'   [fs_structure()], [fs_beta()], [fs_partition()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 2, 4), shape = c(0, 1, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)                    # level "unit"
#'
#' comm <- rbind(A = c(sp1 = 2, sp2 = 1, sp3 = 0),
#'               B = c(sp1 = 0, sp2 = 1, sp3 = 3))
#' tp <- fs_aggregate(tp, comm, name = "plot")    # matrix form
#'
#' # plots into one region, weighted by plot area
#' tp <- fs_aggregate(tp, c(A = "r1", B = "r1"), name = "region",
#'                    weights = c(A = 12, B = 30))
#' summary(tp)
#' fs_level_weights(tp, "region")
#' @export
fs_aggregate <- function(space, groups, name, from = NULL,
                         weights = "groups") {
  stopifnot(is_fspace(space))
  if (is.null(space$tpds) || is.null(space$tpds$levels)) {
    stop("fs_aggregate() needs stored TPDs: run fs_tpd() first.",
         call. = FALSE)
  }
  lv <- space$tpds$levels
  if (is.null(from)) from <- names(lv)[length(lv)]
  if (!is.character(from) || length(from) != 1L || !from %in% names(lv)) {
    stop("`from` must name a stored level (",
         paste(names(lv), collapse = ", "), ").", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(name)) {
    stop("`name` must be a single non-empty string.", call. = FALSE)
  }
  if (name %in% names(lv)) {
    stop("Level '", name, "' already exists. Levels are appended, not ",
         "replaced; pick a new name (or re-run fs_tpd() to start over).",
         call. = FALSE)
  }
  children <- names(lv[[from]]$tpds)
  M <- .as_members(groups, children, weights)
  rule <- attr(M, "weights_rule")
  attr(M, "weights_rule") <- NULL

  tpds <- vector("list", nrow(M))
  names(tpds) <- rownames(M)
  for (g in rownames(M)) {
    tpds[[g]] <- .mix_tpds(lv[[from]]$tpds, M[g, ])
  }
  space$tpds$levels[[name]] <- list(tpds = tpds, from = from,
                                    members = M, weights_rule = rule,
                                    n_children = rowSums(M > 0))
  space
}

# membership matrix ------------------------------------------------------------

# Turns either grouping format into a groups x children matrix whose rows
# are the relative within-level weights (rows sum to 1). The rule actually
# applied is returned as attr(, "weights_rule").
.as_members <- function(groups, children, weights) {
  numeric_w <- is.numeric(weights)
  if (!numeric_w &&
      (!is.character(weights) || length(weights) != 1L ||
       !weights %in% c("groups", "equal"))) {
    stop("`weights` must be \"groups\", \"equal\", or a named numeric ",
         "vector (vector form only).", call. = FALSE)
  }

  if ((is.matrix(groups) || is.data.frame(groups)) && !is.null(dim(groups))) {
    M <- as.matrix(groups)
    if (!is.numeric(M)) {
      stop("A `groups` matrix must be numeric (weight entries, 0 = ",
           "absent).", call. = FALSE)
    }
    if (is.null(colnames(M))) {
      stop("A `groups` matrix must have the children as column names.",
           call. = FALSE)
    }
    unknown <- setdiff(colnames(M), children)
    if (length(unknown)) {
      stop("Unknown child(ren) in `groups`: ",
           paste(utils::head(unknown, 5L), collapse = ", "),
           if (length(unknown) > 5L) ", ...", call. = FALSE)
    }
    if (is.null(rownames(M))) {
      rownames(M) <- paste0("group", seq_len(nrow(M)))
    }
    if (anyNA(M) || any(M < 0)) {
      stop("Weight entries must be non-negative and complete.",
           call. = FALSE)
    }
    if (numeric_w) {
      stop("Named numeric `weights` apply to the vector form only; the ",
           "matrix entries already carry the weights.", call. = FALSE)
    }
    # complete missing children with zeros, in canonical order
    full <- matrix(0, nrow(M), length(children),
                   dimnames = list(rownames(M), children))
    full[, colnames(M)] <- M
    M <- full
    if (identical(weights, "equal")) {
      M <- (M > 0) * 1
    }
  } else {
    # named vector / factor: child -> single parent
    parents <- as.character(groups)
    ch <- names(groups)
    if (is.null(ch) || any(!nzchar(ch))) {
      stop("The vector form of `groups` must be named (names are the ",
           "children of the level to aggregate).", call. = FALSE)
    }
    if (anyDuplicated(ch)) {
      stop("Duplicated child name(s) in `groups`: each child has a ",
           "single parent in the vector form (use the matrix form for ",
           "partial membership).", call. = FALSE)
    }
    unknown <- setdiff(ch, children)
    if (length(unknown)) {
      stop("Unknown child(ren) in `groups`: ",
           paste(utils::head(unknown, 5L), collapse = ", "),
           if (length(unknown) > 5L) ", ...", call. = FALSE)
    }
    unassigned <- setdiff(children, ch)
    if (length(unassigned)) {
      stop("The vector form must assign every child; missing: ",
           paste(utils::head(unassigned, 5L), collapse = ", "),
           if (length(unassigned) > 5L) ", ...",
           " (use the matrix form for partial membership).", call. = FALSE)
    }
    if (anyNA(parents)) {
      stop("`groups` contains missing assignments.", call. = FALSE)
    }
    g_names <- if (is.factor(groups)) {
      levels(droplevels(groups))
    } else {
      unique(parents)
    }
    M <- matrix(0, length(g_names), length(children),
                dimnames = list(g_names, children))
    w_child <- if (numeric_w) {
      wn <- names(weights)
      if (is.null(wn)) {
        stop("Numeric `weights` must be named by child.", call. = FALSE)
      }
      missing_w <- setdiff(children, wn)
      if (length(missing_w)) {
        stop("`weights` lacks entries for: ",
             paste(utils::head(missing_w, 5L), collapse = ", "),
             if (length(missing_w) > 5L) ", ...", call. = FALSE)
      }
      w <- weights[children]
      if (anyNA(w) || any(!is.finite(w)) || any(w <= 0)) {
        stop("All `weights` must be finite and positive.", call. = FALSE)
      }
      as.numeric(w)
    } else {
      # "groups" carries no information in the vector form: equal weights
      rep(1, length(children))
    }
    names(w_child) <- children
    for (k in seq_along(ch)) {
      M[parents[k], ch[k]] <- w_child[ch[k]]
    }
  }

  rs <- rowSums(M)
  if (any(rs == 0)) {
    stop("Group(s) with no children / zero total weight: ",
         paste(utils::head(rownames(M)[rs == 0], 5L), collapse = ", "),
         call. = FALSE)
  }
  M <- M / rs
  attr(M, "weights_rule") <- if (numeric_w) {
    "user"
  } else if (identical(weights, "equal") ||
             (!is.matrix(groups) && !is.data.frame(groups))) {
    "equal"
  } else {
    "groups"
  }
  M
}

# THE mixture operator ---------------------------------------------------------

# Weighted mixture of sparse TPDs ({cells, probs} lists) on a shared grid.
# `w` is a weight vector aligned with `units` (0 = absent). The body
# reproduces the pre-refactor summation order verbatim (see
# helper-legacy.R); do not "improve" it, the regression pins depend on the
# exact floating-point path. With `counts = TRUE` the number of
# contributing units per cell is returned too (used by redundancy).
.mix_tpds <- function(units, w, counts = FALSE) {
  present <- which(w > 0)
  parts <- lapply(present, function(u) {
    list(cells = units[[u]]$cells, p = units[[u]]$probs * w[u])
  })
  all_cells <- sort(unique(unlist(lapply(parts, `[[`, "cells"))))
  pc <- numeric(length(all_cells))
  n_cell <- integer(length(all_cells))
  pos <- seq_along(all_cells)
  names(pos) <- all_cells
  for (m in parts) {
    j <- pos[as.character(m$cells)]
    pc[j] <- pc[j] + m$p
    n_cell[j] <- n_cell[j] + 1L
  }
  out <- list(cells = all_cells, probs = pc / sum(pc))
  if (counts) out$counts <- n_cell
  out
}

# level helpers ----------------------------------------------------------------

# The bottom-level (estimation) TPD list, wherever the stack stores it.
.unit_tpds <- function(space) {
  space$tpds$levels$unit$tpds
}

# Validates that `level` names a stored level; returns it invisibly.
.check_level <- function(space, level) {
  if (is.null(space$tpds) || is.null(space$tpds$levels)) {
    stop("This space stores no TPDs: run fs_tpd() first.", call. = FALSE)
  }
  if (!is.character(level) || length(level) != 1L ||
      !level %in% names(space$tpds$levels)) {
    stop("`level` must name a stored level (",
         paste(names(space$tpds$levels), collapse = ", "), ").",
         call. = FALSE)
  }
  invisible(level)
}

#' Effective weights of a stored level
#'
#' Returns, for every group of a stored level, the effective weights of
#' the units of a finer level (`to`, by default the bottom `"unit"`
#' level): the product of the within-level relative weights along the
#' aggregation chain. These are the weights with which the finer TPDs mix
#' into each group's TPD, and the weights [fs_structure()] and
#' [fs_partition()] use internally when working on stored levels.
#'
#' @param space An `fspace` with stored TPDs (see [fs_tpd()],
#'   [fs_aggregate()]).
#' @param level Name of the stored level whose groups are the rows.
#' @param to Name of the finer level whose units are the columns (default
#'   `"unit"`); must lie on the aggregation chain of `level`.
#'
#' @return A matrix (groups of `level` x units of `to`) whose rows sum
#'   to 1. When `level == to`, the identity matrix.
#' @seealso [fs_aggregate()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 2, 4), shape = c(0, 1, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)
#' comm <- rbind(A = c(sp1 = 2, sp2 = 1, sp3 = 0),
#'               B = c(sp1 = 0, sp2 = 1, sp3 = 3))
#' tp <- fs_aggregate(tp, comm, name = "plot")
#' tp <- fs_aggregate(tp, c(A = "r1", B = "r1"), name = "region",
#'                    weights = c(A = 12, B = 30))
#' fs_level_weights(tp, "region")   # region x species, rows sum to 1
#' @export
fs_level_weights <- function(space, level, to = "unit") {
  stopifnot(is_fspace(space))
  .check_level(space, level)
  .check_level(space, to)
  lv <- space$tpds$levels
  if (identical(level, to)) {
    nm <- names(lv[[level]]$tpds)
    W <- diag(1, length(nm))
    dimnames(W) <- list(nm, nm)
    return(W)
  }
  if (is.null(lv[[level]]$members)) {
    stop("'", to, "' is not on the aggregation chain of '", level, "'.",
         call. = FALSE)
  }
  W <- lv[[level]]$members
  cur <- lv[[level]]$from
  while (!identical(cur, to)) {
    if (is.null(cur) || is.null(lv[[cur]]$members)) {
      stop("'", to, "' is not on the aggregation chain of '", level, "'.",
           call. = FALSE)
    }
    W <- W %*% lv[[cur]]$members
    cur <- lv[[cur]]$from
  }
  W
}

# extraction and plotting ------------------------------------------------------

#' Extract stored TPDs
#'
#' Pulls the TPDs of a stored level (all groups, or a selection) out of
#' the space, together with the grid, as a compact object with a plot
#' method (see [plot.ftpd()]).
#'
#' @param space An `fspace` with stored TPDs.
#' @param level Name of a stored level (default `"unit"`).
#' @param ids Optional character vector selecting groups of that level;
#'   default all.
#'
#' @return A list of class `ftpd`: `grid`, `level`, and `tpds` (named
#'   sparse list with `cells` and `probs` per group).
#' @seealso [fs_aggregate()], [fs_tpd()], [plot.ftpd()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 2, 4), shape = c(0, 1, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)
#' comm <- rbind(A = c(sp1 = 2, sp2 = 1, sp3 = 0),
#'               B = c(sp1 = 0, sp2 = 1, sp3 = 3))
#' tp <- fs_aggregate(tp, comm, name = "plot")
#'
#' f <- fs_get_tpd(tp, "plot", "A")
#' f
#' plot(f)
#' @export
fs_get_tpd <- function(space, level = "unit", ids = NULL) {
  stopifnot(is_fspace(space))
  .check_level(space, level)
  tpds <- space$tpds$levels[[level]]$tpds
  if (!is.null(ids)) {
    ids <- as.character(ids)
    unknown <- setdiff(ids, names(tpds))
    if (length(unknown)) {
      stop("Unknown id(s) at level '", level, "': ",
           paste(utils::head(unknown, 5L), collapse = ", "),
           if (length(unknown) > 5L) ", ...", call. = FALSE)
    }
    tpds <- tpds[ids]
  }
  out <- list(grid = space$tpds$grid, level = level, tpds = tpds)
  class(out) <- "ftpd"
  out
}

#' @export
print.ftpd <- function(x, ...) {
  g <- x$grid
  cat("<ftpd> level '", x$level, "': ", length(x$tpds), " TPD(s) on a ",
      g$res, "^", g$d, " grid\n", sep = "")
  occ <- vapply(x$tpds, function(t) length(t$cells), integer(1L))
  show <- utils::head(seq_along(x$tpds), 10L)
  for (i in show) {
    cat("  ", names(x$tpds)[i], ": ", occ[i], " cells (",
        round(100 * occ[i] / g$n_cells, 1), "% of the grid)\n", sep = "")
  }
  if (length(x$tpds) > 10L) {
    cat("  ... and ", length(x$tpds) - 10L, " more\n", sep = "")
  }
  invisible(x)
}

#' Plot a stored TPD
#'
#' Displays one TPD extracted with [fs_get_tpd()]: a density profile in
#' one dimension, or a density image with contours in two. As everywhere
#' in the package, `par()` is not modified, so further elements can be
#' added to the plot afterwards.
#'
#' @param x An `ftpd` object (see [fs_get_tpd()]).
#' @param id Which TPD to plot (a name of `x$tpds`); default the first.
#' @param asp Aspect ratio for two-dimensional plots (default 1: equal
#'   axis scaling, so distances read true).
#' @param ... Further arguments passed to [graphics::plot()] (1D) or
#'   [graphics::image()] (2D).
#'
#' @return Invisibly, the plotted density vector (1D) or matrix (2D).
#' @seealso [fs_get_tpd()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 2, 4), shape = c(0, 1, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)
#' plot(fs_get_tpd(tp, "unit", "sp2"))
#' @export
plot.ftpd <- function(x, id = NULL, asp = 1, ...) {
  g <- x$grid
  if (g$d > 2L) {
    stop("Plotting is available for 1- and 2-dimensional grids only.",
         call. = FALSE)
  }
  if (is.null(id)) id <- names(x$tpds)[1L]
  id <- as.character(id)[1L]
  if (!id %in% names(x$tpds)) {
    stop("'", id, "' is not among the extracted TPDs.", call. = FALSE)
  }
  td <- x$tpds[[id]]
  dens <- numeric(g$n_cells)
  dens[td$cells] <- td$probs / g$cell_volume
  lab <- colnames(g$cells)

  if (g$d == 1L) {
    graphics::plot(g$mids[[1L]], dens, type = "l", lwd = 2,
                   col = "#08306B", xlab = lab[1L], ylab = "Density",
                   las = 1, ...)
    return(invisible(dens))
  }
  z <- matrix(dens, g$res, g$res)
  z[z == 0] <- NA
  # fixed hex ramp: named hcl palettes have proven platform-dependent
  ramp <- grDevices::colorRampPalette(
    c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(100)
  graphics::image(g$mids[[1L]], g$mids[[2L]], z, col = ramp,
                  xlab = lab[1L], ylab = lab[2L], las = 1, asp = asp, ...)
  graphics::contour(g$mids[[1L]], g$mids[[2L]], z, add = TRUE,
                    col = "grey40", drawlabels = FALSE, nlevels = 5)
  invisible(z)
}
