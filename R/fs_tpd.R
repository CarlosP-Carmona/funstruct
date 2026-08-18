# fs_tpd(): trait probability densities on a grid -----------------------------

#' Define the evaluation grid of a trait space
#'
#' Builds the regular grid on which trait probability densities (TPDs) are
#' evaluated. Default resolution per dimension decreases with
#' dimensionality (200, 50, 25, 15, 8 cells per axis for 1-5 dimensions),
#' keeping the total number of cells manageable.
#'
#' @param space An `fspace` object (typically already reduced with
#'   [fs_reduce()]).
#' @param res Cells per axis (single integer). Default depends on the
#'   dimensionality; see Description.
#' @param padding Fraction of each axis range added on both sides, so that
#'   densities are not clipped at the data limits.
#'
#' @return A list of class `fs_grid` describing the grid (cell midpoints,
#'   steps, cell volume, and the full midpoint matrix in `$cells`).
#' @seealso [fs_tpd()]
#' @export
fs_grid <- function(space, res = NULL, padding = 0.05) {
  stopifnot(is_fspace(space))
  co <- space$coords
  d <- ncol(co)
  .check_tpd_dims(d)
  if (is.null(res)) res <- c(200L, 50L, 25L, 15L, 8L)[d]
  res <- as.integer(res)
  if (length(res) != 1L || is.na(res) || res < 2L) {
    stop("`res` must be a single integer >= 2.", call. = FALSE)
  }
  rng <- apply(co, 2L, range)
  pad <- (rng[2L, ] - rng[1L, ]) * padding
  lo <- rng[1L, ] - pad
  hi <- rng[2L, ] + pad
  step <- (hi - lo) / res
  mids <- lapply(seq_len(d), function(j) {
    lo[j] + step[j] * (seq_len(res) - 0.5)
  })
  cells <- as.matrix(expand.grid(mids, KEEP.OUT.ATTRS = FALSE))
  colnames(cells) <- colnames(co)
  out <- list(d = d, res = res, lo = lo, hi = hi, step = step,
              mids = mids, cells = cells,
              cell_volume = prod(step), n_cells = nrow(cells))
  class(out) <- "fs_grid"
  out
}

#' Estimate trait probability densities
#'
#' Estimates one TPD per unit on the grid of the space, following the TPD
#' framework (Carmona et al. 2016, 2019). Two input routes are supported:
#'
#' **Own observations** (`ids` supplied): the rows of the space are
#' individual observations and `ids` assigns each row to a unit (e.g.
#' individuals to species). Each unit's TPD is a multivariate kernel
#' density estimate of its own observations, with the bandwidth selected
#' per unit from its full data (plug-in by default) and never re-selected
#' on subsets. Alternatively, `obs` can supply new individual-level trait
#' values to be projected into a PCA space built on other units.
#'
#' **Means with imposed spread** (`sds` supplied, no `ids`): the rows of
#' the space are unit means (e.g. species mean traits) and `sds` gives the
#' standard deviation around each mean on each axis (the TPDsMean logic).
#' The kernel width here represents intraspecific spread, a property of
#' the unit: it is attached to the unit and identical wherever the unit
#' appears, never derived from assemblage composition.
#'
#' Bandwidth handling follows a strict rule: the bandwidth is attached to
#' the entity whose spread it represents and stays constant wherever that
#' entity appears. All selected or imposed bandwidths are stored in the
#' returned space and displayed by `summary()`.
#'
#' @param space An `fspace`, typically reduced to few dimensions (see
#'   [fs_adequacy()] for how sample sizes limit dimensionality).
#' @param ids Vector assigning each row of the space (or of `obs`) to a
#'   unit. Triggers the own-observations route.
#' @param obs Optional data.frame of individual observations (same trait
#'   columns as the space's traits) to be projected into a PCA space;
#'   requires `ids` with one entry per row of `obs`.
#' @param sds Standard deviations for the means route: a single value
#'   (all units, all axes), a vector of length `d` (per axis), or a
#'   units x d matrix.
#' @param bw Bandwidth override: `"unit"` (per-unit selection; default for
#'   the observations route), `"common"` (one bandwidth for all units:
#'   the element-wise average of the per-unit selections), a numeric
#'   vector of per-axis standard deviations, a d x d matrix (used directly
#'   as the kernel covariance), or a named list of per-unit matrices.
#' @param bw_selector Bandwidth selector for the observations route:
#'   `"plugin"` (default; `ks::Hpi`/`ks::hpi`) or `"silverman"`
#'   (normal-scale rule; `ks::Hns`/[stats::bw.nrd0()]).
#' @param grid An `fs_grid` object; default `fs_grid(space)`.
#' @param min_n Units with fewer observations than this receive an imputed
#'   bandwidth (the average of the bandwidths selected for well-sampled
#'   units) with a warning.
#' @param alpha Probability mass retained per unit: cells are kept in
#'   decreasing density order until `alpha` is reached, then probabilities
#'   are renormalized (TPD convention; default 0.99). `alpha = 1` keeps
#'   the full grid.
#'
#' @return The `fspace` with two new elements: `tpds` (grid, alpha, and a
#'   sparse per-unit representation: kept cell indices and renormalized
#'   probabilities summing to 1) and `bw` (attachment, selector, per-unit
#'   bandwidths, imputation flags).
#' @references Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J.
#'   (2016) Traits without borders: integrating functional diversity
#'   across scales. *Trends in Ecology & Evolution*, 31, 382-394.
#'   Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J. (2019) Trait
#'   probability density (TPD): measuring functional diversity across
#'   scales based on TPD with R. *Ecology*, 100, e02876.
#' @seealso [fs_grid()], [fs_adequacy()]
#' @export
fs_tpd <- function(space, ids = NULL, obs = NULL, sds = NULL,
                   bw = NULL, bw_selector = c("plugin", "silverman"),
                   grid = NULL, min_n = 5L, alpha = 0.99) {
  stopifnot(is_fspace(space))
  bw_selector <- match.arg(bw_selector)
  d <- ncol(space$coords)
  .check_tpd_dims(d)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha > 1) {
    stop("`alpha` must be a single value in (0, 1].", call. = FALSE)
  }
  if (is.null(grid)) grid <- fs_grid(space)
  if (!inherits(grid, "fs_grid") || grid$d != d) {
    stop("`grid` must be an fs_grid built on a space with the same ",
         "dimensionality.", call. = FALSE)
  }

  if (!is.null(ids)) {
    est <- .tpd_from_obs(space, ids, obs, bw, bw_selector, min_n, d)
  } else if (!is.null(sds) || !is.null(bw)) {
    est <- .tpd_from_means(space, sds, bw, d)
  } else {
    stop("Supply `ids` (rows are observations grouped into units) or ",
         "`sds`/`bw` (rows are unit means with imposed spread). ",
         "See ?fs_tpd.", call. = FALSE)
  }

  units_tpd <- vector("list", length(est$X))
  names(units_tpd) <- names(est$X)
  for (u in names(est$X)) {
    dens <- .kde_at_cells(grid$cells, est$X[[u]], est$H[[u]])
    p <- dens * grid$cell_volume
    s <- sum(p)
    if (s <= 0) {
      stop("Unit '", u, "' has zero probability mass on the grid; ",
           "its observations may lie outside the space. Check inputs.",
           call. = FALSE)
    }
    p <- p / s
    units_tpd[[u]] <- .alpha_trim(p, alpha)
  }

  space$tpds <- list(grid = grid, alpha = alpha, units = units_tpd,
                     n_obs = vapply(est$X, nrow, integer(1L)),
                     route = est$route)
  space$bw <- list(attachment = est$attachment, selector = est$selector,
                   values = est$H, imputed = est$imputed)
  space
}

# routes ----------------------------------------------------------------------

.tpd_from_obs <- function(space, ids, obs, bw, bw_selector, min_n, d) {
  ids <- as.character(ids)
  if (!is.null(obs)) {
    if (space$method != "pca" || is.null(space$traits)) {
      stop("Projecting external observations requires a PCA space with ",
           "stored traits. Alternatively, build the space directly on ",
           "the individual observations and use `ids` alone.",
           call. = FALSE)
    }
    if (length(ids) != nrow(obs)) {
      stop("`ids` must have one entry per row of `obs`.", call. = FALSE)
    }
    X_all <- .project_obs(space, obs)
  } else {
    if (length(ids) != nrow(space$coords)) {
      stop("`ids` must have one entry per row of the space.", call. = FALSE)
    }
    X_all <- space$coords
  }
  groups <- split(seq_len(nrow(X_all)), ids)
  X <- lapply(groups, function(i) X_all[i, , drop = FALSE])
  n_u <- vapply(X, nrow, integer(1L))
  .warn_adequacy(n_u, d)

  # user-supplied bandwidths take precedence
  user <- .parse_user_bw(bw, names(X), d)
  if (!is.null(user)) {
    return(list(X = X, H = user$H, imputed = rep(FALSE, length(X)),
                attachment = user$attachment, selector = "user",
                route = "obs"))
  }

  H <- vector("list", length(X))
  names(H) <- names(X)
  ok <- n_u >= min_n
  for (u in names(X)[ok]) {
    H[[u]] <- .select_bw(X[[u]], bw_selector)
    if (is.null(H[[u]])) ok[u] <- FALSE  # selector failed; impute below
  }
  if (!any(ok)) {
    stop("No unit has enough observations (min_n = ", min_n, ") for ",
         "bandwidth selection. Supply bandwidths via `bw` or use the ",
         "means route with `sds`.", call. = FALSE)
  }
  imputed <- !ok
  if (any(imputed)) {
    H_mean <- Reduce(`+`, H[ok]) / sum(ok)
    for (u in names(X)[imputed]) H[[u]] <- H_mean
    warning(sum(imputed), " unit(s) with fewer than ", min_n,
            " observations received the average bandwidth of the ",
            "well-sampled units: ",
            paste(utils::head(names(X)[imputed], 5L), collapse = ", "),
            if (sum(imputed) > 5L) ", ...", call. = FALSE)
  }
  if (identical(bw, "common")) {
    H_common <- Reduce(`+`, H) / length(H)
    for (u in names(H)) H[[u]] <- H_common
    attachment <- "common"
  } else {
    attachment <- "unit"
  }
  list(X = X, H = H, imputed = imputed, attachment = attachment,
       selector = bw_selector, route = "obs")
}

.tpd_from_means <- function(space, sds, bw, d) {
  co <- space$coords
  units <- rownames(co)
  X <- lapply(seq_len(nrow(co)), function(i) co[i, , drop = FALSE])
  names(X) <- units

  user <- .parse_user_bw(bw, units, d)
  if (!is.null(user)) {
    H <- user$H
    attachment <- user$attachment
  } else {
    if (is.null(sds)) {
      stop("The means route needs `sds` (or an explicit `bw`).",
           call. = FALSE)
    }
    S <- .expand_sds(sds, units, d)
    H <- lapply(units, function(u) diag(S[u, ]^2, nrow = d))
    names(H) <- units
    attachment <- "entity"
  }
  list(X = X, H = H, imputed = rep(FALSE, length(X)),
       attachment = attachment, selector = "user", route = "means")
}

# helpers ---------------------------------------------------------------------

.check_tpd_dims <- function(d) {
  if (d > 5L) {
    stop("TPDs cannot be estimated in more than 5 dimensions: no ",
         "realistic sample size supports density estimation there (see ",
         "fs_adequacy()), and the grid becomes computationally ",
         "intractable. Reduce the space (fs_reduce()) or use the ",
         "point-based engine.", call. = FALSE)
  }
  invisible(d)
}

.project_obs <- function(space, obs) {
  M <- as.matrix(space$traits)
  missing_tr <- setdiff(colnames(M), colnames(obs))
  if (length(missing_tr)) {
    stop("`obs` lacks trait column(s): ",
         paste(missing_tr, collapse = ", "), call. = FALSE)
  }
  O <- as.matrix(as.data.frame(obs)[, colnames(M), drop = FALSE])
  if (anyNA(O)) {
    stop("`obs` contains missing values; impute first.", call. = FALSE)
  }
  scl <- !identical(space$scale, FALSE)
  ctr <- colMeans(M)
  sc <- if (scl) apply(M, 2L, stats::sd) else rep(1, ncol(M))
  Z <- sweep(sweep(O, 2L, ctr), 2L, sc, "/")
  Z %*% space$loadings
}

.select_bw <- function(Xu, selector) {
  d <- ncol(Xu)
  res <- tryCatch({
    if (d == 1L) {
      h <- if (selector == "plugin") ks::hpi(Xu[, 1L]) else
        stats::bw.nrd0(Xu[, 1L])
      matrix(h^2, 1L, 1L)
    } else {
      if (selector == "plugin") ks::Hpi(Xu) else ks::Hns(Xu)
    }
  }, error = function(e) NULL)
  res
}

.parse_user_bw <- function(bw, units, d) {
  if (is.null(bw) || identical(bw, "unit") || identical(bw, "common")) {
    return(NULL)
  }
  if (is.character(bw)) {
    stop("Unknown `bw` option '", bw, "'. Use \"unit\", \"common\", a ",
         "numeric vector of per-axis SDs, a d x d matrix, or a named ",
         "list of per-unit matrices.", call. = FALSE)
  }
  if (is.list(bw)) {
    missing_u <- setdiff(units, names(bw))
    if (length(missing_u)) {
      stop("`bw` list lacks entries for: ",
           paste(utils::head(missing_u, 5L), collapse = ", "),
           call. = FALSE)
    }
    H <- lapply(bw[units], .as_H, d = d)
    names(H) <- units
    return(list(H = H, attachment = "user"))
  }
  H1 <- .as_H(bw, d)
  H <- rep(list(H1), length(units))
  names(H) <- units
  list(H = H, attachment = "user")
}

.as_H <- function(b, d) {
  if (is.matrix(b)) {
    if (nrow(b) != d || ncol(b) != d) {
      stop("Bandwidth matrices must be ", d, " x ", d, ".", call. = FALSE)
    }
    return(b)
  }
  if (is.numeric(b)) {
    if (length(b) == 1L) b <- rep(b, d)
    if (length(b) != d) {
      stop("Numeric bandwidths must have length 1 or ", d,
           " (per-axis standard deviations).", call. = FALSE)
    }
    return(diag(b^2, nrow = d))
  }
  stop("Cannot interpret `bw` entry of class ", class(b)[1L], ".",
       call. = FALSE)
}

.expand_sds <- function(sds, units, d) {
  if (is.matrix(sds) || is.data.frame(sds)) {
    S <- as.matrix(sds)
    if (ncol(S) != d) {
      stop("`sds` must have ", d, " columns (one per dimension).",
           call. = FALSE)
    }
    if (is.null(rownames(S))) {
      if (nrow(S) != length(units)) {
        stop("`sds` must have one row per unit.", call. = FALSE)
      }
      rownames(S) <- units
    }
    missing_u <- setdiff(units, rownames(S))
    if (length(missing_u)) {
      stop("`sds` lacks rows for: ",
           paste(utils::head(missing_u, 5L), collapse = ", "),
           call. = FALSE)
    }
    S <- S[units, , drop = FALSE]
  } else {
    if (length(sds) == 1L) sds <- rep(sds, d)
    if (length(sds) != d) {
      stop("`sds` must be a single value, a vector of length ", d,
           ", or a units x ", d, " matrix.", call. = FALSE)
    }
    S <- matrix(sds, length(units), d, byrow = TRUE,
                dimnames = list(units, NULL))
  }
  if (any(!is.finite(S)) || any(S <= 0)) {
    stop("All `sds` must be finite and positive.", call. = FALSE)
  }
  S
}

.kde_at_cells <- function(G, X, H) {
  d <- ncol(G)
  n <- nrow(X)
  ch <- chol(H)
  Hinv <- chol2inv(ch)
  norm_const <- (2 * pi)^(-d / 2) / prod(diag(ch))
  dens <- numeric(nrow(G))
  for (j in seq_len(n)) {
    dev <- sweep(G, 2L, X[j, ])
    q <- rowSums((dev %*% Hinv) * dev)
    dens <- dens + exp(-0.5 * q)
  }
  dens * norm_const / n
}

.alpha_trim <- function(p, alpha) {
  if (alpha >= 1) {
    keep <- which(p > 0)
    return(list(cells = keep, probs = p[keep] / sum(p[keep])))
  }
  ord <- order(p, decreasing = TRUE)
  cum <- cumsum(p[ord])
  k <- which(cum >= alpha - 1e-12)[1L]
  if (is.na(k)) k <- length(ord)
  keep <- sort(ord[seq_len(k)])
  probs <- p[keep]
  list(cells = keep, probs = probs / sum(probs))
}

.warn_adequacy <- function(n_u, d) {
  req <- if (d <= length(.SILVERMAN_N)) .SILVERMAN_N[d] else Inf
  low <- n_u < req
  if (any(low)) {
    warning(sum(low), " of ", length(n_u), " unit(s) have fewer ",
            "observations than density estimation in ", d,
            " dimension(s) requires (~", format(req, big.mark = ","),
            "; see fs_adequacy()). Interpret their TPDs cautiously or ",
            "reduce the dimensionality.", call. = FALSE)
  }
  invisible(low)
}
