# fs_dist(): Gower-style functional dissimilarities, all bounded 0-1 ----------

#' Functional dissimilarities between units (bounded 0-1)
#'
#' Computes pairwise functional dissimilarities from multiple types of
#' trait information, combining them Gower-style: every trait contributes
#' a distance bounded between 0 and 1, and the combined dissimilarity is
#' their (weighted) mean, so it is itself bounded between 0 and 1. This is
#' the recommended input for [fs_space()] with `method = "pcoa"` or
#' `"nmds"`, and the bounded scale carries through to
#' dissimilarity-based indices (Rao, MPD): a Rao value of 0.25 means "the
#' expected pair of individuals differs in a quarter of the trait range",
#' regardless of how many traits were used.
#'
#' Beyond the classical Gower ingredients (quantitative traits as range-
#' scaled differences, categorical traits as 0/1 matching), `fs_dist` can
#' represent intraspecific variability directly in the distances,
#' following Carmona et al. (2015): when two units overlap in their trait
#' distributions, their distance is the proportion of non-overlap rather
#' than the difference between their means.
#'
#' @section Trait types:
#'
#' * **Quantitative, mean only** (numeric columns of `traits` without an
#'   entry in `sds` or `obs`): classical Gower, `|x_i - x_j| / range`.
#' * **Categorical** (factor/character/logical columns of `traits`):
#'   0 if the two units share the level, 1 otherwise.
#' * **Quantitative with known spread** (numeric columns of `traits` with
#'   a matching column in `sds`): each unit is a normal distribution
#'   (mean from `traits`, sd from `sds`); the distance is 1 minus the
#'   overlap of the two densities. Zero or missing sds are replaced by
#'   the mean sd of that trait, with a warning.
#' * **Quantitative with raw observations** (columns of `obs`, grouped by
#'   `ids`): each unit is a Gaussian kernel density estimate of its
#'   observations (bandwidth `1.06 * sd * n^-0.2`); the distance is
#'   1 minus the overlap of the two densities. Units with fewer than
#'   `min_obs` observations fall back to a normal with the across-unit
#'   mean sd, with a warning.
#' * **Category proportions** (`proportions`): units described by the
#'   proportion of individuals in each level of a trait; the distance is
#'   1 minus the summed pairwise minima of the proportions.
#' * **Periods** (`periods`): units described by a start and an end on a
#'   cycle (e.g. flowering period in days of the year, `cycle = 365`,
#'   wrap-around allowed); the distance is 1 minus the shared fraction of
#'   the shorter period.
#'
#' All types can be mixed freely; each trait contributes one 0-1 distance
#' matrix and the combined dissimilarity is their weighted mean. Missing
#' values in the numeric/categorical columns of `traits` are tolerated:
#' that trait is simply skipped for pairs where either unit is `NA`, and
#' the mean is taken over the traits available to each pair (classical
#' Gower behaviour).
#'
#' @param traits Optional units x traits matrix or data.frame with unique
#'   row names. Numeric columns are treated as quantitative traits,
#'   factor/character/logical columns as categorical.
#' @param sds Optional matrix or data.frame of standard deviations for a
#'   subset of the numeric columns of `traits` (same row names, column
#'   names matching the trait names). Switches those traits to the
#'   normal-overlap distance.
#' @param obs Optional data.frame of raw observations (one row per
#'   observation, numeric columns; column names must not repeat names
#'   already used in `traits`). Requires `ids`.
#' @param ids Vector (length `nrow(obs)`) assigning each observation to a
#'   unit.
#' @param proportions Optional named list; each element is a units x
#'   categories matrix of proportions (rows summing to 1) for one trait.
#' @param periods Optional named list; each element is a units x 2 matrix
#'   (start, end) on `1..cycle` for one trait; `start > end` wraps around
#'   the cycle.
#' @param cycle Length of the cycle for `periods` (default 365 days).
#' @param weights Optional named numeric vector of trait weights (default:
#'   equal weights). Names must match the trait names across all inputs.
#' @param min_obs Minimum observations per unit for the kernel route
#'   (default 4); units below it fall back to the normal-overlap route.
#'
#' @return A [stats::dist] object bounded in `[0, 1]`, with attributes
#'   `per_trait` (named list of the per-trait 0-1 `dist` objects),
#'   `types` (named character vector of trait types) and `weights`.
#' @references Carmona, C.P., Rota, C., Azcárate, F.M. & Peco, B. (2015)
#'   More for less: sampling strategies of plant functional traits across
#'   local environmental gradients. *Functional Ecology*, 29, 579-588.
#' @seealso [fs_space()] (use the result as its `traits` input for
#'   `method = "pcoa"` or `"nmds"`).
#' @examples
#' data(gspff)
#' d <- fs_dist(gspff[1:80, ])          # range-scaled Gower, all numeric
#' range(d)                             # bounded 0-1
#' sp <- fs_space(d, method = "pcoa")   # the intended workflow
#'
#' # intraspecific variability from raw observations:
#' data(grassland)
#' tr <- data.frame(height = log10(grassland$height),
#'                  sla = log10(grassland$sla))
#' d2 <- fs_dist(obs = tr, ids = grassland$species)
#' attr(d2, "types")
#' @export
fs_dist <- function(traits = NULL, sds = NULL, obs = NULL, ids = NULL,
                    proportions = NULL, periods = NULL, cycle = 365,
                    weights = NULL, min_obs = 4L) {
  if (!is.null(sds) && is.null(traits)) {
    stop("`sds` requires `traits` (the corresponding means).",
         call. = FALSE)
  }
  if (is.null(traits) && is.null(obs) && is.null(proportions) &&
      is.null(periods)) {
    stop("Provide at least one of `traits`, `obs` (+ `ids`), ",
         "`proportions` or `periods`.", call. = FALSE)
  }
  if (!is.null(obs) && is.null(ids)) {
    stop("`obs` requires `ids` assigning each observation to a unit.",
         call. = FALSE)
  }

  # -- establish the common set of units --------------------------------------
  units <- NULL
  if (!is.null(traits)) {
    traits <- as.data.frame(traits)
    if (is.null(rownames(traits)) ||
        anyDuplicated(rownames(traits)) > 0L ||
        identical(rownames(traits), as.character(seq_len(nrow(traits))))) {
      stop("`traits` must have unique, informative row names.",
           call. = FALSE)
    }
    units <- rownames(traits)
  }
  if (!is.null(obs)) {
    ids <- as.character(ids)
    if (length(ids) != nrow(as.data.frame(obs))) {
      stop("`ids` must have one entry per row of `obs`.", call. = FALSE)
    }
    if (is.null(units)) units <- unique(ids)
    else if (!setequal(unique(ids), units)) {
      stop("The units in `ids` do not match the row names of `traits`.",
           call. = FALSE)
    }
  }
  for (nm in c("proportions", "periods")) {
    lst <- get(nm)
    if (is.null(lst)) next
    if (!is.list(lst) || is.null(names(lst)) || any(!nzchar(names(lst)))) {
      stop("`", nm, "` must be a named list (one element per trait).",
           call. = FALSE)
    }
    for (tn in names(lst)) {
      m <- as.matrix(lst[[tn]])
      if (is.null(rownames(m))) {
        stop("Element '", tn, "' of `", nm, "` needs row names.",
             call. = FALSE)
      }
      if (is.null(units)) units <- rownames(m)
      else if (!setequal(rownames(m), units)) {
        stop("Element '", tn, "' of `", nm, "` does not cover the same ",
             "units as the other inputs.", call. = FALSE)
      }
    }
  }
  n <- length(units)
  if (n < 3L) stop("At least three units are required.", call. = FALSE)

  # -- per-trait distance matrices, each bounded 0-1 --------------------------
  per_trait <- list()
  types <- character()

  if (!is.null(traits)) {
    sd_names <- character()
    if (!is.null(sds)) {
      sds <- as.data.frame(sds)
      sd_names <- colnames(sds)
      bad <- setdiff(sd_names, colnames(traits))
      if (length(bad)) {
        stop("Columns of `sds` not found among `traits`: ",
             paste(bad, collapse = ", "), call. = FALSE)
      }
      if (!setequal(rownames(sds), units)) {
        stop("`sds` must cover the same units as `traits`.", call. = FALSE)
      }
      sds <- sds[units, , drop = FALSE]
    }
    for (tn in colnames(traits)) {
      x <- traits[units, tn]
      if (is.numeric(x) && tn %in% sd_names) {
        s <- as.numeric(sds[[tn]])
        if (anyNA(x)) {
          stop("NA means are not allowed for normal-overlap traits ('",
               tn, "').", call. = FALSE)
        }
        if (anyNA(s) || any(s <= 0, na.rm = TRUE)) {
          fill <- mean(s[!is.na(s) & s > 0])
          if (!is.finite(fill)) {
            stop("All sds are missing or zero for trait '", tn, "'.",
                 call. = FALSE)
          }
          s[is.na(s) | s <= 0] <- fill
          warning("Zero or missing sds for trait '", tn,
                  "' replaced by the mean sd (", signif(fill, 3), ").",
                  call. = FALSE)
        }
        per_trait[[tn]] <- .dist_gaussian(x, s, units)
        types[tn] <- "gaussian"
      } else if (is.numeric(x)) {
        per_trait[[tn]] <- .dist_gower_num(x, units)
        types[tn] <- "quantitative"
      } else {
        per_trait[[tn]] <- .dist_categorical(as.character(x), units)
        types[tn] <- "categorical"
      }
    }
  }

  if (!is.null(obs)) {
    obs <- as.data.frame(obs)
    num_ok <- vapply(obs, is.numeric, logical(1L))
    if (!all(num_ok)) {
      stop("All columns of `obs` must be numeric.", call. = FALSE)
    }
    clash <- intersect(colnames(obs), names(per_trait))
    if (length(clash)) {
      stop("Trait names repeated across inputs: ",
           paste(clash, collapse = ", "), call. = FALSE)
    }
    for (tn in colnames(obs)) {
      per_trait[[tn]] <- .dist_kernel(obs[[tn]], ids, units, min_obs)
      types[tn] <- "kernel"
    }
  }

  if (!is.null(proportions)) {
    for (tn in names(proportions)) {
      if (tn %in% names(per_trait)) {
        stop("Trait names repeated across inputs: ", tn, call. = FALSE)
      }
      m <- as.matrix(proportions[[tn]])[units, , drop = FALSE]
      if (anyNA(m) || any(m < 0) || any(m > 1)) {
        stop("Proportions for '", tn, "' must be in [0, 1] with no NAs.",
             call. = FALSE)
      }
      if (any(abs(rowSums(m) - 1) > 1e-6)) {
        stop("Proportions for '", tn, "' must sum to 1 within each unit.",
             call. = FALSE)
      }
      per_trait[[tn]] <- .dist_proportions(m, units)
      types[tn] <- "proportions"
    }
  }

  if (!is.null(periods)) {
    for (tn in names(periods)) {
      if (tn %in% names(per_trait)) {
        stop("Trait names repeated across inputs: ", tn, call. = FALSE)
      }
      m <- as.matrix(periods[[tn]])[units, , drop = FALSE]
      if (ncol(m) != 2L || anyNA(m) ||
          any(m < 1) || any(m > cycle)) {
        stop("Periods for '", tn, "' must be a units x 2 (start, end) ",
             "matrix with values in 1..", cycle, ".", call. = FALSE)
      }
      per_trait[[tn]] <- .dist_periods(m, units, cycle)
      types[tn] <- "periods"
    }
  }

  # -- weighted Gower combination ---------------------------------------------
  tn_all <- names(per_trait)
  if (is.null(weights)) {
    weights <- stats::setNames(rep(1, length(tn_all)), tn_all)
  } else {
    if (is.null(names(weights)) || !setequal(names(weights), tn_all)) {
      stop("`weights` must be a named vector covering exactly these ",
           "traits: ", paste(tn_all, collapse = ", "), call. = FALSE)
    }
    if (any(weights < 0) || all(weights == 0)) {
      stop("`weights` must be non-negative and not all zero.",
           call. = FALSE)
    }
    weights <- weights[tn_all]
  }

  wsum <- num <- matrix(0, n, n)
  for (tn in tn_all) {
    D <- as.matrix(per_trait[[tn]])
    ok <- !is.na(D)
    num[ok] <- num[ok] + weights[tn] * D[ok]
    wsum[ok] <- wsum[ok] + weights[tn]
  }
  if (any(wsum[upper.tri(wsum)] == 0)) {
    stop("Some pairs of units share no non-missing trait, so their ",
         "dissimilarity is undefined.", call. = FALSE)
  }
  G <- num / wsum
  G <- pmin(pmax(G, 0), 1)
  diag(G) <- 0
  dimnames(G) <- list(units, units)
  out <- stats::as.dist(G)
  attr(out, "per_trait") <- per_trait
  attr(out, "types") <- types
  attr(out, "weights") <- weights
  attr(out, "method") <- "gower"
  attr(out, "call") <- match.call()
  out
}

# per-trait engines (all return `dist`, values in [0, 1] or NA) ---------------

.dist_gower_num <- function(x, units) {
  if (all(is.na(x))) {
    stop("A quantitative trait has no non-missing values.", call. = FALSE)
  }
  rng <- diff(range(x, na.rm = TRUE))
  D <- abs(outer(x, x, "-"))
  if (rng > 0) D <- D / rng else D[] <- 0
  dimnames(D) <- list(units, units)
  stats::as.dist(D)
}

.dist_categorical <- function(x, units) {
  D <- 1 - outer(x, x, "==")  # logical -> numeric; NA propagates
  dimnames(D) <- list(units, units)
  stats::as.dist(D)
}

.dist_gaussian <- function(m, s, units) {
  n <- length(m)
  D <- matrix(0, n, n, dimnames = list(units, units))
  for (j in seq_len(n - 1L)) {
    for (k in (j + 1L):n) {
      D[k, j] <- 1 - .overlap_normals(m[j], s[j], m[k], s[k])
    }
  }
  stats::as.dist(D)
}

.dist_kernel <- function(x, ids, units, min_obs) {
  splits <- split(x, factor(ids, levels = units))
  splits <- lapply(splits, function(v) v[!is.na(v)])
  n_obs <- lengths(splits)
  small <- n_obs < min_obs
  if (any(n_obs == 0L)) {
    stop("Units with no observations: ",
         paste(units[n_obs == 0L], collapse = ", "), call. = FALSE)
  }
  pooled_sd <- mean(vapply(splits[!small & n_obs > 1L], stats::sd,
                           numeric(1L)))
  if (any(small)) {
    if (!is.finite(pooled_sd)) {
      stop("Too few observations overall to apply the kernel route ",
           "(no unit reaches `min_obs`).", call. = FALSE)
    }
    warning(sum(small), " unit(s) with fewer than ", min_obs,
            " observations use a normal with the across-unit mean sd: ",
            paste(units[small], collapse = ", "), call. = FALSE)
  }
  n <- length(units)
  D <- matrix(0, n, n, dimnames = list(units, units))
  for (j in seq_len(n - 1L)) {
    for (k in (j + 1L):n) {
      D[k, j] <- 1 - .overlap_units(splits[[j]], splits[[k]],
                                    small[j], small[k], pooled_sd)
    }
  }
  stats::as.dist(D)
}

.dist_proportions <- function(m, units) {
  n <- nrow(m)
  D <- matrix(0, n, n, dimnames = list(units, units))
  for (j in seq_len(n - 1L)) {
    for (k in (j + 1L):n) {
      D[k, j] <- 1 - sum(pmin(m[j, ], m[k, ]))
    }
  }
  stats::as.dist(D)
}

.dist_periods <- function(m, units, cycle) {
  expand <- function(a, b) {
    if (a <= b) seq.int(a, b) else c(seq_len(b), seq.int(a, cycle))
  }
  days <- lapply(seq_len(nrow(m)), function(i) expand(m[i, 1L], m[i, 2L]))
  n <- nrow(m)
  D <- matrix(0, n, n, dimnames = list(units, units))
  for (j in seq_len(n - 1L)) {
    for (k in (j + 1L):n) {
      shared <- sum(days[[j]] %in% days[[k]])
      D[k, j] <- 1 - shared / min(length(days[[j]]), length(days[[k]]))
    }
  }
  stats::as.dist(D)
}

# overlap of two densities on a shared grid -----------------------------------

.overlap_grid <- function(f1, f2, dx) {
  o <- pmin(f1, f2)
  # trapezoidal rule
  ov <- dx * (sum(o) - (o[1L] + o[length(o)]) / 2)
  min(max(ov, 0), 1)
}

.overlap_normals <- function(m1, s1, m2, s2, n_grid = 512L) {
  r <- range(m1 - 5 * s1, m1 + 5 * s1, m2 - 5 * s2, m2 + 5 * s2)
  g <- seq(r[1L], r[2L], length.out = n_grid)
  .overlap_grid(stats::dnorm(g, m1, s1), stats::dnorm(g, m2, s2),
                g[2L] - g[1L])
}

.overlap_units <- function(v1, v2, small1, small2, pooled_sd,
                           n_grid = 512L) {
  m1 <- mean(v1); m2 <- mean(v2)
  s1 <- if (small1) pooled_sd else stats::sd(v1)
  s2 <- if (small2) pooled_sd else stats::sd(v2)
  r <- range(m1 - 5 * s1, m1 + 5 * s1, m2 - 5 * s2, m2 + 5 * s2)
  g <- seq(r[1L], r[2L], length.out = n_grid)
  f1 <- if (small1) stats::dnorm(g, m1, s1) else {
    h <- 1.06 * s1 * length(v1)^-0.2
    stats::density(v1, bw = h, kernel = "gaussian",
                   from = r[1L], to = r[2L], n = n_grid)$y
  }
  f2 <- if (small2) stats::dnorm(g, m2, s2) else {
    h <- 1.06 * s2 * length(v2)^-0.2
    stats::density(v2, bw = h, kernel = "gaussian",
                   from = r[1L], to = r[2L], n = n_grid)$y
  }
  .overlap_grid(f1, f2, g[2L] - g[1L])
}
