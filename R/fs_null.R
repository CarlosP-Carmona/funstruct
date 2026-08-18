# fs_null(): null models and standardized effect sizes ------------------------

#' Null models for functional structure
#'
#' Computes functional structure ([fs_structure()]) together with null
#' expectations, standardized effect sizes (SES) and two-sided p-values.
#' Unit TPDs and bandwidths are estimated once and never re-selected
#' inside the randomizations: null models permute identities or
#' abundances, not the smoothing.
#'
#' Null models:
#' * `"shuffle.units"`: unit identities are permuted globally (the
#'   classical taxa-labels null): each assemblage keeps its abundance
#'   structure, but which trait profile carries each abundance is
#'   randomized consistently across assemblages.
#' * `"shuffle.abund"`: abundances are permuted independently within each
#'   assemblage across all units.
#'
#' @param space An `fspace` (with TPDs for the probabilistic engine).
#' @param comm Community matrix or long data.frame; see [fs_structure()].
#' @param engine,indices Passed to [fs_structure()].
#' @param model Null model; see Details.
#' @param n Number of randomizations (default 499).
#' @param seed Random seed.
#' @param cores Number of cores (unix only; Windows runs serially).
#'
#' @return The observed `fstructure` with additional columns `ses_*` and
#'   `p_*` for every index, and an attribute `nulls` with the null means
#'   and standard deviations.
#' @seealso [fs_structure()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 2, 4), shape = c(0, 1, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)
#' comm <- rbind(A = c(sp1 = 0.7, sp2 = 0.3, sp3 = 0),
#'               B = c(sp1 = 0.2, sp2 = 0.3, sp3 = 0.5))
#'
#' st <- fs_null(tp, comm, indices = c("richness", "dispersion"),
#'               n = 49, seed = 1)
#' st[, c("richness", "ses_richness", "p_richness")]
#' @export
fs_null <- function(space, comm = NULL, engine = c("prob", "points"),
                    indices = NULL,
                    model = c("shuffle.units", "shuffle.abund"),
                    n = 499L, seed = NULL, cores = 1L) {
  engine <- match.arg(engine)
  model <- match.arg(model)
  obs <- fs_structure(space, comm, engine = engine, indices = indices)
  unit_names <- if (engine == "prob") names(space$tpds$units) else
    rownames(space$coords)
  W <- .as_comm(comm, unit_names, relative = TRUE)

  if (.Platform$OS.type == "windows") cores <- 1L
  if (!is.null(seed)) {
    if (cores > 1L) RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }
  one_null <- function(k) {
    Wn <- W
    if (model == "shuffle.units") {
      colnames(Wn) <- sample(colnames(W))
      Wn <- Wn[, unit_names, drop = FALSE]
    } else {
      for (a in seq_len(nrow(Wn))) Wn[a, ] <- sample(Wn[a, ])
    }
    as.matrix(fs_structure(space, Wn, engine = engine, indices = indices))
  }
  nulls <- if (cores > 1L) {
    parallel::mclapply(seq_len(n), one_null, mc.cores = cores)
  } else {
    lapply(seq_len(n), one_null)
  }
  arr <- array(unlist(nulls),
               dim = c(nrow(obs), ncol(obs), n),
               dimnames = list(rownames(obs), colnames(obs), NULL))
  mu <- apply(arr, c(1L, 2L), mean, na.rm = TRUE)
  sdv <- apply(arr, c(1L, 2L), stats::sd, na.rm = TRUE)
  obs_m <- as.matrix(obs)
  ses <- (obs_m - mu) / sdv
  p_hi <- (apply(sweep(arr, c(1L, 2L), obs_m, ">="), c(1L, 2L), sum,
                 na.rm = TRUE) + 1) / (n + 1)
  p_lo <- (apply(sweep(arr, c(1L, 2L), obs_m, "<="), c(1L, 2L), sum,
                 na.rm = TRUE) + 1) / (n + 1)
  # first argument must be the matrix: pmin() keeps the attributes
  # (including dim) of its first argument only
  pval <- pmin(2 * pmin(p_hi, p_lo), 1)

  out <- obs
  for (j in colnames(obs)) {
    out[[paste0("ses_", j)]] <- ses[, j]
    out[[paste0("p_", j)]] <- pval[, j]
  }
  attr(out, "nulls") <- list(model = model, n = n,
                             null_mean = mu, null_sd = sdv)
  attr(out, "engine") <- attr(obs, "engine")
  attr(out, "settings") <- attr(obs, "settings")
  class(out) <- class(obs)
  out
}
