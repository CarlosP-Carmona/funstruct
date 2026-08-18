# fs_itv(): effect of intraspecific variability on the trait space ------------

#' Swapping procedure: intraspecific variability and trait-space properties
#'
#' Implements the swapping procedure of Puglielli et al. (2025) to quantify
#' how intraspecific trait variability (ITV) affects the main properties of
#' a between-species trait space. For every individual observation, the
#' corresponding species' mean is replaced by that observation, the PCA is
#' recomputed, and two target parameters are extracted per principal
#' component: the angle between the original and the swapped eigenvector
#' (degrees; 0 = axis unchanged) and the eigenvalue ratio
#' (swapped/original; 1 = variance structure unchanged).
#'
#' The procedure is exhaustive and deterministic: every observation is
#' swapped once. It applies to PCA spaces built on species mean traits.
#'
#' @param space A PCA `fspace` built on species mean traits (species as row
#'   names).
#' @param obs Data.frame of individual-level observations with the same
#'   trait columns as the space.
#' @param ids Vector (length `nrow(obs)`) assigning each observation to a
#'   species; all values must match species names in the space.
#'
#' @return A list of class `fs_itv` with elements `results` (data.frame:
#'   `species`, `obs`, `pc`, `angle`, `eig_ratio`) and `summary`
#'   (per-component means and quantiles across swaps).
#' @references Puglielli, G., Carmona, C.P., Bissi, A. & Tordoni, E. (2025)
#'   Quantifying the influence of intraspecific variability in trait
#'   spaces. *npj Biodiversity*, 4, 36.
#' @seealso [fs_space()]
#' @export
fs_itv <- function(space, obs, ids) {
  stopifnot(is_fspace(space))
  if (space$method != "pca" || is.null(space$traits)) {
    stop("The swapping procedure requires a PCA space with stored traits ",
         "(species means).", call. = FALSE)
  }
  M <- as.matrix(space$traits)
  species <- rownames(M)
  ids <- as.character(ids)
  if (length(ids) != nrow(obs)) {
    stop("`ids` must have one entry per row of `obs`.", call. = FALSE)
  }
  unknown <- setdiff(unique(ids), species)
  if (length(unknown)) {
    stop("Observations refer to species absent from the space: ",
         paste(utils::head(unknown, 5L), collapse = ", "),
         if (length(unknown) > 5L) ", ...", call. = FALSE)
  }
  missing_tr <- setdiff(colnames(M), colnames(obs))
  if (length(missing_tr)) {
    stop("`obs` lacks trait column(s): ",
         paste(missing_tr, collapse = ", "), call. = FALSE)
  }
  O <- as.matrix(as.data.frame(obs)[, colnames(M), drop = FALSE])
  if (anyNA(O)) {
    stop("`obs` contains missing trait values; complete or impute first.",
         call. = FALSE)
  }
  scl <- !identical(space$scale, FALSE)

  ref <- stats::prcomp(M, center = TRUE, scale. = scl)
  npc <- ncol(ref$rotation)

  res <- vector("list", nrow(O))
  for (i in seq_len(nrow(O))) {
    M2 <- M
    M2[ids[i], ] <- O[i, ]
    sw <- stats::prcomp(M2, center = TRUE, scale. = scl)
    kk <- min(npc, ncol(sw$rotation))
    dots <- abs(colSums(ref$rotation[, seq_len(kk), drop = FALSE] *
                        sw$rotation[, seq_len(kk), drop = FALSE]))
    angle <- acos(pmin(1, dots)) * 180 / pi
    ratio <- (sw$sdev[seq_len(kk)]^2) / (ref$sdev[seq_len(kk)]^2)
    res[[i]] <- data.frame(species = ids[i], obs = i, pc = seq_len(kk),
                           angle = angle, eig_ratio = ratio)
  }
  results <- do.call(rbind, res)
  rownames(results) <- NULL

  smr <- do.call(rbind, lapply(split(results, results$pc), function(z) {
    data.frame(pc = z$pc[1L],
               n_swaps = nrow(z),
               angle_mean = mean(z$angle),
               angle_q95 = stats::quantile(z$angle, 0.95, names = FALSE),
               eig_ratio_mean = mean(z$eig_ratio),
               eig_ratio_q05 = stats::quantile(z$eig_ratio, 0.05,
                                               names = FALSE),
               eig_ratio_q95 = stats::quantile(z$eig_ratio, 0.95,
                                               names = FALSE))
  }))
  rownames(smr) <- NULL

  out <- list(results = results, summary = smr,
              n_species = length(unique(ids)), n_obs = nrow(O))
  class(out) <- "fs_itv"
  out
}

#' @export
print.fs_itv <- function(x, ...) {
  cat("<fs_itv> swapping procedure: ", x$n_obs, " observations across ",
      x$n_species, " species\n", sep = "")
  cat("Effect of ITV per principal component:\n")
  print.data.frame(cbind(x$summary["pc"],
                         round(x$summary[, -1L], 3)),
                   row.names = FALSE)
  cat("(angle in degrees: 0 = axis unchanged; eig_ratio: 1 = variance ",
      "unchanged)\n", sep = "")
  invisible(x)
}
