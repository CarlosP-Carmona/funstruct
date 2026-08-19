# fs_means(): unit means from individual observations -------------------------

#' Aggregate individual observations into unit means
#'
#' Computes one mean trait value per unit (e.g. per species) from
#' individual-level observations. This is the entry point of the
#' mean-based workflow for unbalanced sampling: when units differ in
#' replication, an ordination built directly on the individuals lets the
#' well-sampled units pull the axes. Building the space on one mean per
#' unit weights every unit equally, and the individuals are then
#' projected into that space afterwards -- with [fs_project()] for the
#' coordinates, or in one step with [fs_tpd()]`(space, obs = , ids = )`
#' for trait probability densities.
#'
#' @param obs Data.frame or matrix of individual observations (rows) by
#'   numeric traits (columns).
#' @param ids Vector (length `nrow(obs)`) assigning each observation to a
#'   unit.
#'
#' @return A data.frame of per-unit trait means (units as row names, in
#'   order of first appearance in `ids`), with attributes `n` (named
#'   integer vector of observations per unit) and `sds` (data.frame of
#'   per-unit, per-trait standard deviations; `NA` for units with a
#'   single observation). Missing values are ignored within each unit.
#' @seealso [fs_space()], [fs_project()], [fs_tpd()]
#' @examples
#' data(grassland)
#' tr <- data.frame(height = log10(grassland$height),
#'                  sla = log10(grassland$sla),
#'                  row.names = paste0("ind", seq_len(nrow(grassland))))
#'
#' trM <- fs_means(tr, grassland$species)   # one row per species
#' head(trM)
#' head(attr(trM, "n"))
#'
#' # mean-based space, individuals projected afterwards:
#' sp <- fs_space(trM, method = "pca")
#' co <- fs_project(sp, tr)                 # individual coordinates
#' head(co)
#' @export
fs_means <- function(obs, ids) {
  obs <- as.data.frame(obs)
  if (ncol(obs) < 1L || nrow(obs) < 1L) {
    stop("`obs` must have at least one row and one column.", call. = FALSE)
  }
  num_ok <- vapply(obs, is.numeric, logical(1L))
  if (!all(num_ok)) {
    stop("All columns of `obs` must be numeric; found: ",
         paste(names(obs)[!num_ok], collapse = ", "), call. = FALSE)
  }
  ids <- as.character(ids)
  if (length(ids) != nrow(obs)) {
    stop("`ids` must have one entry per row of `obs`.", call. = FALSE)
  }
  f <- factor(ids, levels = unique(ids))
  M <- as.data.frame(lapply(obs, function(v) {
    as.vector(tapply(v, f, mean, na.rm = TRUE))
  }), row.names = levels(f))
  S <- as.data.frame(lapply(obs, function(v) {
    as.vector(tapply(v, f, stats::sd, na.rm = TRUE))
  }), row.names = levels(f))
  attr(M, "n") <- stats::setNames(as.integer(table(f)), levels(f))
  attr(M, "sds") <- S
  M
}
