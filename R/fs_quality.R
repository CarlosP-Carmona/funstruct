# fs_quality(): distance-preservation quality of a trait space ----------------

#' Quality of a trait space across dimensionalities
#'
#' Measures how faithfully the first `k` axes of a space preserve the
#' original functional dissimilarities between units, following the logic of
#' Maire et al. (2015). For each candidate dimensionality, the Euclidean
#' distances in the reduced space are compared with the initial
#' dissimilarities (both scaled to unit maximum); the mean absolute
#' deviation (mAD) and mean squared deviation (mSD) summarize the mismatch
#' (0 = perfect preservation).
#'
#' @param space An `fspace` object. The initial dissimilarities are taken
#'   from the distance matrix stored in the space (PCoA/NMDS) or computed
#'   from the stored traits (PCA/raw).
#' @param dims Candidate dimensionalities to evaluate. Default: `1` up to
#'   `min(10, available axes)`.
#'
#' @return A data.frame of class `fs_quality` with columns `dims`, `mAD`
#'   and `mSD`, and an attribute `best` (dimensionality with lowest mSD).
#' @references Maire, E., Grenouillet, G., Brosse, S. & Villeger, S. (2015)
#'   How many dimensions are needed to accurately assess functional
#'   diversity? A pragmatic approach for assessing the quality of
#'   functional spaces. *Global Ecology and Biogeography*, 24, 728-740.
#' @seealso [fs_dimensionality()], [fs_reduce()]
#' @examples
#' data(gspff)
#' sp <- fs_space(gspff[1:300, ], method = "pca")
#' fs_quality(sp)
#' @export
fs_quality <- function(space, dims = NULL) {
  stopifnot(is_fspace(space))
  d_init <- .initial_dist(space)
  kmax <- ncol(space$coords)
  if (is.null(dims)) dims <- seq_len(min(kmax, 10L))
  dims <- sort(unique(as.integer(dims)))
  if (any(is.na(dims)) || any(dims < 1L) || any(dims > kmax)) {
    stop("`dims` must be integers between 1 and ", kmax, ".", call. = FALSE)
  }
  di <- as.vector(d_init)
  di <- di / max(di)
  res <- vapply(dims, function(k) {
    dk <- as.vector(stats::dist(space$coords[, seq_len(k), drop = FALSE]))
    dk <- dk / max(dk)
    dev <- dk - di
    c(mean(abs(dev)), mean(dev^2))
  }, numeric(2L))
  out <- data.frame(dims = dims, mAD = res[1L, ], mSD = res[2L, ])
  attr(out, "best") <- dims[which.min(out$mSD)]
  class(out) <- c("fs_quality", "data.frame")
  out
}

#' @export
print.fs_quality <- function(x, ...) {
  cat("<fs_quality> distance-preservation quality (0 = perfect)\n")
  print.data.frame(cbind(x["dims"], round(x[c("mAD", "mSD")], 4)),
                   row.names = FALSE)
  cat("Lowest mSD at", attr(x, "best"), "dimension(s).\n")
  invisible(x)
}

# initial dissimilarities used as the reference -------------------------------

.initial_dist <- function(space) {
  if (!is.null(space$dist)) return(space$dist)
  if (is.null(space$traits)) {
    stop("This space stores neither a distance matrix nor the trait table, ",
         "so the initial dissimilarities cannot be recovered. Rebuild with ",
         "fs_space() or supply traits to as_fspace().", call. = FALSE)
  }
  numeric_ok <- vapply(as.data.frame(space$traits), is.numeric, logical(1L))
  if (!all(numeric_ok) || anyNA(space$traits)) {
    return(cluster::daisy(as.data.frame(space$traits), metric = "gower"))
  }
  scl <- !identical(space$scale, FALSE)
  stats::dist(base::scale(as.matrix(space$traits), center = TRUE,
                          scale = scl))
}
