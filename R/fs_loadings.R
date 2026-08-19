# fs_loadings(): post-hoc trait-axis loadings for PCoA/NMDS spaces ------------

#' Post-hoc trait loadings for dissimilarity-based spaces
#'
#' Computes correlations between the original traits and the axes of a
#' PCoA or NMDS space, and stores them as the space's `loadings` so that
#' [plot.fspace()] can draw trait arrows. This formalizes standard
#' practice for dissimilarity-based ordinations (compare
#' `mFD::traits.faxes.cor()` and `vegan::envfit()`).
#'
#' The convention is continuous with PCA: when the dissimilarity is the
#' Euclidean distance on scaled traits, PCoA reproduces the PCA
#' configuration and these correlations coincide (up to axis sign) with
#' the loadings [fs_space()] reports for PCA.
#'
#' Two caveats, both consequences of the axes being functions of
#' dissimilarities rather than linear combinations of traits:
#'
#' * These loadings are *descriptive*, computed after the fact. They do
#'   not define the axes and cannot be used to project new units into the
#'   space ([fs_project()] remains PCA-only).
#' * A correlation only captures the linear part of a trait-axis
#'   relationship. Gower-based axes often relate to traits monotonically
#'   but not linearly; `method = "spearman"` is the more honest choice in
#'   that case.
#'
#' Numeric traits enter the loadings matrix as (signed) correlations.
#' Categorical traits have no direction, so they are not drawn as
#' arrows; their association with each axis is reported separately as the
#' square root of the correlation ratio (eta, in `[0, 1]`, the
#' multiple-correlation analogue used by mFD) in
#' `attr(space$loadings, "categorical")`.
#'
#' @param space A `pcoa` or `nmds` `fspace` (PCA and raw spaces already
#'   carry construction-time loadings, or are the traits themselves).
#' @param traits Matrix or data.frame of the trait values the
#'   dissimilarities were computed from, with row names covering all
#'   units of the space. Missing values are tolerated (pairwise-complete
#'   correlations).
#' @param method Correlation type for the numeric traits: `"pearson"`
#'   (default) or `"spearman"`.
#'
#' @return The `fspace` with `$loadings` filled: a numeric-traits x axes
#'   matrix of correlations, with attributes `method` and (when
#'   categorical traits are present) `categorical` (traits x axes matrix
#'   of eta values).
#' @seealso [fs_space()], [fs_dist()], [plot.fspace()]
#' @examples
#' data(gspff)
#' x <- gspff[1:150, ]
#' sp <- fs_space(fs_dist(x), method = "pcoa")
#' sp <- fs_loadings(sp, x)
#' round(sp$loadings[, 1:2], 2)
#' plot(sp)   # trait arrows, now also for PCoA
#' @export
fs_loadings <- function(space, traits, method = c("pearson", "spearman")) {
  stopifnot(is_fspace(space))
  method <- match.arg(method)
  if (!space$method %in% c("pcoa", "nmds")) {
    stop("Post-hoc loadings are for dissimilarity-based spaces ",
         "(pcoa/nmds). PCA spaces already carry loadings from ",
         "construction, and raw axes are the traits themselves.",
         call. = FALSE)
  }
  traits <- as.data.frame(traits)
  if (is.null(rownames(traits))) {
    stop("`traits` must have row names identifying the units.",
         call. = FALSE)
  }
  units <- rownames(space$coords)
  missing_u <- setdiff(units, rownames(traits))
  if (length(missing_u)) {
    stop("`traits` lacks rows for unit(s): ",
         paste(utils::head(missing_u, 5L), collapse = ", "),
         if (length(missing_u) > 5L) ", ...", call. = FALSE)
  }
  traits <- traits[units, , drop = FALSE]
  co <- space$coords
  num_ok <- vapply(traits, is.numeric, logical(1L))
  if (!any(num_ok) && !any(!num_ok)) {
    stop("`traits` has no usable columns.", call. = FALSE)
  }

  L <- NULL
  if (any(num_ok)) {
    L <- stats::cor(as.matrix(traits[, num_ok, drop = FALSE]), co,
                    use = "pairwise.complete.obs", method = method)
    colnames(L) <- colnames(co)
  } else {
    warning("No numeric traits: the loadings matrix is empty and no ",
            "arrows will be drawn; see the `categorical` attribute.",
            call. = FALSE)
    L <- matrix(numeric(0), 0L, ncol(co),
                dimnames = list(NULL, colnames(co)))
  }
  attr(L, "method") <- method
  attr(L, "posthoc") <- TRUE

  if (any(!num_ok)) {
    cat_names <- names(traits)[!num_ok]
    E <- matrix(NA_real_, length(cat_names), ncol(co),
                dimnames = list(cat_names, colnames(co)))
    for (tn in cat_names) {
      g <- as.character(traits[[tn]])
      for (j in seq_len(ncol(co))) {
        E[tn, j] <- .eta_assoc(co[, j], g)
      }
    }
    attr(L, "categorical") <- E
  }

  space$loadings <- L
  space
}

# square root of the correlation ratio (eta): between-group share of the
# axis variance explained by the levels of a categorical trait
.eta_assoc <- function(z, g) {
  ok <- !is.na(z) & !is.na(g)
  z <- z[ok]
  g <- g[ok]
  if (length(z) < 3L || length(unique(g)) < 2L) return(NA_real_)
  mu <- mean(z)
  m <- tapply(z, g, mean)
  n <- tapply(z, g, length)
  ss_b <- sum(n * (m - mu)^2)
  ss_t <- sum((z - mu)^2)
  if (ss_t <= 0) return(NA_real_)
  sqrt(ss_b / ss_t)
}
