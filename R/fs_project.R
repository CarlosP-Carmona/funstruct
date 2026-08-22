# fs_project(): project new units into an existing PCA space -----------------

#' Project new units into an existing trait space
#'
#' Places new units (e.g. species with imputed traits) into a PCA space
#' built from other units (e.g. only species with complete, empirical
#' trait measurements). The new trait values are centred and scaled with
#' the parameters of the data the space was built from -- stored in the
#' space at construction time -- so both sets of units end up in exactly
#' the same coordinate system. Works for reduced and varimax-rotated
#' spaces.
#'
#' This supports the recommended imputation workflow: build the space
#' from empirical measurements only, then project imputed species into it
#' rather than letting imputed values shape the axes (see [fs_impute()]).
#' The same logic handles unbalanced individual-level sampling: build the
#' space on one mean per unit ([fs_means()]) so every unit weighs equally
#' in the ordination, then project the individual observations into it.
#'
#' @param space A PCA `fspace` (from [fs_space()], possibly after
#'   [fs_reduce()] / [fs_rotate()]).
#' @param newdata Data.frame or matrix of new units (rows) with the same
#'   trait columns as the space, and unique row names. If `newdata` comes
#'   from [fs_impute()] (it carries the `imputed` attribute), only the
#'   units with at least one imputed value are projected: the fully
#'   empirical units are the material the space should be built from, not
#'   projected into it (a message reports how many were skipped). To
#'   override and project every row, remove the attribute first:
#'   `attr(newdata, "imputed") <- NULL`.
#' @param add Logical; if `TRUE`, return the `fspace` with the projected
#'   units appended to its coordinates (flagged `projected` in the unit
#'   metadata) instead of the bare coordinate matrix. Units whose names
#'   are already present in the space are not appended: the coordinates
#'   the space was built with take precedence (with a warning listing the
#'   skipped units).
#'
#' @return A matrix of coordinates for the new units (default), or the
#'   extended `fspace` when `add = TRUE`.
#' @seealso [fs_impute()], [fs_means()], [fs_space()]
#' @examples
#' data(gspff)
#' sp <- fs_reduce(fs_space(gspff[1:500, ], method = "pca"), 2)
#' co <- fs_project(sp, gspff[501:520, ])
#' head(co)
#'
#' # append them to the space instead:
#' sp2 <- fs_project(sp, gspff[501:520, ], add = TRUE)
#' tail(sp2$units)
#' @export
fs_project <- function(space, newdata, add = FALSE) {
  stopifnot(is_fspace(space))
  if (space$method != "pca" || is.null(space$proj) ||
      is.null(space$center)) {
    stop("Projection requires a PCA space built with fs_space() (or ",
         "imported from prcomp/princomp with its centring information).",
         call. = FALSE)
  }
  traits_ref <- names(space$center)
  nd <- as.data.frame(newdata)
  imp <- attr(newdata, "imputed")
  if (!is.null(imp) && !is.null(imp$imputed)) {
    keep_imp <- intersect(rownames(nd), imp$imputed)
    n_skip <- nrow(nd) - length(keep_imp)
    if (!length(keep_imp)) {
      stop("`newdata` comes from fs_impute() but contains no units with ",
           "imputed values; there is nothing to project. (Fully empirical ",
           "units belong in the space itself, not projected into it.)",
           call. = FALSE)
    }
    if (n_skip > 0L) {
      message("Projecting the ", length(keep_imp), " unit(s) with imputed ",
              "values; ", n_skip, " fully empirical unit(s) skipped (see ",
              "attr(newdata, \"imputed\")$complete).")
    }
    nd <- nd[keep_imp, , drop = FALSE]
  }
  missing_tr <- setdiff(traits_ref, colnames(nd))
  if (length(missing_tr)) {
    stop("`newdata` lacks trait column(s): ",
         paste(missing_tr, collapse = ", "), call. = FALSE)
  }
  O <- as.matrix(nd[, traits_ref, drop = FALSE])
  if (anyNA(O)) {
    stop("`newdata` contains missing values; impute first ",
         "(see fs_impute()).", call. = FALSE)
  }
  if (is.null(rownames(O)) || anyDuplicated(rownames(O)) > 0L) {
    stop("`newdata` must have unique row names.", call. = FALSE)
  }
  scl <- space$scale_values
  if (is.null(scl)) scl <- rep(1, length(traits_ref))
  Z <- sweep(sweep(O, 2L, space$center), 2L, scl, "/")
  co <- Z %*% space$proj
  colnames(co) <- colnames(space$coords)
  if (!add) return(co)

  clash <- intersect(rownames(co), rownames(space$coords))
  if (length(clash)) {
    warning(length(clash), " unit(s) already present in the space were ",
            "not appended (the coordinates the space was built with are ",
            "kept): ", paste(utils::head(clash, 5L), collapse = ", "),
            if (length(clash) > 5L) ", ...", call. = FALSE)
    keep <- setdiff(rownames(co), clash)
    if (!length(keep)) return(space)
    co <- co[keep, , drop = FALSE]
    O <- O[keep, , drop = FALSE]
  }
  if (!is.null(space$tpds)) {
    warning("Appending units invalidates the estimated TPDs and any ",
            "aggregated levels built from them; all have been removed. ",
            "Re-run fs_tpd() (and fs_aggregate()).", call. = FALSE)
    space$tpds <- NULL
    space$bw <- NULL
  }
  space$coords <- rbind(space$coords, co)
  if (is.null(space$units$projected)) space$units$projected <- FALSE
  new_units <- data.frame(id = rownames(co), n_obs = 1L,
                          has_own_obs = FALSE, imputed_traits = FALSE,
                          projected = TRUE, stringsAsFactors = FALSE)
  if (!is.null(imp)) {
    new_units$imputed_traits <- new_units$id %in% imp$imputed
  }
  space$units <- rbind(space$units, new_units)
  if (!is.null(space$traits)) {
    space$traits <- rbind(space$traits,
                          as.data.frame(O)[, colnames(space$traits),
                                           drop = FALSE])
  }
  space
}
