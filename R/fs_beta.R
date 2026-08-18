# fs_beta(): functional beta diversity between assemblages --------------------

#' Functional beta diversity between assemblages
#'
#' Pairwise functional dissimilarity between assemblages.
#'
#' * `engine = "prob"`: overlap-based dissimilarity between assemblage
#'   TPDs (1 minus the shared probability mass), with the decomposition of
#'   the TPD framework (ported from `TPD::dissim`; Carmona et al. 2019):
#'   the dissimilarity is split into the proportion attributable to
#'   probability differences inside the co-occupied region (`P_shared`)
#'   and the proportion attributable to regions occupied by only one
#'   assemblage (`P_non_shared`).
#' * `engine = "points"`: Rao-based dissimilarity (Rao's DISC): the
#'   quadratic entropy of the pooled pair minus the mean of the two
#'   within-assemblage entropies, on max-scaled Euclidean distances. No
#'   decomposition is defined for this engine.
#'
#' @param space An `fspace`; the probabilistic engine requires TPDs
#'   ([fs_tpd()]).
#' @param comm Community matrix (assemblages x units) or long data.frame;
#'   see [fs_structure()].
#' @param engine `"prob"` or `"points"`.
#' @param decompose Logical; compute the shared/non-shared decomposition
#'   (probabilistic engine only).
#'
#' @return A list of class `fbeta` with `dissimilarity` (assemblages x
#'   assemblages) and, for the probabilistic engine with
#'   `decompose = TRUE`, `P_shared` and `P_non_shared` (NA on pairs with
#'   zero dissimilarity).
#' @references Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J.
#'   (2019) Trait probability density (TPD): measuring functional
#'   diversity across scales based on TPD with R. *Ecology*, 100, e02876.
#' @seealso [fs_partition()], [fs_structure()]
#' @examples
#' sp <- fs_space(data.frame(size = c(0, 2, 4), shape = c(0, 1, 0),
#'                           row.names = c("sp1", "sp2", "sp3")),
#'                method = "raw", scale = FALSE)
#' tp <- fs_tpd(sp, sds = 0.4)
#' comm <- rbind(A = c(sp1 = 1, sp2 = 1, sp3 = 0),
#'               B = c(sp1 = 0, sp2 = 1, sp3 = 1),
#'               C = c(sp1 = 1, sp2 = 1, sp3 = 1))
#'
#' b <- fs_beta(tp, comm)
#' b$dissimilarity
#' b$P_non_shared   # part of the dissimilarity from non-shared space
#' @export
fs_beta <- function(space, comm, engine = c("prob", "points"),
                    decompose = TRUE) {
  stopifnot(is_fspace(space))
  engine <- match.arg(engine)
  if (engine == "prob" && is.null(space$tpds)) {
    stop("The probabilistic engine needs TPDs: run fs_tpd() first, or ",
         "use engine = \"points\".", call. = FALSE)
  }
  unit_names <- if (engine == "prob") names(space$tpds$units) else
    rownames(space$coords)
  W <- .as_comm(comm, unit_names, relative = TRUE)
  n <- nrow(W)
  if (n < 2L) {
    stop("Beta diversity needs at least two assemblages.", call. = FALSE)
  }
  nm <- rownames(W)
  Dm <- matrix(0, n, n, dimnames = list(nm, nm))
  Psh <- Pns <- if (engine == "prob" && decompose) {
    matrix(NA_real_, n, n, dimnames = list(nm, nm))
  } else NULL

  if (engine == "prob") {
    tpdc <- lapply(seq_len(n), function(a) {
      .tpdc_sparse(space$tpds$units, W[a, ])
    })
    for (i in seq_len(n - 1L)) {
      ci <- tpdc[[i]]$cells; pi_ <- tpdc[[i]]$probs
      for (j in seq(i + 1L, n)) {
        cj <- tpdc[[j]]$cells; pj_ <- tpdc[[j]]$probs
        common <- intersect(ci, cj)
        ii <- match(common, ci); jj <- match(common, cj)
        O <- if (length(common)) sum(pmin(pi_[ii], pj_[jj])) else 0
        dis <- 1 - O
        Dm[i, j] <- Dm[j, i] <- dis
        if (!is.null(Psh)) {
          if (dis > 1e-12) {
            A <- if (length(common)) sum(pmax(pi_[ii], pj_[jj])) - O else 0
            B <- 1 - if (length(common)) sum(pi_[ii]) else 0
            C <- 1 - if (length(common)) sum(pj_[jj]) else 0
            pns <- (2 * min(B, C)) / (A + 2 * min(B, C))
            Pns[i, j] <- Pns[j, i] <- pns
            Psh[i, j] <- Psh[j, i] <- 1 - pns
          }
        }
      }
    }
  } else {
    D <- as.matrix(stats::dist(space$coords))
    D <- D / max(D)
    rao_of <- function(w) as.numeric(t(w) %*% D %*% w)
    rao_a <- apply(W, 1L, rao_of)
    for (i in seq_len(n - 1L)) {
      for (j in seq(i + 1L, n)) {
        pooled <- (W[i, ] + W[j, ]) / 2
        Dm[i, j] <- Dm[j, i] <- rao_of(pooled) - (rao_a[i] + rao_a[j]) / 2
      }
    }
    if (decompose) {
      message("No decomposition is defined for the points engine; ",
              "returning total dissimilarities only.")
    }
  }

  out <- list(engine = engine, dissimilarity = Dm,
              P_shared = Psh, P_non_shared = Pns, call = match.call())
  class(out) <- "fbeta"
  out
}

#' @export
print.fbeta <- function(x, ...) {
  n <- nrow(x$dissimilarity)
  off <- x$dissimilarity[upper.tri(x$dissimilarity)]
  cat("<fbeta> engine: ", x$engine, " | assemblages: ", n, "\n", sep = "")
  cat("Pairwise dissimilarity: mean = ", round(mean(off), 4),
      ", range = [", round(min(off), 4), ", ", round(max(off), 4),
      "]\n", sep = "")
  k <- min(n, 6L)
  print(round(x$dissimilarity[seq_len(k), seq_len(k)], 3))
  if (n > k) cat("... (", n, " x ", n, " matrix)\n", sep = "")
  invisible(x)
}

# assemblage TPD as a sparse abundance-weighted mixture ------------------------

.tpdc_sparse <- function(units, w) {
  present <- which(w > 0)
  parts <- lapply(present, function(u) {
    list(cells = units[[u]]$cells, p = units[[u]]$probs * w[u])
  })
  all_cells <- sort(unique(unlist(lapply(parts, `[[`, "cells"))))
  pc <- numeric(length(all_cells))
  pos <- seq_along(all_cells)
  names(pos) <- all_cells
  for (m in parts) {
    j <- pos[as.character(m$cells)]
    pc[j] <- pc[j] + m$p
  }
  list(cells = all_cells, probs = pc / sum(pc))
}
