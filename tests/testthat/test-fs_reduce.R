# fs_reduce() and fs_rotate() -------------------------------------------------

toy_traits <- function(n = 30L, p = 5L, seed = 1L) {
  set.seed(seed)
  m <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("sp", seq_len(n)),
                              paste0("tr", seq_len(p))))
  as.data.frame(m)
}

test_that("fs_reduce truncates PCA/PCoA spaces", {
  sp <- fs_space(toy_traits(), method = "pca")
  r <- fs_reduce(sp, 3L)
  expect_identical(ncol(r$coords), 3L)
  expect_identical(ncol(r$loadings), 3L)
  expect_length(r$eig, 3L)
  expect_true(r$reduced)
  expect_identical(r$dims_full, 5L)
  # truncation preserves the retained axes exactly
  expect_equal(r$coords, sp$coords[, 1:3])
})

test_that("fs_reduce refits NMDS at the requested k", {
  sp <- fs_space(fs_dist(toy_traits()), method = "nmds", k = 3L, seed = 1L)
  r <- fs_reduce(sp, 2L)
  expect_identical(ncol(r$coords), 2L)
  expect_true(is.numeric(r$stress))
})

test_that("fs_reduce validates dims", {
  sp <- fs_space(toy_traits(), method = "pca")
  expect_error(fs_reduce(sp, 99L), "exceeds")
  expect_error(fs_reduce(sp, 0L), "positive")
  expect_error(fs_reduce(sp, c(1L, 2L)), "single")
})

test_that("fs_rotate requires a reduced PCA space", {
  sp <- fs_space(toy_traits(), method = "pca")
  expect_error(fs_rotate(sp), "fs_reduce")
  pco <- fs_space(fs_dist(toy_traits()), method = "pcoa")
  expect_error(fs_rotate(fs_reduce(pco, 3L)), "PCA")
})

test_that("fs_rotate rigid (default) preserves geometry exactly", {
  sp <- fs_reduce(fs_space(toy_traits(), method = "pca"), 3L)
  rot <- fs_rotate(sp)
  expect_identical(rot$rotation_type, "rigid")
  # orthogonal rotation matrix
  expect_equal(crossprod(rot$rotmat), diag(3L), ignore_attr = TRUE)
  # all pairwise distances are untouched
  expect_equal(dist(rot$coords), dist(sp$coords), ignore_attr = TRUE)
  # total variance preserved (trace invariance)
  expect_equal(sum(rot$eig_rotated), sum(unname(sp$eig)),
               tolerance = 1e-10)
  # loadings are exactly the trait-score correlations (scale = TRUE)
  expect_equal(unname(stats::cor(as.matrix(sp$traits), rot$coords)),
               unname(rot$loadings), tolerance = 1e-10)
  # projecting the original traits reproduces the rotated coordinates
  pr <- fs_project(rot, sp$traits)
  expect_equal(pr, rot$coords, tolerance = 1e-8)
  # a second rotation is refused
  expect_error(fs_rotate(rot), "already rotated")
})

test_that("fs_rotate rescaled follows the factor-analytic convention", {
  sp <- fs_reduce(fs_space(toy_traits(), method = "pca"), 3L)
  rot <- fs_rotate(sp, type = "rescaled")
  expect_identical(rot$rotation, "varimax")
  # the rotation matrix is orthogonal
  expect_equal(crossprod(rot$rotmat), diag(3L), ignore_attr = TRUE)
  # total variance is redistributed, not created or lost
  expect_equal(sum(rot$eig_rotated), sum(unname(sp$eig)),
               tolerance = 1e-10)
  # the loading structure (communalities) is invariant under rotation
  expect_equal(rot$loadings %*% t(rot$loadings),
               sp$loadings %*% t(sp$loadings),
               tolerance = 1e-10, ignore_attr = TRUE)
  # rotated scores have variance equal to the rotated SS loadings
  cvar <- apply(rot$coords, 2L, function(z) mean((z - mean(z))^2))
  expect_equal(unname(cvar), unname(rot$eig_rotated), tolerance = 1e-8)
  # projecting the original traits reproduces the rotated coordinates
  pr <- fs_project(rot, sp$traits)
  expect_equal(pr, rot$coords, tolerance = 1e-8)
})

test_that("reduction drops stale TPDs with a warning", {
  sp <- fs_space(toy_traits(), method = "pca")
  sp$tpds <- list(alpha = 0.99)  # simulate estimated TPDs
  expect_warning(r <- fs_reduce(sp, 2L), "invalidates")
  expect_null(r$tpds)
})
