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
  sp <- fs_space(toy_traits(), method = "nmds", seed = 1L)
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
  pco <- fs_space(toy_traits(), method = "pcoa")
  expect_error(fs_rotate(fs_reduce(pco, 3L)), "PCA")
})

test_that("fs_rotate returns an orthogonal rotation of the scores", {
  sp <- fs_reduce(fs_space(toy_traits(), method = "pca"), 3L)
  rot <- fs_rotate(sp)
  expect_identical(rot$rotation, "varimax")
  # rotation preserves total variance of the retained axes
  expect_equal(sum(apply(rot$coords, 2L, var)),
               sum(apply(sp$coords, 2L, var)))
  # and preserves pairwise distances between units
  expect_equal(dist(rot$coords), dist(sp$coords), ignore_attr = TRUE)
})

test_that("reduction drops stale TPDs with a warning", {
  sp <- fs_space(toy_traits(), method = "pca")
  sp$tpds <- list(alpha = 0.99)  # simulate estimated TPDs
  expect_warning(r <- fs_reduce(sp, 2L), "invalidates")
  expect_null(r$tpds)
})
