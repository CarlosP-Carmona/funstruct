# fs_angles() ------------------------------------------------------------------

toy_traits <- function(n = 30L, p = 5L, seed = 1L) {
  set.seed(seed)
  m <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("sp", seq_len(n)),
                              paste0("tr", seq_len(p))))
  as.data.frame(m)
}

blocky_traits <- function(n = 400L, seed = 2L) {
  set.seed(seed)
  f1 <- rnorm(n); f2 <- rnorm(n)
  data.frame(a1 = f1, a2 = f1 + rnorm(n, sd = 0.05),
             b1 = f2, b2 = f2 + rnorm(n, sd = 0.05),
             row.names = paste0("sp", seq_len(n)))
}

test_that("full-space PCA angles equal acos of the trait correlations", {
  x <- toy_traits()
  sp <- fs_space(x, method = "pca")
  an <- fs_angles(sp)
  ref <- acos(pmin(pmax(stats::cor(as.matrix(x)), -1), 1)) * 180 / pi
  diag(ref) <- 0
  expect_equal(unname(an$angles), unname(ref), tolerance = 1e-8)
  expect_true(isSymmetric(an$angles, tol = 1e-10))
})

test_that("reduced-space angles read the trait structure", {
  x <- blocky_traits()
  sp2 <- fs_reduce(fs_space(x, method = "pca"), 2L)
  an <- fs_angles(sp2)
  A <- an$angles
  expect_lt(A["a1", "a2"], 15)                    # duplicates align
  expect_gt(A["a1", "b1"], 75)                    # blocks orthogonal
  expect_lt(A["a1", "b1"], 105)
  expect_true(all(A >= 0 & A <= 180))
})

test_that("dims, center and radians options work", {
  x <- toy_traits()
  sp <- fs_space(x, method = "pca")
  a2 <- fs_angles(sp, dims = 1:2)
  expect_identical(a2$dims, 1:2)
  expect_error(fs_angles(sp, dims = 1L), "two axes")
  expect_error(fs_angles(sp, dims = 99L), "between 1 and")
  ar <- fs_angles(sp, degrees = FALSE)
  expect_true(all(ar$angles <= pi + 1e-12))
  ac <- fs_angles(sp, center = TRUE)
  expect_false(isTRUE(all.equal(ac$angles, fs_angles(sp)$angles)))
})

test_that("profile measures angle preservation across dimensionalities", {
  x <- blocky_traits()
  sp <- fs_space(x, method = "pca")
  ap <- fs_angles(sp, profile = TRUE)
  pr <- ap$profile
  expect_identical(pr$dims, 2:4)
  # at the full dimensionality the reference is reproduced exactly
  expect_equal(pr$cor[pr$dims == 4L], 1, tolerance = 1e-12)
  expect_equal(pr$mad[pr$dims == 4L], 0, tolerance = 1e-12)
  # the structure is 2-dimensional, so k = 2 already preserves it well
  expect_gt(pr$cor[pr$dims == 2L], 0.95)
  expect_lt(pr$mad[pr$dims == 2L], 10)
  # per-dimensionality matrices stored, and consistent with the profile
  expect_named(ap$by_dims, c("d2", "d3", "d4"))
  expect_equal(ap$by_dims$d4, ap$angles, tolerance = 1e-12)
  # dims and profile cannot be combined
  expect_error(fs_angles(sp, dims = 1:2, profile = TRUE), "combined")
  # profile plot needs a stored profile
  expect_error(plot(fs_angles(sp), which = "profile"), "profile = TRUE")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(ap, which = "profile"))
  expect_output(print(ap), "preservation")
})

test_that("PCoA route needs fs_loadings first, then works", {
  x <- toy_traits(n = 40L)
  spo <- fs_space(fs_dist(x), method = "pcoa")
  expect_error(fs_angles(spo), "fs_loadings")
  spo <- fs_loadings(spo, x)
  an <- fs_angles(fs_reduce(spo, 3L))
  expect_identical(rownames(an$angles), colnames(x))
  expect_true(all(an$angles >= 0 & an$angles <= 180))
  # print and plot run
  expect_output(print(an), "pairwise trait angles")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(an))
})
