# as_fspace() -----------------------------------------------------------------

toy_traits <- function(n = 30L, p = 5L, seed = 1L) {
  set.seed(seed)
  m <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("sp", seq_len(n)),
                              paste0("tr", seq_len(p))))
  as.data.frame(m)
}

test_that("prcomp and princomp imports agree with fs_space", {
  x <- toy_traits()
  ref <- fs_space(x, method = "pca")
  # princomp import is exactly fs_space's backend
  imp2 <- as_fspace(princomp(as.matrix(x), cor = TRUE))
  expect_equal(unname(imp2$coords), unname(ref$coords))
  expect_equal(unname(imp2$loadings), unname(ref$loadings))
  # prcomp uses the n-1 divisor: same configuration up to a constant
  imp1 <- as_fspace(prcomp(as.matrix(x), center = TRUE, scale. = TRUE))
  cmp <- fs_compare(imp1, ref)
  expect_lt(cmp$m2, 1e-10)
  # both importers expose true loadings (squared sums = eigenvalues)
  expect_equal(unname(colSums(imp1$loadings^2)), unname(imp1$eig))
})

test_that("as_fspace.dist performs PCoA", {
  x <- toy_traits()
  d <- dist(scale(as.matrix(x)))
  imp <- as_fspace(d, correction = "cailliez")
  ref <- fs_space(d, method = "pcoa", correction = "cailliez")
  expect_equal(imp$coords, ref$coords)
})

test_that("as_fspace.matrix wraps coordinates", {
  co <- matrix(rnorm(20), 10, 2,
               dimnames = list(paste0("sp", 1:10), c("A1", "A2")))
  imp <- as_fspace(co)
  expect_s3_class(imp, "fspace")
  expect_identical(imp$method, "raw")
  expect_error(as_fspace(unname(co)), "row names")
})

test_that("as_fspace.default rejects unknown classes", {
  expect_error(as_fspace(lm(y ~ x, data.frame(x = 1:5, y = 1:5))),
               "Don't know")
})

test_that("vegan and ade4 objects convert when available", {
  skip_if_not_installed("vegan")
  x <- toy_traits()
  fit <- suppressWarnings(
    vegan::metaMDS(dist(scale(as.matrix(x))), k = 2, trace = 0)
  )
  imp <- as_fspace(fit)
  expect_identical(imp$method, "nmds")
  expect_identical(ncol(imp$coords), 2L)

  skip_if_not_installed("ade4")
  du <- ade4::dudi.pca(x, scannf = FALSE, nf = 3)
  imp2 <- as_fspace(du)
  expect_identical(imp2$method, "pca")
  expect_identical(ncol(imp2$coords), 3L)
})
