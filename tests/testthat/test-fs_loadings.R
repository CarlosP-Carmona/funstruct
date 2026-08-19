# fs_loadings() ---------------------------------------------------------------

toy_traits <- function(n = 30L, p = 5L, seed = 1L) {
  set.seed(seed)
  m <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("sp", seq_len(n)),
                              paste0("tr", seq_len(p))))
  as.data.frame(m)
}

test_that("Euclidean PCoA loadings match PCA loadings up to axis sign", {
  x <- toy_traits()
  spP <- fs_space(x, method = "pca")
  d <- dist(scale(as.matrix(x)))
  spo <- fs_loadings(fs_space(d, method = "pcoa", correction = "none"), x)
  k <- min(ncol(spP$loadings), ncol(spo$loadings))
  expect_equal(abs(unname(spo$loadings[, 1:k])),
               abs(unname(spP$loadings[, 1:k])), tolerance = 1e-8)
  expect_identical(attr(spo$loadings, "method"), "pearson")
  expect_true(isTRUE(attr(spo$loadings, "posthoc")))
})

test_that("categorical traits get eta associations, not arrows", {
  x <- toy_traits(n = 24L, p = 2L)
  # a factor aligned with tr1: high vs low
  x$grp <- ifelse(x$tr1 > stats::median(x$tr1), "hi", "lo")
  sp <- fs_loadings(fs_space(fs_dist(x), method = "pcoa"), x)
  expect_identical(rownames(sp$loadings), c("tr1", "tr2"))
  E <- attr(sp$loadings, "categorical")
  expect_identical(rownames(E), "grp")
  expect_true(all(E >= 0 & E <= 1, na.rm = TRUE))
  # the factor was built from tr1, so it must associate strongly with the
  # axis that tr1 loads on most
  ax <- which.max(abs(sp$loadings["tr1", ]))
  expect_gt(E["grp", ax], 0.5)
})

test_that("spearman option and NA tolerance work", {
  x <- toy_traits(n = 20L, p = 3L)
  x[1L, 2L] <- NA
  sp <- fs_space(fs_dist(x), method = "pcoa")
  sp <- fs_loadings(sp, x, method = "spearman")
  expect_identical(attr(sp$loadings, "method"), "spearman")
  expect_false(anyNA(sp$loadings))
})

test_that("fs_loadings validates its inputs", {
  x <- toy_traits(n = 15L, p = 3L)
  spP <- fs_space(x, method = "pca")
  expect_error(fs_loadings(spP, x), "PCA")
  spo <- fs_space(fs_dist(x), method = "pcoa")
  expect_error(fs_loadings(spo, x[1:10, ]), "lacks rows")
  bad <- x; rownames(bad) <- NULL
  expect_error(fs_loadings(spo, as.matrix(bad)), "row names")
})

test_that("fs_reduce trims post-hoc loadings and their attributes", {
  x <- toy_traits(n = 20L, p = 4L)
  x$cat <- factor(rep(letters[1:2], 10L))
  sp <- fs_loadings(fs_space(fs_dist(x), method = "pcoa"), x)
  r <- fs_reduce(sp, 2L)
  expect_identical(ncol(r$loadings), 2L)
  expect_true(isTRUE(attr(r$loadings, "posthoc")))
  expect_identical(ncol(attr(r$loadings, "categorical")), 2L)
})

test_that("NMDS refit drops stale post-hoc loadings with a warning", {
  x <- toy_traits(n = 20L, p = 4L)
  sp <- fs_space(fs_dist(x), method = "nmds", k = 3L, seed = 1L)
  sp <- fs_loadings(sp, x)
  expect_identical(ncol(sp$loadings), 3L)
  expect_warning(r <- fs_reduce(sp, 2L), "fs_loadings")
  expect_null(r$loadings)
})
