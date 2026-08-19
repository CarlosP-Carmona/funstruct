# fs_trait_dim() and fs_dimensionality(method = "end") ------------------------

# two independent factors, each measured by two (near-)duplicate traits
blocky_traits <- function(n = 400L, seed = 2L) {
  set.seed(seed)
  f1 <- rnorm(n); f2 <- rnorm(n)
  data.frame(a1 = f1, a2 = f1 + rnorm(n, sd = 0.05),
             b1 = f2, b2 = f2 + rnorm(n, sd = 0.05),
             row.names = paste0("sp", seq_len(n)))
}

test_that("END recovers the effective structure of the eigenvalues", {
  x <- blocky_traits()
  sp <- fs_space(x, method = "pca")
  d <- fs_dimensionality(sp, method = "end")
  # two near-duplicated factors -> about two effective dimensions
  expect_gt(d$end, 1.8)
  expect_lt(d$end, 2.3)
  expect_identical(d$suggested, max(1L, as.integer(round(d$end))))
  expect_identical(colnames(d$curve), c("axis", "eigenvalue", "share"))
  # equal shares -> END equals the number of axes (exactly, by formula)
  expect_equal(.end_hill(rep(1, 5L), q = 2), 5)
  expect_equal(.end_hill(rep(0.2, 4L), q = 1), 4)
  # one dominant axis -> END approaches 1
  expect_lt(.end_hill(c(100, 0.01, 0.01), q = 2), 1.01)
  # non-PCA spaces are refused
  spo <- fs_space(fs_dist(x), method = "pcoa")
  expect_error(fs_dimensionality(spo, method = "end"), "pca")
})

test_that("fs_trait_dim separates redundant from independent traits", {
  x <- blocky_traits()
  # add one trait independent of both factors
  set.seed(9)
  x$c1 <- rnorm(nrow(x))
  td <- fs_trait_dim(x, n_null = 49L, seed = 1L)
  expect_s3_class(td, "fs_trait_dim")
  tb <- td$table
  expect_setequal(tb$trait, colnames(x))
  # the duplicate traits add almost nothing; the independent one adds
  # about as much as a random trait
  expect_lt(tb$ratio[tb$trait == "a2"], 0.35)
  expect_gt(tb$ratio[tb$trait == "c1"], 0.7)
  expect_lt(tb$p_random[tb$trait == "a2"], 0.05)
  # null contributions are positive and near one effective dimension
  expect_true(all(tb$null_mean > 0.3))
  # print and plot run
  expect_output(print(td), "END of the full space")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(td))
})

test_that("fs_trait_dim single-trait mode and fspace input", {
  x <- blocky_traits(n = 200L)
  td1 <- fs_trait_dim(x, trait = "b2", n_null = 29L, seed = 3L)
  expect_identical(nrow(td1$table), 1L)
  expect_identical(td1$table$trait, "b2")
  sp <- fs_space(x, method = "pca")
  td2 <- fs_trait_dim(sp, trait = "b2", n_null = 29L, seed = 3L)
  expect_equal(td1$table$delta, td2$table$delta, tolerance = 1e-12)
})

test_that("fs_trait_dim validates inputs", {
  x <- blocky_traits(n = 50L)
  expect_error(fs_trait_dim(x[, 1:2]), "three traits")
  expect_error(fs_trait_dim(x, trait = "nope"), "not found")
  xna <- x; xna[1L, 1L] <- NA
  expect_error(fs_trait_dim(xna), "impute")
  xm <- x; xm$cat <- "a"
  expect_error(fs_trait_dim(xm), "numeric")
  spo <- fs_space(fs_dist(x), method = "pcoa")
  expect_error(fs_trait_dim(spo), "PCA space")
})
