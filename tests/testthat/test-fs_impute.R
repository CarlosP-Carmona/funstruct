# fs_impute() ------------------------------------------------------------------

toy_gappy <- function(n = 40L, p = 4L, prop_na = 0.1, seed = 3L) {
  set.seed(seed)
  f <- rnorm(n)
  m <- sapply(seq_len(p), function(j) f + rnorm(n, sd = 0.3))
  dimnames(m) <- list(paste0("sp", seq_len(n)), paste0("tr", seq_len(p)))
  df <- as.data.frame(m)
  idx <- cbind(sample(n, round(n * p * prop_na), replace = TRUE),
               sample(p, round(n * p * prop_na), replace = TRUE))
  idx <- unique(idx)
  df[idx] <- NA
  df
}

test_that("fs_impute completes traits and records what was imputed", {
  skip_if_not_installed("missForest")
  x <- toy_gappy()
  out <- fs_impute(x, seed = 1L)
  expect_false(anyNA(out))
  imp <- attr(out, "imputed")
  expect_identical(dim(imp$cells), dim(as.matrix(x)))
  expect_identical(sum(imp$cells), sum(is.na(x)))
  expect_setequal(imp$imputed, rownames(x)[rowSums(is.na(x)) > 0L])
  expect_setequal(imp$complete, rownames(x)[rowSums(is.na(x)) == 0L])
  # non-missing values are untouched
  keep <- !is.na(x)
  expect_equal(as.matrix(out)[keep], as.matrix(x)[keep])
})

test_that("imputation flags propagate into fs_space unit metadata", {
  skip_if_not_installed("missForest")
  x <- toy_gappy()
  out <- fs_impute(x, seed = 1L)
  sp <- fs_space(out, method = "pca")
  flagged <- sp$units$id[sp$units$imputed_traits]
  expect_setequal(flagged, attr(out, "imputed")$imputed)
})

test_that("fs_impute handles the no-NA case and validates inputs", {
  skip_if_not_installed("missForest")
  x <- toy_gappy()
  x[is.na(x)] <- 0
  expect_message(out <- fs_impute(x), "No missing values")
  expect_identical(attr(out, "imputed")$imputed, character(0))
  bad <- x
  rownames(bad) <- NULL
  expect_error(fs_impute(bad), "row names")
})

test_that("phylogenetic eigenvectors route runs and drops helper columns", {
  skip_if_not_installed("missForest")
  skip_if_not_installed("ape")
  x <- toy_gappy(n = 20L)
  set.seed(7)
  tree <- ape::rtree(20L, tip.label = rownames(x))
  out <- fs_impute(x, phylo = tree, n_eigen = 5L, seed = 1L)
  expect_false(anyNA(out))
  expect_identical(colnames(out), colnames(x))
  # a tree lacking tips errors clearly
  tree2 <- ape::drop.tip(tree, "sp1")
  expect_error(fs_impute(x, phylo = tree2), "absent from the phylogeny")
})
