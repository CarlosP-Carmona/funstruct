# fs_space() ------------------------------------------------------------------

toy_traits <- function(n = 30L, p = 5L, seed = 1L) {
  set.seed(seed)
  m <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("sp", seq_len(n)),
                              paste0("tr", seq_len(p))))
  as.data.frame(m)
}

test_that("PCA space matches prcomp exactly", {
  x <- toy_traits()
  sp <- fs_space(x, method = "pca")
  ref <- prcomp(as.matrix(x), center = TRUE, scale. = TRUE)
  expect_s3_class(sp, "fspace")
  expect_equal(unname(sp$coords), unname(ref$x))
  expect_equal(unname(sp$loadings), unname(ref$rotation))
  expect_equal(sp$eig, ref$sdev^2)
  expect_identical(rownames(sp$coords), rownames(x))
})

test_that("PCA respects scale = FALSE", {
  x <- toy_traits()
  sp <- fs_space(x, method = "pca", scale = FALSE)
  ref <- prcomp(as.matrix(x), center = TRUE, scale. = FALSE)
  expect_equal(unname(sp$coords), unname(ref$x))
})

test_that("PCoA reproduces cmdscale with Cailliez correction", {
  x <- toy_traits()
  d <- dist(scale(as.matrix(x)))
  sp <- fs_space(x, method = "pcoa", correction = "cailliez")
  ref <- cmdscale(d, k = nrow(x) - 1L, eig = TRUE, add = TRUE)
  keep <- which(ref$eig > sqrt(.Machine$double.eps))
  keep <- keep[keep <= ncol(ref$points)]
  expect_equal(unname(sp$coords), unname(ref$points[, keep, drop = FALSE]))
  expect_true(all(sp$eig > 0))
})

test_that("PCoA accepts a user-supplied dist and mixed traits use Gower", {
  x <- toy_traits()
  d <- dist(as.matrix(x))
  sp <- fs_space(x, method = "pcoa", dist = d, correction = "none")
  expect_s3_class(sp, "fspace")
  xm <- x
  xm$cat <- factor(rep(letters[1:3], length.out = nrow(x)))
  sp2 <- fs_space(xm, method = "pcoa")
  expect_s3_class(sp2, "fspace")
  expect_identical(rownames(sp2$coords), rownames(xm))
})

test_that("NMDS fits at default k and stores stress and dist", {
  x <- toy_traits()
  sp <- fs_space(x, method = "nmds", seed = 42L)
  expect_identical(ncol(sp$coords), 4L)
  expect_true(is.numeric(sp$stress) && sp$stress >= 0)
  expect_s3_class(sp$dist, "dist")
})

test_that("raw method returns scaled traits", {
  x <- toy_traits()
  sp <- fs_space(x, method = "raw")
  expect_equal(unname(sp$coords),
               unname(scale(as.matrix(x))),
               ignore_attr = TRUE)
})

test_that("input validation catches bad inputs", {
  x <- toy_traits()
  bad <- x
  rownames(bad) <- NULL
  expect_error(fs_space(bad), "row names")
  expect_error(fs_space(x[, 1, drop = FALSE]), "two traits")
  xm <- x
  xm$cat <- factor("a")
  expect_error(fs_space(xm, method = "pca"), "Non-numeric")
  xna <- x
  xna[1, 1] <- NA
  expect_error(fs_space(xna, method = "pca"), "Missing")
})

test_that("print and summary run silently and return invisibly", {
  sp <- fs_space(toy_traits())
  expect_output(print(sp), "fspace")
  expect_invisible(print(sp))
  expect_output(summary(sp), "TPDs")
})
