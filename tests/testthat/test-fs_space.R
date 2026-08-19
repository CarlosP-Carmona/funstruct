# fs_space() ------------------------------------------------------------------

toy_traits <- function(n = 30L, p = 5L, seed = 1L) {
  set.seed(seed)
  m <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("sp", seq_len(n)),
                              paste0("tr", seq_len(p))))
  as.data.frame(m)
}

test_that("PCA space matches princomp and returns true loadings", {
  x <- toy_traits()
  sp <- fs_space(x, method = "pca")
  ref <- princomp(as.matrix(x), cor = TRUE)
  expect_s3_class(sp, "fspace")
  expect_equal(unname(sp$coords), unname(ref$scores))
  expect_equal(unname(sp$eig), unname(ref$sdev^2))
  expect_equal(unname(sp$eigenvectors), unname(unclass(ref$loadings)))
  # loadings = eigenvectors scaled by component sdev (funspace convention)
  expect_equal(unname(sp$loadings),
               unname(t(ref$sdev * t(unclass(ref$loadings)))))
  # squared loadings sum to the eigenvalue of each axis
  expect_equal(unname(colSums(sp$loadings^2)), unname(sp$eig))
  # with cor = TRUE, loadings are exactly the trait-axis correlations
  expect_equal(unname(stats::cor(as.matrix(x), sp$coords)),
               unname(sp$loadings), tolerance = 1e-10)
  expect_identical(rownames(sp$coords), rownames(x))
})

test_that("PCA respects scale = FALSE and requires n > p", {
  x <- toy_traits()
  sp <- fs_space(x, method = "pca", scale = FALSE)
  ref <- princomp(as.matrix(x), cor = FALSE)
  expect_equal(unname(sp$coords), unname(ref$scores))
  wide <- toy_traits(n = 4L, p = 6L)
  expect_error(fs_space(wide, method = "pca"), "more units than traits")
})

test_that("fs_project reproduces the space's own coordinates", {
  x <- toy_traits()
  sp <- fs_reduce(fs_space(x, method = "pca"), 2L)
  pr <- fs_project(sp, x)
  expect_equal(pr, sp$coords, tolerance = 1e-10)
  # appending new units flags them in the metadata
  new <- toy_traits(n = 5L, seed = 99L)
  rownames(new) <- paste0("new", 1:5)
  sp2 <- fs_project(sp, new, add = TRUE)
  expect_identical(nrow(sp2$coords), 35L)
  expect_identical(sum(sp2$units$projected), 5L)
  expect_error(fs_project(sp, x[, 1:3]), "lacks trait")
})

test_that("PCoA reproduces cmdscale with Cailliez correction", {
  x <- toy_traits()
  d <- dist(scale(as.matrix(x)))
  sp <- fs_space(d, method = "pcoa", correction = "cailliez")
  ref <- cmdscale(d, k = nrow(x) - 1L, eig = TRUE, add = TRUE)
  keep <- which(ref$eig > sqrt(.Machine$double.eps))
  keep <- keep[keep <= ncol(ref$points)]
  expect_equal(unname(sp$coords), unname(ref$points[, keep, drop = FALSE]))
  expect_true(all(sp$eig > 0))
  expect_identical(rownames(sp$coords), rownames(x))
})

test_that("pcoa/nmds require a dist and pca/raw refuse one", {
  x <- toy_traits()
  expect_error(fs_space(x, method = "pcoa"), "fs_dist")
  expect_error(fs_space(x, method = "nmds"), "fs_dist")
  d <- dist(scale(as.matrix(x)))
  expect_error(fs_space(d, method = "pca"), "pcoa")
  expect_error(fs_space(d, method = "raw"), "pcoa")
  # the intended workflow, including mixed trait types via fs_dist
  xm <- x
  xm$cat <- factor(rep(letters[1:3], length.out = nrow(x)))
  sp2 <- fs_space(fs_dist(xm), method = "pcoa")
  expect_s3_class(sp2, "fspace")
  expect_identical(rownames(sp2$coords), rownames(xm))
})

test_that("NMDS fits at default k = 2 and stores stress and dist", {
  x <- toy_traits()
  d <- fs_dist(x)
  sp <- fs_space(d, method = "nmds", seed = 42L)
  expect_identical(ncol(sp$coords), 2L)
  expect_true(is.numeric(sp$stress) && sp$stress >= 0)
  expect_s3_class(sp$dist, "dist")
  sp3 <- fs_space(d, method = "nmds", k = 3L, seed = 42L)
  expect_identical(ncol(sp3$coords), 3L)
})

test_that("covariance PCA warns when trait variances differ", {
  x <- toy_traits()
  x$tr1 <- x$tr1 * 10          # inflate one variance
  expect_warning(fs_space(x, method = "pca", scale = FALSE),
                 "covariance")
  # equal variances: no warning
  xs <- as.data.frame(scale(as.matrix(toy_traits())))
  expect_silent(fs_space(xs, method = "pca", scale = FALSE))
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

test_that("fs_project add = TRUE skips units already in the space", {
  x <- toy_traits()
  sp <- fs_reduce(fs_space(x, method = "pca"), 2L)
  new <- toy_traits(n = 5L, seed = 99L)
  rownames(new) <- c("sp1", "sp2", "newA", "newB", "newC")  # 2 clashes
  expect_warning(sp2 <- fs_project(sp, new, add = TRUE),
                 "already present")
  expect_identical(nrow(sp2$coords), 33L)          # only the 3 new ones
  expect_setequal(sp2$units$id[sp2$units$projected],
                  c("newA", "newB", "newC"))
  # the clashing units keep the coordinates the space was built with
  expect_equal(sp2$coords["sp1", ], sp$coords["sp1", ])
  # all duplicates -> space returned unchanged (plus the warning)
  expect_warning(sp3 <- fs_project(sp, x[1:3, ], add = TRUE),
                 "already present")
  expect_identical(nrow(sp3$coords), nrow(sp$coords))
})

test_that("fs_means aggregates observations into unit means", {
  set.seed(5)
  ids <- rep(c("spA", "spB", "spC"), times = c(4L, 3L, 1L))
  obs <- data.frame(t1 = rnorm(8), t2 = rnorm(8))
  m <- fs_means(obs, ids)
  expect_identical(rownames(m), c("spA", "spB", "spC"))
  expect_equal(m["spA", "t1"], mean(obs$t1[1:4]), tolerance = 1e-12)
  expect_identical(attr(m, "n"),
                   c(spA = 4L, spB = 3L, spC = 1L))
  sds <- attr(m, "sds")
  expect_equal(sds["spB", "t2"], sd(obs$t2[5:7]), tolerance = 1e-12)
  expect_true(is.na(sds["spC", "t1"]))     # single observation
  # NAs are ignored within units
  obs$t1[1L] <- NA
  m2 <- fs_means(obs, ids)
  expect_equal(m2["spA", "t1"], mean(obs$t1[2:4]), tolerance = 1e-12)
  # validation
  expect_error(fs_means(obs, ids[-1L]), "one entry per row")
  expect_error(fs_means(data.frame(a = letters[1:3]), c("u", "u", "v")),
               "numeric")
})

test_that("mean-based space + projected observations round-trips", {
  set.seed(6)
  ids <- rep(paste0("sp", 1:6), each = 10L)
  obs <- data.frame(t1 = rnorm(60) + rep(1:6, each = 10L),
                    t2 = rnorm(60),
                    row.names = paste0("i", 1:60))
  trM <- fs_means(obs, ids)
  sp <- fs_space(trM, method = "pca")
  co <- fs_project(sp, obs)
  expect_identical(dim(co), c(60L, 2L))
  # each species' projected individuals average to its mean's position
  agg <- rowsum(co, ids) / 10
  expect_equal(agg[rownames(sp$coords), ], sp$coords, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("print and summary run silently and return invisibly", {
  sp <- fs_space(toy_traits())
  expect_output(print(sp), "fspace")
  expect_invisible(print(sp))
  expect_output(summary(sp), "TPDs")
})
