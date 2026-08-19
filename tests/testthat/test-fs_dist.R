# fs_dist() -------------------------------------------------------------------

toy_traits <- function(n = 30L, p = 5L, seed = 1L) {
  set.seed(seed)
  m <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("sp", seq_len(n)),
                              paste0("tr", seq_len(p))))
  as.data.frame(m)
}

test_that("numeric-only fs_dist is classical range-scaled Gower", {
  x <- toy_traits(n = 12L, p = 3L)
  d <- fs_dist(x)
  expect_s3_class(d, "dist")
  expect_true(all(d >= 0 & d <= 1))
  expect_identical(labels(d), rownames(x))
  # manual Gower for one pair
  ranges <- vapply(x, function(v) diff(range(v)), numeric(1L))
  manual <- mean(abs(unlist(x[1, ]) - unlist(x[2, ])) / ranges)
  expect_equal(as.matrix(d)["sp2", "sp1"], manual, tolerance = 1e-12)
  # parity with cluster::daisy where available
  skip_if_not_installed("cluster")
  expect_equal(unname(as.matrix(d)),
               unname(as.matrix(cluster::daisy(x, metric = "gower"))),
               tolerance = 1e-10)
})

test_that("categorical traits contribute 0/1 distances", {
  x <- toy_traits(n = 9L, p = 2L)
  x$cat <- rep(c("a", "b", "c"), each = 3L)
  d <- fs_dist(x)
  pt <- attr(d, "per_trait")
  expect_named(attr(d, "types"), c("tr1", "tr2", "cat"))
  expect_identical(attr(d, "types")[["cat"]], "categorical")
  m <- as.matrix(pt$cat)
  expect_identical(m["sp1", "sp2"], 0)   # same level
  expect_identical(m["sp1", "sp4"], 1)   # different level
})

test_that("NAs are skipped pairwise, Gower style", {
  x <- toy_traits(n = 8L, p = 3L)
  x[1L, 1L] <- NA
  d <- fs_dist(x)
  expect_false(anyNA(as.vector(d)))
  # the sp1-sp2 value is the mean over the two traits both have
  ranges <- vapply(x, function(v) diff(range(v, na.rm = TRUE)), numeric(1L))
  manual <- mean(abs(unlist(x[1L, 2:3]) - unlist(x[2L, 2:3])) / ranges[2:3])
  expect_equal(as.matrix(d)["sp2", "sp1"], manual, tolerance = 1e-12)
  # a pair sharing no trait is an error
  x2 <- toy_traits(n = 5L, p = 2L)
  x2[1L, 1L] <- NA; x2[2L, 2L] <- NA
  expect_error(fs_dist(x2), "no non-missing trait")
})

test_that("gaussian route: overlap of normal distributions", {
  m <- data.frame(t1 = c(0, 0, 100), row.names = c("a", "b", "c"))
  m$t2 <- c(1, 2, 3)  # second trait so the units differ somewhere
  s <- data.frame(t1 = c(1, 1, 1), row.names = c("a", "b", "c"))
  d <- fs_dist(m, sds = s)
  expect_identical(attr(d, "types")[["t1"]], "gaussian")
  pt <- as.matrix(attr(d, "per_trait")$t1)
  expect_lt(pt["b", "a"], 1e-3)     # identical distributions overlap fully
  expect_gt(pt["c", "a"], 0.999)    # 100 sds apart: no overlap
  # zero/missing sds fall back to the mean sd, with a warning
  s2 <- s; s2$t1[2L] <- 0
  expect_warning(fs_dist(m, sds = s2), "mean sd")
})

test_that("kernel route: observations, min_obs fallback, bounds", {
  set.seed(7)
  ids <- rep(c("a", "b", "c"), times = c(20L, 20L, 2L))
  obs <- data.frame(t1 = c(rnorm(20, 0), rnorm(20, 0), rnorm(2, 8)))
  expect_warning(d <- fs_dist(obs = obs, ids = ids), "fewer than")
  expect_identical(attr(d, "types")[["t1"]], "kernel")
  expect_true(all(d >= 0 & d <= 1))
  m <- as.matrix(d)
  expect_lt(m["b", "a"], 0.25)      # same distribution: high overlap
  expect_gt(m["c", "a"], m["b", "a"])
})

test_that("proportions and periods routes", {
  pr <- list(colour = rbind(a = c(1, 0), b = c(1, 0), c = c(0, 1)))
  colnames(pr$colour) <- c("red", "blue")
  pe <- list(flowering = rbind(a = c(100, 200), b = c(100, 200),
                               c = c(300, 60)))
  d <- fs_dist(proportions = pr, periods = pe)
  pc <- as.matrix(attr(d, "per_trait")$colour)
  pf <- as.matrix(attr(d, "per_trait")$flowering)
  expect_identical(pc["b", "a"], 0)
  expect_identical(pc["c", "a"], 1)
  expect_identical(pf["b", "a"], 0)
  expect_identical(pf["c", "a"], 1)   # 300->60 wraps; no overlap with 100-200
  # nested period: fully contained -> distance 0
  pe2 <- list(f = rbind(a = c(100, 200), b = c(120, 150), c = c(1, 50)))
  d2 <- fs_dist(periods = pe2)
  expect_identical(as.matrix(attr(d2, "per_trait")$f)["b", "a"], 0)
  # validation
  bad <- list(colour = rbind(a = c(0.5, 0.4), b = c(1, 0), c = c(0, 1)))
  expect_error(fs_dist(proportions = bad), "sum to 1")
})

test_that("weights reweight the Gower mean", {
  x <- toy_traits(n = 10L, p = 2L)
  d <- fs_dist(x, weights = c(tr1 = 1, tr2 = 0))
  expect_equal(as.vector(d), as.vector(attr(d, "per_trait")$tr1),
               tolerance = 1e-12)
  expect_error(fs_dist(x, weights = c(bad = 1, tr2 = 1)), "named")
})

test_that("input validation", {
  expect_error(fs_dist(), "at least one")
  expect_error(fs_dist(obs = data.frame(t = 1:6)), "ids")
  x <- toy_traits(n = 6L, p = 2L)
  expect_error(fs_dist(sds = x), "requires `traits`")
  # trait name collisions across inputs
  obs <- data.frame(tr1 = rnorm(12))
  ids <- rep(rownames(x), 2L)
  expect_error(fs_dist(x, obs = obs, ids = ids), "repeated")
})

test_that("fs_dist output feeds fs_space directly", {
  x <- toy_traits(n = 15L, p = 3L)
  x$cat <- factor(rep(letters[1:3], 5L))
  sp <- fs_space(fs_dist(x), method = "pcoa")
  expect_s3_class(sp, "fspace")
  expect_identical(rownames(sp$coords), rownames(x))
  expect_true(all(sp$eig > 0))
})
