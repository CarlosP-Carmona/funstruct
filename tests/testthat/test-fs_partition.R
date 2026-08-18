# fs_partition(): equivalent numbers and Rao decomposition --------------------

toy_individuals <- function(n_per = 40L, seed = 11L) {
  set.seed(seed)
  centers <- rbind(spA = c(0, 0), spB = c(4, 0), spC = c(0, 4))
  X <- do.call(rbind, lapply(rownames(centers), function(s) {
    sweep(matrix(rnorm(n_per * 2L, sd = 0.7), n_per, 2L), 2L,
          centers[s, ], "+")
  }))
  ids <- rep(rownames(centers), each = n_per)
  rownames(X) <- paste0("ind", seq_len(nrow(X)))
  colnames(X) <- c("tr1", "tr2")
  list(traits = as.data.frame(X), ids = ids)
}

tpd_space <- function() {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  fs_tpd(sp, ids = ti$ids)
}

toy_comm <- function() {
  rbind(A1 = c(spA = 1, spB = 0, spC = 0),
        A2 = c(spA = 0, spB = 1, spC = 0),
        A3 = c(spA = 0, spB = 0, spC = 1),
        A4 = c(spA = 1, spB = 1, spC = 1))
}

test_that("fs_hierarchy validates and names its levels", {
  h <- fs_hierarchy(site = c("s1", "s1", "s2", "s2"),
                    region = c("r1", "r1", "r1", "r1"))
  expect_s3_class(h, "fs_hierarchy")
  expect_identical(colnames(h), c("site", "region"))
  expect_error(fs_hierarchy(site = c("a", "b"), region = "x"),
               "same length")
})

test_that("tpd_eqv: elementary properties hold", {
  sp <- tpd_space()
  comm <- toy_comm()
  pt <- fs_partition(sp, comm, method = "tpd_eqv", q = 0)
  eqv <- pt$values$assemblage
  # single-unit assemblages have equivalent number 1
  expect_equal(unname(eqv[c("A1", "A2", "A3")]), rep(1, 3L),
               tolerance = 1e-8)
  # three equally abundant, fully separated units -> close to 3
  expect_gt(eqv[["A4"]], 2.6)
  expect_lt(eqv[["A4"]], 3 + 1e-8)
  # gamma equals the union of everything (~3) and components sum to it
  expect_gt(pt$values$total, 2.6)
  expect_equal(sum(pt$table$value), pt$values$total, tolerance = 1e-10)
  expect_true(all(pt$table$value >= -1e-10))
})

test_that("tpd_eqv: identical units yield an equivalent number of 1", {
  m <- data.frame(tr1 = c(0, 0), tr2 = c(0, 0),
                  row.names = c("u1", "u2"))
  spm <- fs_space(m, method = "raw", scale = FALSE)
  expect_warning(tp <- fs_tpd(spm, sds = c(0.5, 0.5)), "zero variation")
  comm <- rbind(A = c(u1 = 1, u2 = 1))
  pt <- fs_partition(tp, comm, method = "tpd_eqv", q = 0)
  expect_equal(pt$values$assemblage[["A"]], 1, tolerance = 1e-8)
})

test_that("tpd_eqv: dominance reduces the equivalent number as q grows", {
  sp <- tpd_space()
  comm <- rbind(U = c(spA = 0.9, spB = 0.05, spC = 0.05))
  e0 <- fs_partition(sp, comm, method = "tpd_eqv", q = 0)
  e2 <- fs_partition(sp, comm, method = "tpd_eqv", q = 2)
  expect_gt(e0$values$assemblage[["U"]], e2$values$assemblage[["U"]])
  # and the profile is non-increasing in q
  pr <- fs_partition(sp, comm, method = "tpd_eqv", profile = TRUE)
  z <- pr$profiles[pr$profiles$scale == "assemblage", ]
  expect_true(all(diff(z$value) <= 1e-8))
})

test_that("tpd_eqv: hierarchical partition is additive and non-negative", {
  sp <- tpd_space()
  comm <- toy_comm()
  h <- fs_hierarchy(site = c("s1", "s1", "s2", "s2"))
  pt <- fs_partition(sp, comm, hierarchy = h, method = "tpd_eqv", q = 0)
  expect_identical(pt$levels, c("assemblage", "site", "total"))
  expect_identical(nrow(pt$table), 3L)
  expect_equal(sum(pt$table$value), pt$values$total, tolerance = 1e-10)
  expect_true(all(pt$table$value >= -1e-10))
  expect_length(pt$values$site, 2L)
})

test_that("crossed hierarchies are rejected", {
  sp <- tpd_space()
  comm <- toy_comm()
  h <- fs_hierarchy(site = c("s1", "s1", "s2", "s2"),
                    region = c("r1", "r2", "r1", "r2"))
  expect_error(fs_partition(sp, comm, hierarchy = h, method = "tpd_eqv"),
               "nested")
})

test_that("rao partition: alpha, beta and corrections behave", {
  sp <- tpd_space()
  comm <- toy_comm()
  pt <- fs_partition(sp, comm, method = "rao")
  # single-species assemblages have zero Rao
  expect_equal(unname(pt$values$assemblage[c("A1", "A2", "A3")]),
               rep(0, 3L), tolerance = 1e-10)
  # pooled gamma is positive and beta components non-negative
  expect_gt(pt$values$total, 0)
  expect_true(all(pt$table$value >= -1e-10))
  expect_equal(sum(pt$table$value), pt$values$total, tolerance = 1e-10)
  # proportional corrected components sum to 1
  expect_equal(sum(pt$table$prop_eqv), 1, tolerance = 1e-10)
})

test_that("rao works without TPDs using scaled Euclidean distances", {
  ti <- toy_individuals(n_per = 5L)
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  comm <- matrix(1, 1L, nrow(sp$coords),
                 dimnames = list("all", rownames(sp$coords)))
  pt <- fs_partition(sp, comm, method = "rao")
  expect_gt(pt$values$assemblage[[1L]], 0)
  expect_lt(pt$values$assemblage[[1L]], 1)
})

test_that("hill is a clear, planned error and prints work", {
  sp <- tpd_space()
  comm <- toy_comm()
  expect_error(fs_partition(sp, comm, method = "hill"),
               "not yet implemented")
  pt <- fs_partition(sp, comm, method = "tpd_eqv")
  expect_output(print(pt), "fpartition")
})
