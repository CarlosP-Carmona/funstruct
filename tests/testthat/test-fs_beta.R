# fs_beta() --------------------------------------------------------------------

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

test_that("prob beta: overlap dissimilarity and TPD decomposition", {
  sp <- tpd_space()
  comm <- rbind(A1 = c(spA = 1, spB = 0, spC = 0),
                A2 = c(spA = 0, spB = 1, spC = 0),
                A4 = c(spA = 1, spB = 1, spC = 1))
  b <- fs_beta(sp, comm)
  expect_s3_class(b, "fbeta")
  D <- b$dissimilarity
  expect_equal(unname(diag(D)), rep(0, 3L))
  expect_equal(D, t(D))
  # disjoint single-species assemblages: complete, non-shared dissimilarity
  expect_gt(D["A1", "A2"], 0.99)
  expect_equal(b$P_non_shared["A1", "A2"], 1, tolerance = 1e-6)
  # nested case: A1 subset of A4 -> dissimilarity ~ 2/3, entirely from the
  # shared region (A1 occupies no cells outside A4)
  expect_gt(D["A1", "A4"], 0.55)
  expect_lt(D["A1", "A4"], 0.75)
  expect_equal(b$P_shared["A1", "A4"], 1, tolerance = 1e-6)
  expect_equal(b$P_non_shared["A1", "A4"], 0, tolerance = 1e-6)
})

test_that("prob beta: identical assemblages give zero and NA components", {
  sp <- tpd_space()
  comm <- rbind(X = c(spA = 1, spB = 1, spC = 0),
                Y = c(spA = 1, spB = 1, spC = 0))
  b <- fs_beta(sp, comm)
  expect_equal(b$dissimilarity["X", "Y"], 0, tolerance = 1e-12)
  expect_true(is.na(b$P_shared["X", "Y"]))
})

test_that("points beta: exact Rao DISC on a two-point toy", {
  co <- matrix(c(0, 0, 1, 0), 2L, 2L, byrow = TRUE,
               dimnames = list(c("u1", "u2"), c("A1", "A2")))
  sp <- as_fspace(co)
  comm <- rbind(P = c(u1 = 1, u2 = 0),
                Q = c(u1 = 0, u2 = 1))
  expect_message(b <- fs_beta(sp, comm, engine = "points"),
                 "No decomposition")
  # within-Rao = 0; pooled (0.5, 0.5) with scaled distance 1 -> Q = 0.5
  expect_equal(b$dissimilarity["P", "Q"], 0.5)
  expect_null(b$P_shared)
})

test_that("fs_beta validates inputs", {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  comm1 <- matrix(1, 1L, 3L,
                  dimnames = list("only", c("spA", "spB", "spC")))
  expect_error(fs_beta(sp, comm1, engine = "prob"), "fs_tpd")
  tp <- tpd_space()
  expect_error(fs_beta(tp, comm1), "at least two")
})
