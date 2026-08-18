# fs_null, fs_sensitivity, fs_pool, fs_redundancy, fs_itv_contribution --------

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

toy_comm2 <- function() {
  rbind(A = c(spA = 0.6, spB = 0.4, spC = 0),
        B = c(spA = 0.2, spB = 0.3, spC = 0.5))
}

test_that("fs_null attaches SES and p-values to the observed structure", {
  sp <- tpd_space()
  st <- fs_null(sp, toy_comm2(), indices = c("richness", "dispersion"),
                n = 19L, seed = 1L)
  expect_s3_class(st, "fstructure")
  expect_true(all(c("ses_richness", "p_richness",
                    "ses_dispersion", "p_dispersion") %in% colnames(st)))
  expect_true(all(st$p_richness > 0 & st$p_richness <= 1))
  expect_true(all(is.finite(st$ses_dispersion)))
  # observed values match a direct fs_structure call
  direct <- fs_structure(sp, toy_comm2(),
                         indices = c("richness", "dispersion"))
  expect_equal(st$richness, direct$richness)
  nl <- attr(st, "nulls")
  expect_identical(nl$n, 19L)
})

test_that("fs_null shuffle.abund permutes within assemblages", {
  sp <- tpd_space()
  st <- fs_null(sp, toy_comm2(), indices = "dispersion",
                model = "shuffle.abund", n = 19L, seed = 2L)
  expect_true(all(is.finite(st$ses_dispersion)))
})

test_that("fs_sensitivity scans bandwidth factors; richness grows with bw", {
  sp <- tpd_space()
  sen <- fs_sensitivity(sp, toy_comm2(), range = c(0.6, 1.4), steps = 3L,
                        indices = "richness")
  expect_s3_class(sen, "fs_sensitivity")
  expect_identical(nrow(sen$values), 6L)  # 3 factors x 2 assemblages
  ag <- tapply(sen$values$richness, sen$values$factor, mean)
  expect_true(all(diff(ag) > 0))  # wider kernels occupy more volume
  expect_identical(nrow(sen$stability), 1L)
})

test_that("fs_sensitivity requires stored observations", {
  m <- data.frame(tr1 = c(0, 1), tr2 = c(1, 0),
                  row.names = c("u1", "u2"))
  spm <- fs_space(m, method = "raw", scale = FALSE)
  tp <- fs_tpd(spm, sds = c(0.3, 0.3))
  # means route stores single-row "observations", so this works too
  sen <- fs_sensitivity(tp, range = c(0.8, 1.2), steps = 3L,
                        indices = "richness")
  expect_s3_class(sen, "fs_sensitivity")
})

test_that("fs_pool warns, pools, and yields assemblage-level TPDs", {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  comm <- toy_comm2()
  expect_warning(
    expect_message(pool <- fs_pool(sp, ti$ids, comm), "plug-in bandwidth"),
    "pseudoreplication")
  expect_identical(names(pool$tpds$units), c("A", "B"))
  for (a in names(pool$tpds$units)) {
    expect_equal(sum(pool$tpds$units[[a]]$probs), 1, tolerance = 1e-12)
  }
  expect_identical(pool$tpds$route, "pooled")
  expect_identical(pool$bw$attachment, "common")
  # the pooled space feeds fs_beta with an identity community matrix
  idc <- diag(2L)
  dimnames(idc) <- list(c("A", "B"), c("A", "B"))
  b <- fs_beta(pool, idc)
  expect_true(b$dissimilarity["A", "B"] > 0)
})

test_that("fs_redundancy: separated vs identical units", {
  sp <- tpd_space()
  red <- fs_redundancy(sp, q = c(0, 1, 2))
  expect_s3_class(red, "fs_redundancy")
  # three fully separated species: redundancy ~ 0 at all q
  expect_true(all(red$values$redundancy < 0.15))

  m <- data.frame(tr1 = c(0, 0), tr2 = c(0, 0),
                  row.names = c("u1", "u2"))
  spm <- fs_space(m, method = "raw", scale = FALSE)
  suppressWarnings(tpm <- fs_tpd(spm, sds = c(0.5, 0.5)))
  red2 <- fs_redundancy(tpm, q = c(0, 1, 2))
  expect_true(all(abs(red2$values$redundancy - 1) < 1e-6))
})

test_that("fs_itv_contribution is positive for spread-type indices", {
  sp <- tpd_space()
  itc <- fs_itv_contribution(sp, toy_comm2())
  expect_s3_class(itc, "fs_itv_contribution")
  expect_identical(colnames(itc$contribution),
                   c("richness", "dispersion"))
  # individual-level spread must add volume beyond the mean-based TPDs
  expect_true(all(itc$contribution$richness > 0))
  expect_true(all(itc$contribution$richness < 1))
  # and it errors without observation-based TPDs
  m <- data.frame(tr1 = c(0, 1), tr2 = c(1, 0),
                  row.names = c("u1", "u2"))
  spm <- fs_space(m, method = "raw", scale = FALSE)
  tpm <- fs_tpd(spm, sds = c(0.3, 0.3))
  expect_error(fs_itv_contribution(tpm), "individual")
})
