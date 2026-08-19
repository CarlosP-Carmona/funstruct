# Evaluation module: fs_quality, fs_dimensionality, fs_adequacy,
# fs_compare, fs_itv ----------------------------------------------------------

toy_traits <- function(n = 30L, p = 5L, seed = 1L) {
  set.seed(seed)
  m <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("sp", seq_len(n)),
                              paste0("tr", seq_len(p))))
  as.data.frame(m)
}

# traits with a clear 2-factor structure, for dimensionality criteria
structured_traits <- function(n = 40L, seed = 2L) {
  set.seed(seed)
  f1 <- rnorm(n); f2 <- rnorm(n)
  m <- cbind(f1 + rnorm(n, sd = 0.1), f1 + rnorm(n, sd = 0.1),
             f1 + rnorm(n, sd = 0.1), f2 + rnorm(n, sd = 0.1),
             f2 + rnorm(n, sd = 0.1), f2 + rnorm(n, sd = 0.1))
  dimnames(m) <- list(paste0("sp", seq_len(n)), paste0("tr", 1:6))
  as.data.frame(m)
}

test_that("fs_quality: PCA preserves distances perfectly at full dims", {
  sp <- fs_space(toy_traits(), method = "pca")
  q <- fs_quality(sp)
  expect_s3_class(q, "fs_quality")
  expect_identical(q$dims, 1:5)
  expect_lt(q$mSD[5L], 1e-10)
  expect_lte(q$mSD[5L], q$mSD[1L])
  expect_identical(attr(q, "best"), 5L)
})

test_that("fs_quality works on PCoA spaces via the stored dist", {
  sp <- fs_space(fs_dist(toy_traits()), method = "pcoa",
                 correction = "cailliez")
  q <- fs_quality(sp, dims = 1:4)
  expect_identical(nrow(q), 4L)
  expect_true(all(q$mSD >= 0))
})

test_that("fs_quality validates dims and uses the stored Gower dist", {
  sp <- fs_space(toy_traits(), method = "pca")
  expect_error(fs_quality(sp, dims = 99), "between 1 and")
  xm <- toy_traits()
  xm$cat <- factor(rep(letters[1:3], length.out = nrow(xm)))
  spm <- fs_space(fs_dist(xm), method = "pcoa")
  expect_s3_class(fs_quality(spm, dims = 1:3), "fs_quality")
})

test_that("fs_dimensionality auc and elbow return sensible suggestions", {
  sp <- fs_space(structured_traits(), method = "pca")
  d1 <- fs_dimensionality(sp, method = "auc")
  expect_s3_class(d1, "fs_dimensionality")
  expect_true(d1$suggested >= 1L && d1$suggested <= 6L)
  d2 <- fs_dimensionality(sp, method = "elbow")
  expect_true(d2$suggested >= 1L && d2$suggested <= 6L)
})

test_that("fs_dimensionality parallel recovers the 2-factor structure", {
  sp <- fs_space(structured_traits(), method = "pca")
  d <- fs_dimensionality(sp, method = "parallel", n_perm = 99L, seed = 7L)
  expect_true(!is.na(d$suggested))
  expect_true(d$suggested >= 1L && d$suggested <= 3L)
  expect_identical(colnames(d$curve), c("axis", "eigenvalue", "null95"))
})

test_that("fs_dimensionality stress scans NMDS dimensionalities", {
  sp <- fs_space(fs_dist(toy_traits()), method = "nmds", seed = 3L)
  d <- fs_dimensionality(sp, method = "stress", k_max = 3L)
  expect_identical(nrow(d$curve), 3L)
  expect_true(all(is.finite(d$curve$stress)))
  # stress should broadly decrease with k
  expect_lt(d$curve$stress[3L], d$curve$stress[1L])
})

test_that("fs_dimensionality errors informatively where undefined", {
  spn <- fs_space(fs_dist(toy_traits()), method = "nmds", seed = 1L)
  expect_error(fs_dimensionality(spn, method = "elbow"), "eigenvalues")
  spp <- fs_space(toy_traits(), method = "pca")
  expect_error(fs_dimensionality(spp, method = "stress"), "distance matrix")
})

test_that("fs_adequacy applies the Silverman requirements", {
  a <- fs_adequacy(n = 100, d = 2)
  expect_true(a$adequate)
  expect_identical(a$n_required, 19)
  b <- fs_adequacy(n = 50, d = 6)
  expect_false(b$adequate)
  sp <- fs_space(toy_traits(), method = "pca")
  cc <- fs_adequacy(fs_reduce(sp, 3L), n = c(10, 100, 1000))
  expect_identical(cc$d, 3L)
  expect_identical(cc$adequate, c(FALSE, TRUE, TRUE))
  expect_warning(fs_adequacy(n = 1e6, d = 12), "10 dimensions")
})

test_that("fs_compare: a space equals itself and differs from noise", {
  sp <- fs_space(toy_traits(), method = "pca")
  self <- fs_compare(sp, sp)
  expect_lt(self$m2, 1e-10)
  expect_equal(self$correlation, 1, tolerance = 1e-6)

  co <- sp$coords
  set.seed(9)
  co[] <- co[sample(nrow(co)), ]
  scr <- as_fspace(co)
  other <- fs_compare(sp, scr)
  expect_gt(other$m2, 0.05)

  r <- fs_compare(sp, sp, method = "correlation")
  expect_equal(r$correlation, 1, tolerance = 1e-6)
})

test_that("fs_itv: swapping means with themselves changes nothing", {
  x <- toy_traits(n = 15L)
  sp <- fs_space(x, method = "pca")
  obs <- x[c(1, 5, 9), ]
  ids <- rownames(x)[c(1, 5, 9)]
  res <- fs_itv(sp, obs, ids)
  expect_s3_class(res, "fs_itv")
  expect_true(all(res$results$angle < 1e-6))
  expect_true(all(abs(res$results$eig_ratio - 1) < 1e-10))
})

test_that("fs_itv: perturbed observations produce finite, positive effects", {
  x <- toy_traits(n = 15L)
  sp <- fs_space(x, method = "pca")
  set.seed(4)
  obs <- x[rep(1:3, each = 2L), ] + rnorm(2 * 3 * ncol(x), sd = 0.5)
  ids <- rownames(x)[rep(1:3, each = 2L)]
  res <- fs_itv(sp, obs, ids)
  expect_identical(nrow(res$results), 6L * ncol(x))
  expect_true(all(is.finite(res$results$angle)))
  expect_true(all(res$results$angle >= 0 & res$results$angle <= 90))
  expect_true(all(res$results$eig_ratio > 0))
  expect_identical(nrow(res$summary), ncol(x))
})

test_that("fs_itv validates its inputs", {
  x <- toy_traits(n = 10L)
  sp <- fs_space(x, method = "pca")
  expect_error(fs_itv(sp, x[1:2, ], c("spX", "spY")), "absent")
  expect_error(fs_itv(sp, x[1:2, -1], rownames(x)[1:2]), "lacks trait")
  spo <- fs_space(fs_dist(x), method = "pcoa")
  expect_error(fs_itv(spo, x[1:2, ], rownames(x)[1:2]), "PCA")
})
