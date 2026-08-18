# fs_structure(): both engines ------------------------------------------------

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

test_that("prob engine: indices behave sensibly on separated species", {
  sp <- tpd_space()
  comm <- rbind(A = c(spA = 0.5, spB = 0.5, spC = 0),
                B = c(spA = 1, spB = 1, spC = 1))
  st <- fs_structure(sp, comm)
  expect_s3_class(st, "fstructure")
  expect_identical(rownames(st), c("A", "B"))
  # adding a third, separated species increases occupied volume
  expect_gt(st["B", "richness"], st["A", "richness"])
  expect_true(all(st$evenness > 0 & st$evenness <= 1))
  expect_true(all(st$divergence >= 0 & st$divergence <= 1))
  expect_true(all(st$dispersion > 0))
  # overlap-based dissimilarities: separated species are nearly maximally
  # dissimilar
  expect_gt(st["A", "mpd"], 0.9)
  expect_true(all(st$rao >= 0 & st$rao <= 1))
  # no overlap between species -> essentially no redundancy
  expect_lt(st["B", "redundancy"], 0.15)
})

test_that("prob engine: identical units give redundancy ~ 1", {
  m <- data.frame(tr1 = c(0, 0), tr2 = c(0, 0),
                  row.names = c("u1", "u2"))
  spm <- fs_space(m, method = "raw", scale = FALSE)
  expect_warning(tp <- fs_tpd(spm, sds = c(0.5, 0.5)), "zero variation")
  st <- fs_structure(tp)
  expect_equal(st[1L, "redundancy"], 1, tolerance = 1e-6)
})

test_that("prob engine: single-species assemblage", {
  sp <- tpd_space()
  comm <- rbind(onlyA = c(spA = 1, spB = 0, spC = 0))
  st <- fs_structure(sp, comm)
  expect_true(is.na(st$mpd))
  expect_true(is.na(st$originality))
  # CWM of a single species sits near its centroid (0, 0)
  expect_lt(abs(st$cwm_1), 0.3)
  expect_lt(abs(st$cwm_2), 0.3)
})

test_that("long-format comm equals the matrix input", {
  sp <- tpd_space()
  comm <- rbind(A = c(spA = 0.5, spB = 0.5, spC = 0),
                B = c(spA = 1, spB = 1, spC = 1))
  long <- data.frame(
    assemblage = c("A", "A", "B", "B", "B"),
    unit = c("spA", "spB", "spA", "spB", "spC"),
    abundance = c(0.5, 0.5, 1, 1, 1))
  st1 <- fs_structure(sp, comm)
  st2 <- fs_structure(sp, long)
  expect_equal(as.data.frame(st1), as.data.frame(st2))
})

test_that("comm validation catches problems", {
  sp <- tpd_space()
  bad <- rbind(A = c(spA = 1, spX = 1))
  expect_error(fs_structure(sp, bad), "Unknown unit")
  neg <- rbind(A = c(spA = -1, spB = 1, spC = 0))
  expect_error(fs_structure(sp, neg), "non-negative")
  zero <- rbind(A = c(spA = 0, spB = 0, spC = 0))
  expect_error(fs_structure(sp, zero), "zero total abundance")
})

test_that("prob engine requires TPDs", {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  expect_error(fs_structure(sp), "fs_tpd")
})

test_that("points engine: exact values on a two-point assemblage", {
  co <- matrix(c(0, 0, 1, 0), 2L, 2L, byrow = TRUE,
               dimnames = list(c("u1", "u2"), c("A1", "A2")))
  sp <- as_fspace(co)
  st <- fs_structure(sp, engine = "points")
  # distances: d(u1, u2) = 1; equal weights 0.5
  expect_equal(st$rao, 0.5)
  expect_equal(st$mpd, 1)
  expect_equal(st$dispersion, 0.5)
  expect_equal(st$originality, 1)
  expect_equal(st$cwm_1, 0.5)
  expect_equal(st$cwm_2, 0)
  expect_true(is.na(st$richness))   # 2 points cannot span a 2D hull
  expect_true(is.na(st$evenness))   # FEve needs S >= 3
  expect_true(is.na(st$redundancy)) # undefined in the points engine
})

test_that("points engine: hull richness matches geometry::convhulln", {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  st <- fs_structure(sp, engine = "points")
  ref <- geometry::convhulln(sp$coords, options = "FA")$vol
  expect_equal(st$richness, ref, tolerance = 1e-10)
  expect_true(st$evenness >= 0 && st$evenness <= 1)
  expect_true(st$divergence >= 0 && st$divergence <= 1)
})

test_that("points engine FEve is within [0, 1] and reacts to evenness", {
  set.seed(5)
  co <- matrix(rnorm(20L), 10L, 2L,
               dimnames = list(paste0("u", 1:10), c("A1", "A2")))
  sp <- as_fspace(co)
  even <- matrix(rep(1, 10L), 1L, dimnames = list("e", rownames(co)))
  st <- fs_structure(sp, comm = even, engine = "points",
                     indices = "evenness")
  expect_true(st$evenness >= 0 && st$evenness <= 1)
})
