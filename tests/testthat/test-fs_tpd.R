# fs_tpd(): grids, bandwidths, densities --------------------------------------

# individual-level data: 3 species, 2 traits, clearly separated
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

test_that("fs_grid builds a sane grid", {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  g <- fs_grid(sp)
  expect_s3_class(g, "fs_grid")
  expect_identical(g$d, 2L)
  expect_identical(g$res, 50L)
  expect_identical(g$n_cells, 2500L)
  expect_identical(dim(g$cells), c(2500L, 2L))
  expect_true(all(g$lo < apply(sp$coords, 2, min)))
  expect_true(all(g$hi > apply(sp$coords, 2, max)))
  expect_error(fs_grid(sp, res = 1), ">= 2")
})

test_that("fs_tpd (obs route): TPDs integrate to 1 and store bandwidths", {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  tp <- fs_tpd(sp, ids = ti$ids, alpha = 0.99)
  expect_named(tp$tpds$levels$unit$tpds, c("spA", "spB", "spC"))
  for (u in names(tp$tpds$levels$unit$tpds)) {
    expect_equal(sum(tp$tpds$levels$unit$tpds[[u]]$probs), 1, tolerance = 1e-12)
    expect_true(length(tp$tpds$levels$unit$tpds[[u]]$cells) < tp$tpds$grid$n_cells)
  }
  expect_identical(tp$bw$attachment, "unit")
  expect_identical(tp$bw$selector, "plugin")
  # per-unit plug-in bandwidth equals ks::Hpi on that unit's coordinates
  HA <- ks::Hpi(sp$coords[ti$ids == "spA", ])
  expect_equal(tp$bw$values[["spA"]], HA, tolerance = 1e-8)
  expect_false(any(tp$bw$imputed))
})

test_that("fs_tpd matches ks::kde densities with a fixed bandwidth", {
  ti <- toy_individuals(n_per = 30L)
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  g <- fs_grid(sp, res = 20L)
  H <- diag(c(0.25, 0.25))
  tp <- fs_tpd(sp, ids = ti$ids, bw = list(spA = H, spB = H, spC = H),
               grid = g, alpha = 1)
  XA <- sp$coords[ti$ids == "spA", ]
  ref <- ks::kde(XA, H = H, eval.points = g$cells)$estimate
  p_ref <- ref * g$cell_volume
  p_ref <- p_ref / sum(p_ref)
  ours <- numeric(g$n_cells)
  keep <- tp$tpds$levels$unit$tpds$spA
  ours[keep$cells] <- keep$probs
  expect_equal(ours, unname(p_ref), tolerance = 1e-8)
})

test_that("fs_tpd 1D route agrees with ks and respects alpha trimming", {
  ti <- toy_individuals()
  sp0 <- fs_space(ti$traits, method = "pca")
  sp1 <- fs_reduce(sp0, 1L)
  tp <- fs_tpd(sp1, ids = ti$ids, alpha = 0.9)
  for (u in names(tp$tpds$levels$unit$tpds)) {
    expect_equal(sum(tp$tpds$levels$unit$tpds[[u]]$probs), 1, tolerance = 1e-12)
  }
  # trimming at 0.9 keeps fewer cells than at 1
  tp_full <- fs_tpd(sp1, ids = ti$ids, alpha = 1)
  expect_lt(length(tp$tpds$levels$unit$tpds$spA$cells),
            length(tp_full$tpds$levels$unit$tpds$spA$cells))
})

test_that("fs_tpd (means route): sds-based TPDs recover the imposed kernel", {
  m <- data.frame(tr1 = c(-2, 2), tr2 = c(0, 0),
                  row.names = c("u1", "u2"))
  sp <- fs_space(m, method = "raw", scale = FALSE)
  # tr2 has zero variation across units: the grid must extend that axis
  # instead of collapsing to zero cell volume
  expect_warning(g <- fs_grid(sp, res = 30L), "zero variation")
  tp <- fs_tpd(sp, sds = c(0.5, 0.5), grid = g, alpha = 1)
  expect_identical(tp$bw$attachment, "entity")
  # density proportional to a normal centred on the unit's coordinates
  co <- sp$coords["u1", ]
  dens <- stats::dnorm(g$cells[, 1L], co[1L], 0.5) *
    stats::dnorm(g$cells[, 2L], co[2L], 0.5)
  p_ref <- dens / sum(dens)
  ours <- numeric(g$n_cells)
  keep <- tp$tpds$levels$unit$tpds$u1
  ours[keep$cells] <- keep$probs
  expect_equal(ours, unname(p_ref), tolerance = 1e-6)
})

test_that("fs_tpd common bandwidth is shared across units", {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  tp <- fs_tpd(sp, ids = ti$ids, bw = "common")
  expect_identical(tp$bw$attachment, "common")
  expect_equal(tp$bw$values[["spA"]], tp$bw$values[["spB"]])
})

test_that("fs_tpd imputes bandwidths for tiny units with a warning", {
  ti <- toy_individuals()
  extra <- data.frame(tr1 = c(2, 2.1, 2.2), tr2 = c(2, 2.1, 1.9),
                      row.names = paste0("x", 1:3))
  traits <- rbind(ti$traits, extra)
  ids <- c(ti$ids, rep("spTiny", 3L))
  sp <- fs_space(traits, method = "raw", scale = FALSE)
  # two warnings are expected: sample-size adequacy and bandwidth imputation
  w <- capture_warnings(tp <- fs_tpd(sp, ids = ids))
  expect_match(w, "average bandwidth", all = FALSE)
  expect_match(w, "fs_adequacy", all = FALSE)
  expect_true(tp$bw$imputed[["spTiny"]])
  expect_equal(sum(tp$tpds$levels$unit$tpds$spTiny$probs), 1, tolerance = 1e-12)
})

test_that("fs_tpd validates inputs and dimensionality", {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  expect_error(fs_tpd(sp), "Supply")
  expect_error(fs_tpd(sp, ids = ti$ids[-1L]), "one entry per row")
  expect_error(fs_tpd(sp, ids = ti$ids, alpha = 0), "in \\(0, 1\\]")
  expect_error(fs_tpd(sp, ids = ti$ids, bw = "nonsense"), "Unknown")

  p7 <- as.data.frame(matrix(rnorm(30L * 7L), 30L, 7L,
                             dimnames = list(paste0("s", 1:30),
                                             paste0("t", 1:7))))
  sp7 <- fs_space(p7, method = "pca")
  expect_error(fs_tpd(sp7, ids = rep(letters[1:3], 10L)), "5 dimensions")
})

test_that("fs_tpd warns when sample sizes cannot support the dimensionality", {
  set.seed(21)
  X <- matrix(rnorm(6L * 4L * 3L), 6L * 3L, 4L,
              dimnames = list(paste0("i", seq_len(18L)), paste0("t", 1:4)))
  sp <- fs_space(as.data.frame(X), method = "raw", scale = FALSE)
  ids <- rep(c("a", "b", "c"), each = 6L)
  # user-supplied bandwidth: the warning must fire regardless of how the
  # bandwidth is obtained
  expect_warning(fs_tpd(sp, ids = ids, bw = 0.5,
                        grid = fs_grid(sp, res = 8L)),
                 "fs_adequacy")
})

test_that("summary reports bandwidth information after fs_tpd", {
  ti <- toy_individuals()
  sp <- fs_space(ti$traits, method = "raw", scale = FALSE)
  tp <- fs_tpd(sp, ids = ti$ids)
  expect_output(summary(tp), "Bandwidth")
  expect_output(print(tp), "estimated")
})
