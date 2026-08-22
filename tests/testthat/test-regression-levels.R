# Regression pins: grassland outputs must never change ------------------------
#
# fs_structure(), fs_beta() (probabilistic engine) and fs_partition() on the
# grassland data are compared against the frozen legacy implementations in
# helper-legacy.R. These tests pin the numerical behaviour at the time the
# fs_aggregate/level-stack refactor started; the refactored code must
# reproduce it exactly (tolerance 1e-12).

data(grassland, envir = environment())
data(grassland_abund, envir = environment())

tr_g <- data.frame(height = log10(grassland$height),
                   sla = log10(grassland$sla),
                   row.names = paste0("ind", seq_len(nrow(grassland))))
sp_g <- fs_space(tr_g, method = "raw")
tp_g <- suppressWarnings(fs_tpd(sp_g, ids = grassland$species))

units_g <- legacy_get_units(tp_g)
grid_g <- tp_g$tpds$grid

strip_fstructure <- function(x) {
  x <- as.data.frame(x)
  attr(x, "engine") <- NULL
  attr(x, "settings") <- NULL
  class(x) <- "data.frame"
  x
}

test_that("fs_structure reproduces the pinned grassland values", {
  all_idx <- c("richness", "evenness", "divergence", "dispersion",
               "rao", "mpd", "originality", "redundancy", "cwm")
  st <- fs_structure(tp_g, grassland_abund)
  W <- legacy_as_comm(grassland_abund, names(units_g), relative = TRUE)
  ref <- legacy_structure_prob(grid_g, units_g, W, all_idx)
  expect_equal(strip_fstructure(st), ref, tolerance = 1e-12)
})

test_that("fs_beta reproduces the pinned grassland values", {
  b <- fs_beta(tp_g, grassland_abund)
  W <- legacy_as_comm(grassland_abund, names(units_g), relative = TRUE)
  ref <- legacy_beta_prob(units_g, W, decompose = TRUE)
  expect_equal(b$dissimilarity, ref$dissimilarity, tolerance = 1e-12)
  expect_equal(b$P_shared, ref$P_shared, tolerance = 1e-12)
  expect_equal(b$P_non_shared, ref$P_non_shared, tolerance = 1e-12)
})

test_that("fs_partition (tpd_eqv) reproduces the pinned grassland values", {
  blocks <- fs_hierarchy(block = rep(paste0("b", 1:4), each = 10L))
  W <- legacy_as_comm(grassland_abund, names(units_g), relative = TRUE)
  H <- as.data.frame(blocks)
  for (qq in c(0, 1)) {
    pt <- fs_partition(tp_g, grassland_abund, hierarchy = blocks,
                       method = "tpd_eqv", q = qq)
    ref <- legacy_eqv_levels(units_g, W, H, qq)
    expect_equal(pt$table$value, ref$table$value, tolerance = 1e-12)
    expect_identical(pt$table$component, ref$table$component)
    expect_equal(pt$values$assemblage, ref$values$assemblage,
                 tolerance = 1e-12)
    expect_equal(pt$values$block, ref$values$block, tolerance = 1e-12)
    expect_equal(pt$values$total, ref$values$total, tolerance = 1e-12)
  }
})

test_that("fs_partition (rao) reproduces the pinned grassland values", {
  blocks <- fs_hierarchy(block = rep(paste0("b", 1:4), each = 10L))
  W <- legacy_as_comm(grassland_abund, names(units_g), relative = TRUE)
  H <- as.data.frame(blocks)
  pt <- fs_partition(tp_g, grassland_abund, hierarchy = blocks,
                     method = "rao")
  ref <- legacy_partition_rao(units_g, W, H)
  expect_equal(pt$table$value, ref$table$value, tolerance = 1e-12)
  expect_equal(pt$table$value_eqv, ref$table$value_eqv, tolerance = 1e-12)
  expect_equal(pt$table$prop_eqv, ref$table$prop_eqv, tolerance = 1e-12)
  expect_identical(pt$table$component, ref$table$component)
})

test_that("null comm (all units, one assemblage) is pinned too", {
  st <- fs_structure(tp_g, comm = NULL,
                     indices = c("richness", "evenness", "redundancy"))
  W <- legacy_as_comm(NULL, names(units_g), relative = TRUE)
  ref <- legacy_structure_prob(grid_g, units_g, W,
                               c("richness", "evenness", "redundancy"))
  expect_equal(strip_fstructure(st), ref, tolerance = 1e-12)
})
