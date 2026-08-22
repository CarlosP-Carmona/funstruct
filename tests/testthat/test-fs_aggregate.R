# fs_aggregate(): the level stack ---------------------------------------------

# cheap means-route stack: 4 species, 3 plots, 2 blocks, 1 region
toy_stack <- function() {
  m <- data.frame(t1 = c(0, 2, 4, 6), t2 = c(0, 1, 0, 1),
                  row.names = paste0("sp", 1:4))
  sp <- fs_space(m, method = "raw", scale = FALSE)
  g <- fs_grid(sp, res = 12L)
  fs_tpd(sp, sds = 0.6, grid = g, alpha = 1)
}

comm_toy <- rbind(A = c(sp1 = 2, sp2 = 1, sp3 = 0, sp4 = 0),
                  B = c(sp1 = 0, sp2 = 1, sp3 = 1, sp4 = 0),
                  C = c(sp1 = 0, sp2 = 0, sp3 = 1, sp4 = 3))

strip_st <- function(x) {
  x <- as.data.frame(x)
  attr(x, "engine") <- NULL
  attr(x, "settings") <- NULL
  class(x) <- "data.frame"
  x
}

test_that("fs_tpd stores the bottom level of the stack", {
  tp <- toy_stack()
  expect_named(tp$tpds$levels, "unit")
  lv <- tp$tpds$levels$unit
  expect_named(lv$tpds, paste0("sp", 1:4))
  expect_null(lv$from)
  expect_null(lv$members)
})

test_that("fs_aggregate matrix form: members, weights, mixture TPDs", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  expect_named(tp$tpds$levels, c("unit", "plot"))
  lv <- tp$tpds$levels$plot
  expect_identical(lv$from, "unit")
  expect_identical(lv$weights_rule, "groups")
  expect_named(lv$tpds, c("A", "B", "C"))
  expect_equal(unname(rowSums(lv$members)), rep(1, 3), tolerance = 1e-12)
  expect_equal(unname(lv$members["A", ]), c(2 / 3, 1 / 3, 0, 0),
               tolerance = 1e-12)
  expect_equal(unname(lv$n_children), rep(2, 3))
  for (g in names(lv$tpds)) {
    expect_equal(sum(lv$tpds[[g]]$probs), 1, tolerance = 1e-12)
  }
  # each group TPD is exactly the one-step weighted mixture of its units
  units <- tp$tpds$levels$unit$tpds
  for (g in rownames(lv$members)) {
    ref <- legacy_tpdc_sparse(units, lv$members[g, ])
    expect_identical(lv$tpds[[g]]$cells, ref$cells)
    expect_equal(lv$tpds[[g]]$probs, ref$probs, tolerance = 1e-15)
  }
})

test_that("fs_aggregate matrix form completes missing children with zeros", {
  tp <- fs_aggregate(toy_stack(), comm_toy[, 1:3], name = "plot")
  M <- tp$tpds$levels$plot$members
  expect_identical(colnames(M), paste0("sp", 1:4))
  expect_true(all(M[, "sp4"] == 0))
})

test_that("fs_aggregate equal weights ignore the matrix entries", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot",
                     weights = "equal")
  lv <- tp$tpds$levels$plot
  expect_identical(lv$weights_rule, "equal")
  expect_equal(unname(lv$members["A", ]), c(0.5, 0.5, 0, 0),
               tolerance = 1e-12)
})

test_that("fs_aggregate vector form: nesting, user weights, defaults", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  # every child must be assigned
  expect_error(fs_aggregate(tp, c(A = "b1"), name = "block"),
               "assign every child")
  blk <- c(A = "b1", B = "b1", C = "b2")
  area <- c(A = 1, B = 3, C = 5)
  tp <- fs_aggregate(tp, blk, name = "block", weights = area)
  lv <- tp$tpds$levels$block
  expect_identical(lv$from, "plot")   # default from = top level
  expect_identical(lv$weights_rule, "user")
  expect_equal(lv$members["b1", c("A", "B")], c(A = 0.25, B = 0.75),
               tolerance = 1e-12)
  expect_equal(unname(lv$members["b2", ]), c(0, 0, 1), tolerance = 1e-12)
  # vector form without weights: equal, and recorded as such
  tp <- fs_aggregate(tp, c(b1 = "r1", b2 = "r1"), name = "region")
  expect_identical(tp$tpds$levels$region$weights_rule, "equal")
  expect_equal(unname(tp$tpds$levels$region$members), matrix(0.5, 1, 2),
               tolerance = 1e-12)
  # factor groups: row order follows the factor levels
  tp2 <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  fblk <- factor(c(A = "b2", B = "b2", C = "b1"), levels = c("b1", "b2"))
  names(fblk) <- c("A", "B", "C")
  tp2 <- fs_aggregate(tp2, fblk, name = "block")
  expect_identical(rownames(tp2$tpds$levels$block$members), c("b1", "b2"))
})

test_that("fs_aggregate validates inputs", {
  m <- data.frame(t1 = c(0, 2), t2 = c(0, 1), row.names = c("u1", "u2"))
  sp <- fs_space(m, method = "raw", scale = FALSE)
  expect_error(fs_aggregate(sp, comm_toy, name = "plot"), "fs_tpd")
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  expect_error(fs_aggregate(tp, comm_toy, name = "plot"),
               "already exists")
  expect_error(fs_aggregate(tp, comm_toy, name = "x", from = "nope"),
               "stored level")
  bad <- comm_toy
  colnames(bad)[1] <- "spX"
  expect_error(fs_aggregate(toy_stack(), bad, name = "plot"), "Unknown")
  expect_error(fs_aggregate(toy_stack(), comm_toy, name = "plot",
                            weights = c(sp1 = 1)), "vector form")
  expect_error(fs_aggregate(toy_stack(), comm_toy, name = "plot",
                            weights = "nonsense"), "weights")
})

test_that("fs_level_weights composes as the product across levels", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  tp <- fs_aggregate(tp, c(A = "b1", B = "b1", C = "b2"), name = "block",
                     weights = c(A = 1, B = 3, C = 5))
  tp <- fs_aggregate(tp, c(b1 = "r1", b2 = "r1"), name = "region")

  Wb <- fs_level_weights(tp, "block", to = "unit")
  ref <- tp$tpds$levels$block$members %*% tp$tpds$levels$plot$members
  expect_equal(Wb, ref, tolerance = 1e-15)
  Wr <- fs_level_weights(tp, "region", to = "unit")
  expect_equal(unname(rowSums(Wr)), 1, tolerance = 1e-12)
  # identity at the level itself
  I <- fs_level_weights(tp, "unit")
  expect_equal(unname(I), diag(1, 4L))
  expect_identical(rownames(I), paste0("sp", 1:4))
  # intermediate target
  Wrb <- fs_level_weights(tp, "region", to = "block")
  expect_equal(unname(Wrb), matrix(0.5, 1, 2), tolerance = 1e-12)
  # `to` must lie on the chain
  expect_error(fs_level_weights(tp, "unit", to = "block"),
               "aggregation chain")
})

test_that("multi-step aggregation equals the effective one-step mixture", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  tp <- fs_aggregate(tp, c(A = "b1", B = "b1", C = "b2"), name = "block",
                     weights = c(A = 1, B = 3, C = 5))
  tp <- fs_aggregate(tp, c(b1 = "r1", b2 = "r1"), name = "region")
  units <- tp$tpds$levels$unit$tpds
  for (level in c("block", "region")) {
    W <- fs_level_weights(tp, level, to = "unit")
    for (g in rownames(W)) {
      ref <- legacy_tpdc_sparse(units, W[g, ])
      got <- tp$tpds$levels[[level]]$tpds[[g]]
      expect_identical(got$cells, ref$cells)
      expect_equal(got$probs, ref$probs, tolerance = 1e-12)
    }
  }
})

test_that("fs_structure on a stored level matches the comm route", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  st_lvl <- fs_structure(tp, level = "plot")
  st_comm <- fs_structure(tp, comm_toy)
  expect_equal(strip_st(st_lvl), strip_st(st_comm), tolerance = 1e-12)
  # level = "unit": one row per unit
  st_u <- fs_structure(tp, level = "unit", indices = "richness")
  expect_identical(rownames(st_u), paste0("sp", 1:4))
  # errors
  expect_error(fs_structure(tp, comm_toy, level = "plot"), "not both")
  expect_error(fs_structure(tp, level = "plot", engine = "points"),
               "probabilistic")
  expect_error(fs_structure(tp, level = "nope"), "stored level")
})

test_that("fs_beta on a stored level matches the comm route", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  b_lvl <- fs_beta(tp, level = "plot")
  b_comm <- fs_beta(tp, comm_toy)
  expect_equal(b_lvl$dissimilarity, b_comm$dissimilarity,
               tolerance = 1e-12)
  expect_equal(b_lvl$P_shared, b_comm$P_shared, tolerance = 1e-12)
  expect_equal(b_lvl$P_non_shared, b_comm$P_non_shared,
               tolerance = 1e-12)
  expect_error(fs_beta(tp), "Supply")
  expect_error(fs_beta(tp, comm_toy, level = "plot"), "not both")
  expect_error(fs_beta(tp, level = "plot", engine = "points"),
               "probabilistic")
})

test_that("fs_partition on stored levels matches the comm route", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  tp <- fs_aggregate(tp, c(A = "b1", B = "b1", C = "b2"), name = "block")
  h <- fs_hierarchy(block = c("b1", "b1", "b2"))
  for (m in c("rao", "tpd_eqv")) {
    p_comm <- fs_partition(tp, comm_toy, hierarchy = h, method = m)
    p_lvl <- fs_partition(tp, levels = c("unit", "plot", "block"),
                          method = m)
    expect_equal(p_lvl$table$value, p_comm$table$value, tolerance = 1e-12)
    expect_equal(unname(p_lvl$values[["plot"]]),
                 unname(p_comm$values[["assemblage"]]), tolerance = 1e-12)
    expect_equal(p_lvl$values$block, p_comm$values$block,
                 tolerance = 1e-12)
  }
  expect_identical(fs_partition(tp, levels = c("unit", "plot", "block"),
                                method = "rao")$levels,
                   c("plot", "block", "total"))
  # validation
  expect_error(fs_partition(tp, comm_toy,
                            levels = c("unit", "plot")), "not both")
  expect_error(fs_partition(tp, levels = "plot"), "at least two")
  expect_error(fs_partition(tp, levels = c("plot", "unit")),
               "fine to coarse")
  expect_error(fs_partition(tp, method = "rao"), "Supply")
})

test_that("unbalanced hierarchies use size-weighted alpha", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  h <- fs_hierarchy(block = c("b1", "b1", "b2"))   # 2 + 1 assemblages
  for (m in c("rao", "tpd_eqv")) {
    p <- fs_partition(tp, comm_toy, hierarchy = h, method = m)
    g_vals <- p$values$block
    alpha_block <- p$table$value[1L] + p$table$value[2L]
    expect_equal(alpha_block,
                 (2 * g_vals[["b1"]] + g_vals[["b2"]]) / 3,
                 tolerance = 1e-12)
  }
})

test_that("fs_pool output works with fs_beta(level = 'unit')", {
  set.seed(1)
  tr <- data.frame(t1 = rnorm(60, rep(c(0, 3, 6), each = 20), 0.5),
                   t2 = rnorm(60, rep(c(0, 3, 0), each = 20), 0.5),
                   row.names = paste0("i", 1:60))
  ids <- rep(c("a", "b", "c"), each = 20)
  sp <- fs_space(tr, method = "raw", scale = FALSE)
  comm <- rbind(P = c(a = 1, b = 1, c = 0),
                Q = c(a = 0, b = 1, c = 1))
  pooled <- suppressWarnings(suppressMessages(
    fs_pool(sp, ids, comm, grid = fs_grid(sp, res = 15L))))
  expect_named(pooled$tpds$levels, "unit")
  b <- fs_beta(pooled, level = "unit")
  expect_identical(rownames(b$dissimilarity), c("P", "Q"))
  expect_gt(b$dissimilarity["P", "Q"], 0)
})

test_that("fs_get_tpd extracts and plots stored TPDs", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  f <- fs_get_tpd(tp, "plot")
  expect_s3_class(f, "ftpd")
  expect_named(f$tpds, c("A", "B", "C"))
  expect_output(print(f), "plot")
  f1 <- fs_get_tpd(tp, "plot", ids = "B")
  expect_named(f1$tpds, "B")
  expect_error(fs_get_tpd(tp, "plot", ids = "nope"), "Unknown")
  expect_error(fs_get_tpd(tp, "nope"), "stored level")
  grDevices::pdf(NULL)
  plot(f, id = "A")
  grDevices::dev.off()
  # 1D grids plot as density profiles
  m1 <- data.frame(t1 = c(0, 3, 6), t2 = c(0, 1, 0),
                   row.names = paste0("u", 1:3))
  sp1 <- fs_reduce(fs_space(m1, method = "raw", scale = FALSE), 1L)
  tp1 <- fs_tpd(sp1, sds = 0.5, grid = fs_grid(sp1, res = 40L))
  grDevices::pdf(NULL)
  plot(fs_get_tpd(tp1, "unit", "u2"))
  grDevices::dev.off()
})

test_that("summary prints the stack; reducing discards it", {
  tp <- fs_aggregate(toy_stack(), comm_toy, name = "plot")
  tp <- fs_aggregate(tp, c(A = "b1", B = "b1", C = "b2"), name = "block")
  expect_output(print(tp), "unit < plot < block")
  expect_output(summary(tp), "TPD level stack")
  expect_output(summary(tp), "occupancy")
  expect_warning(r <- fs_reduce(tp, 1L), "aggregated levels")
  expect_null(r$tpds)
})
