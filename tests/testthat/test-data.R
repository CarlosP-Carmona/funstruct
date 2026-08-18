# Package datasets (skipped until data-raw/make_data.R has been run) ----------

has_data <- function(name) {
  name %in% data(package = "funstruct")$results[, "Item"]
}

test_that("grassland data are consistent", {
  skip_if_not(has_data("grassland"), "data not built yet")
  data("grassland", package = "funstruct", envir = environment())
  data("grassland_abund", package = "funstruct", envir = environment())
  expect_identical(sort(unique(grassland$plot)), 1:40)
  expect_false(anyNA(grassland))
  expect_true(all(grassland$height > 0))
  expect_true(all(grassland$sla > 0))
  expect_identical(nrow(grassland_abund), 40L)
  expect_equal(unname(rowSums(grassland_abund)), rep(1, 40L),
               tolerance = 1e-6)
  # every measured species appears in the abundance matrix
  expect_true(all(unique(grassland$species) %in%
                    colnames(grassland_abund)))
})

test_that("grassland runs through the core pipeline", {
  skip_if_not(has_data("grassland"), "data not built yet")
  data("grassland", package = "funstruct", envir = environment())
  data("grassland_abund", package = "funstruct", envir = environment())
  tr <- grassland[, c("height", "sla")]
  tr$height <- log10(tr$height)
  tr$sla <- log10(tr$sla)
  rownames(tr) <- paste0("ind", seq_len(nrow(tr)))
  sp <- fs_space(tr, method = "raw")
  tp <- suppressWarnings(fs_tpd(sp, ids = grassland$species))
  st <- fs_structure(tp, grassland_abund[1:5, ])
  expect_identical(nrow(st), 5L)
  expect_true(all(st$richness > 0))
})

test_that("gspff data are consistent", {
  skip_if_not(has_data("gspff"), "data not built yet")
  data("gspff", package = "funstruct", envir = environment())
  data("gspff_tax", package = "funstruct", envir = environment())
  expect_identical(dim(gspff), c(2630L, 6L))
  expect_setequal(colnames(gspff),
                  c("la", "ln", "ph", "sla", "ssd", "sm"))
  expect_identical(nrow(gspff_tax), 2630L)
})
