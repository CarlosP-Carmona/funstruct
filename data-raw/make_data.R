# Build the funstruct package datasets from raw sources.
#
# Run ONCE from the package root, before committing:
#   source("data-raw/make_data.R")
#
# Requires the 'funspace' and 'usethis' packages installed.
# Creates data/*.rda (xz-compressed) for:
#   grassland, grassland_abund, gspff, gspff_missing, gspff_tax, gspff_phylo

# ---- Mediterranean grassland (Valdeloshielos) --------------------------------
grassland <- read.csv("data-raw/grassland_long.csv",
                      stringsAsFactors = FALSE)
grassland$plot <- as.integer(grassland$plot)
grassland$individual <- as.integer(grassland$individual)
stopifnot(!anyNA(grassland), length(unique(grassland$plot)) == 40L)

cov <- read.csv("data-raw/grassland_cover.csv", stringsAsFactors = FALSE)
grassland_abund <- with(cov, tapply(cover, list(plot, species), sum))
grassland_abund[is.na(grassland_abund)] <- 0
grassland_abund <- grassland_abund[
  order(as.integer(rownames(grassland_abund))), , drop = FALSE]
stopifnot(nrow(grassland_abund) == 40L,
          all(abs(rowSums(grassland_abund) - 1) < 1e-6))

# ---- Global spectrum of plant form and function (via funspace) ---------------
e <- new.env()
data(list = c("GSPFF", "GSPFF_missing", "GSPFF_tax", "phylo"),
     package = "funspace", envir = e)
gspff <- e$GSPFF
gspff_missing <- e$GSPFF_missing
gspff_tax <- e$GSPFF_tax
gspff_phylo <- e$phylo
stopifnot(nrow(gspff) == 2630L, nrow(gspff_missing) == 10746L)

usethis::use_data(grassland, grassland_abund, gspff, gspff_missing,
                  gspff_tax, gspff_phylo,
                  overwrite = TRUE, compress = "xz")

cat("Done. data/ now contains the six package datasets.\n")
