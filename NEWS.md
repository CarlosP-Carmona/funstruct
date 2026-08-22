# funstruct 0.0.0.9000

## The TPD level stack

* New `fs_aggregate()`: TPDs are now storable at every level of a
  user-defined nesting (individuals -> species -> communities ->
  habitats -> ... -> biomes) through one closed aggregation operator: a
  weighted mixture of TPDs is a TPD on the same grid, so aggregation can
  be applied repeatedly. Groupings can be a community-style matrix
  (many-to-many, entries = weights) or a named child -> parent vector;
  weighting per level is explicit (`"groups"`, `"equal"`, or a named
  numeric vector such as areas) and recorded in the stored level.
  Probability-mass trimming (`alpha`) happens once, at estimation in
  `fs_tpd()`, never during aggregation.
* Storage layout: `space$tpds$levels` is the level stack;
  `fs_tpd()`/`fs_pool()` create its bottom level `"unit"`. The former
  `space$tpds$units` slot is gone (use
  `space$tpds$levels$unit$tpds`). `print()`/`summary()` show the stack,
  including per-level grid occupancy.
* New `fs_level_weights()`: effective bottom-unit weights of any stored
  level (the product of within-level relative weights along the chain).
* New `fs_get_tpd()` + plot method: extract and plot any stored TPD
  (1D profile or 2D image + contours).
* `fs_structure()` and `fs_beta()` gain a `level` argument (mutually
  exclusive with `comm`): compute on the groups of a stored level.
  `fs_structure(space, comm = )` remains one-shot and stores nothing.
  `fs_beta(pooled, level = "unit")` replaces the identity-community
  trick for `fs_pool()` output.
* `fs_partition()` gains a `levels` argument (stored levels, fine to
  coarse); the `comm` + `hierarchy` route now feeds the exact same
  internal engines.
* Unbalanced hierarchies: at every grouping scale, alpha is now the mean
  of group values weighted by the number of assemblages per group
  (Crist et al. 2003), matching the pooled gamma; balanced designs are
  unchanged.
* Internally, all assemblage TPDs are built by a single mixture operator
  shared by `fs_aggregate()`, `fs_structure()`, `fs_beta()` and
  `fs_redundancy()`; regression tests pin the pre-refactor numerical
  outputs on the grassland data.

## Earlier development

* Initial scaffold: `fspace` class, `fs_space()`, `as_fspace()`, `fs_reduce()`, `fs_rotate()`.
* New `fs_dist()`: Gower-style dissimilarities from multiple trait types
  (quantitative, categorical, distributions from means + sds or raw
  observations, category proportions, periods), every trait bounded 0-1.
* `fs_space()` no longer computes dissimilarities: the `dist` argument is
  gone, and `method = "pcoa"`/`"nmds"` take a `dist` object (typically from
  `fs_dist()`) as the `traits` input. NMDS default dimensionality is now 2.
* PCA on the covariance matrix (`scale = FALSE`) warns when trait variances
  differ noticeably.
* New `fs_means()`: per-unit trait means (with `n` and `sds` attributes)
  from individual observations -- the entry point of the mean-based
  workflow (space from means, individuals projected via `fs_project()` or
  `fs_tpd(obs = , ids = )`).
* `fs_project()` recognises `fs_impute()` output and projects only the
  units with imputed values (fully empirical units are skipped with a
  message; subset the object to override). With `add = TRUE`, units
  already present in the space are skipped with a warning instead of
  erroring; the space's own coordinates win.
* The `imputed` attribute of `fs_impute()` renames its `units` element to
  `imputed`, matching the `complete` element.
* Plot methods no longer save and restore `par()`: restoring it reset the
  coordinate system, so `points()`/`text()` added after `plot()` landed in
  the wrong place. All plots now stay open for annotation.
* `plot.fspace()` gains `asp = 1` as default: both axes are in the same
  units, so one unit measures the same length on x and y and distances in
  the ordination read true. Override with `asp = NA` for free scaling.
* `fs_dimensionality()` gains `method = "end"`: the effective number of
  dimensions (Hill number of order `q` of the eigenvalue shares; q = 2 =
  inverse Simpson, Beccari & Carmona 2024 Nat Commun), PCA spaces only.
* New `fs_trait_dim()`: contribution of each trait (or one target trait)
  to the effective dimensionality, with a permutation null in which the
  target trait is randomized (a random trait adds ~1 effective
  dimension). Supports iterative backward selection of the simplest
  sufficient trait set.
* New `fs_angles()`: pairwise angles between traits in the retained
  space, from their loading vectors across axes (Bueno et al. 2023;
  Beccari & Carmona 2024) -- the multidimensional view of trait
  structure that bivariate projections miss. In the full PCA space the
  angles equal `acos(cor(traits))` exactly. Works for PCA directly and
  for PCoA/NMDS after `fs_loadings()`; heatmap `plot()` method. With
  `profile = TRUE`, angles are recomputed at every dimensionality and
  compared with the full-space angles (correlation + mean absolute
  angular deviation): a trait-side dimensionality diagnostic
  complementing the unit-side `fs_quality()`.
* New `fs_loadings()`: post-hoc trait-axis correlations for PCoA/NMDS
  spaces (Pearson or Spearman; categorical traits reported as eta in a
  `categorical` attribute), stored as `$loadings` so `plot()` draws trait
  arrows. Descriptive only: projection stays PCA-only.
