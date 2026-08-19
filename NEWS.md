# funstruct 0.0.0.9000

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
