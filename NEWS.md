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
