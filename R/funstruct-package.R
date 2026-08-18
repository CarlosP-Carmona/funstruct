#' funstruct: Functional Structure in Trait Spaces
#'
#' A unified workflow for trait-based ecology, built around three classes:
#' `fspace` (a trait space), `fstructure` (functional structure of
#' assemblages) and `fpartition` (diversity partitioned across scales).
#'
#' The intended pipeline is:
#'
#' 1. Build the full space: [fs_space()] (or import one with [as_fspace()]).
#' 2. Evaluate it: `fs_dimensionality()`, `fs_quality()`, `fs_adequacy()`.
#' 3. Reduce (and optionally rotate): [fs_reduce()], [fs_rotate()].
#' 4. Estimate trait probability densities: `fs_tpd()`.
#' 5. Estimate structure and partition diversity: `fs_structure()`,
#'    `fs_partition()`, `fs_beta()`.
#'
#' funstruct supersedes the 'TPD' and 'funspace' packages.
#'
#' @keywords internal
"_PACKAGE"
