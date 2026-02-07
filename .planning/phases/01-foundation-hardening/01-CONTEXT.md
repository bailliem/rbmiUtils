# Phase 1: Foundation Hardening - Context

**Gathered:** 2026-02-07
**Status:** Ready for planning

<domain>
## Phase Boundary

Make existing rbmiUtils functions robust against edge cases, version drift, and malformed inputs. This covers: fixing tidy_pool_obj() parameter parsing, refactoring analyse_mi_data() to delegate to rbmi::analyse(), hardening gcomp_responder()/gcomp_binary() input validation, and ensuring reduce_imputed_data()/expand_imputed_data() preserve full fidelity through round-trips. No new public API surfaces are added.

</domain>

<decisions>
## Implementation Decisions

### Error messaging style
- Use `cli` package formatting (`cli::cli_abort()`, `cli::cli_warn()`) for all user-facing messages
- `cli` is fine to add as a new dependency
- Use custom error class hierarchy (e.g., `rbmiUtils_error_validation`, `rbmiUtils_error_type`) for programmatic catching via `tryCatch()`
- Claude's Discretion: whether individual errors suggest corrective actions or just state what's wrong

### Breaking change policy
- Package follows semver 0.x.y (pre-1.0) convention — minor bumps can include breaking changes
- `analyse_mi_data()` refactor uses a deprecation cycle via `lifecycle::deprecate_warn()` — old internals available for one version before removal
- `tidy_pool_obj()` parsing fix is treated as a breaking change (even though old output was incorrect) — document in NEWS.md
- Add `lifecycle` as a dependency for deprecation warnings with roxygen2 badge integration

### Validation strictness
- Fail early, fail loud — check everything upfront before computation begins
- Single-arm data passed to two-arm functions errors immediately with a clear message
- Do NOT check minimum sample size — let the model decide if it fits
- Stop at the first validation error (no batching of multiple errors)

### Round-trip guarantees
- Exact reproduction required: `identical(original, roundtripped)` must be TRUE
- All column types, factor levels, attributes, and tibble metadata must survive
- Preserve all custom column classes (including any rbmi-specific classes)
- Built-in integrity check: store a digest on reduce, verify on expand
- If integrity check fails on expand, error with details showing which columns/attributes differ

### Claude's Discretion
- Exact digest algorithm for round-trip verification
- Which specific edge cases warrant corrective hints in error messages vs plain error statements
- Internal architecture of validation helpers (shared vs per-function)

</decisions>

<specifics>
## Specific Ideas

- Error messages should follow cli formatting conventions consistent with tidyverse/pharmaverse ecosystem
- Lifecycle badges on deprecated functions for roxygen2 documentation
- Round-trip integrity verification should be transparent — error message should help user diagnose what went wrong

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 01-foundation-hardening*
*Context gathered: 2026-02-07*
