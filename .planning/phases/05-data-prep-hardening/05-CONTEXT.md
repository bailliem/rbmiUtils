# Phase 5: Data Prep Hardening - Context

**Gathered:** 2026-02-08
**Status:** Ready for planning

<domain>
## Phase Boundary

Make `validate_data()` and `prepare_data_ice()` reject bad input with clear, actionable error messages. Handle edge cases gracefully. No new data prep functions -- harden what exists.

</domain>

<decisions>
## Implementation Decisions

### Error message style
- Use cli package for formatted error messages (bullets, color, cross marks) -- tidyverse-style
- cli is a hard dependency (added to Imports), not conditional
- Claude's discretion on verbosity level per situation (concise vs detailed)
- Claude's discretion on whether to echo user's actual data values in messages
- Claude's discretion on batching multiple warnings vs separate messages

### Validation strictness
- Fail fast on first fatal error -- don't collect all errors before stopping
- Don't check for unused columns in the data -- only validate what the model references
- Claude's discretion on whether malformed interaction terms are errors vs warnings (based on whether they'd produce wrong results)
- Claude's discretion on NULL strategy handling in prepare_data_ice() (error vs warn+default)

### Edge case behavior
- All-NA covariate columns: warn and automatically drop from the model
- All-complete data (no missing outcomes): show informational message confirming no ICE imputations needed
- Claude's discretion on single-subject/single-visit thresholds
- Claude's discretion on all-NA outcome severity

### Coercion policy
- Character visit columns: warn but do NOT auto-coerce -- user must convert to factor themselves (preserves explicit ordering control)
- Character-to-numeric columns: hard error, require correct types upfront -- no auto-coercion
- Claude's discretion on batching type-related warnings

### Claude's Discretion
- Verbosity and data context in error messages (per-situation judgment)
- Batching strategy for warnings vs separate messages
- Severity of malformed interaction terms (error vs warning)
- NULL strategy behavior in prepare_data_ice()
- Single-subject/single-visit handling threshold
- All-NA outcome severity level

</decisions>

<specifics>
## Specific Ideas

- Error messages should feel like tidyverse -- cli-formatted with bullets and cross marks
- Fail-fast philosophy: surface the first fatal problem, don't make users wade through a list
- Character visits are a common mistake in clinical trial data -- the warning should guide users to `factor()` with explicit level ordering
- When all data is complete, the informational message helps users confirm they're calling the right function

</specifics>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 05-data-prep-hardening*
*Context gathered: 2026-02-08*
