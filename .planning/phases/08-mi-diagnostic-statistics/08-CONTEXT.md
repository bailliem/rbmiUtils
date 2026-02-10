# Phase 8: MI Diagnostic Statistics - Context

**Gathered:** 2026-02-10
**Status:** Ready for planning

<domain>
## Phase Boundary

Enrich ARD output from `pool_to_ard()` with MI-specific diagnostic metadata -- FMI, lambda, degrees of freedom, and relative efficiency -- so regulatory reviewers can assess imputation quality directly from the ARD without manual recomputation. No new accessor functions; users filter the ARD data frame directly. Backward compatibility with existing `pool_to_ard()` usage is required.

</domain>

<decisions>
## Implementation Decisions

### Stat naming & grouping
- Claude's discretion on stat_name conventions (naming style, notation) -- should align with cards/ARD ecosystem and established R MI packages (mice, mitools)
- Claude's discretion on whether stats appear inline with parameter rows or in a separate group -- follow ARD structural conventions
- Claude's discretion on stat_label inclusion and per-visit vs summary diagnostics -- determine based on rbmi pool object structure
- No convenience accessor function -- users filter the ARD data frame directly (e.g., dplyr::filter)
- No specific regulatory format alignment required (not CDISC-constrained) -- just ensure stats are present and clearly labeled
- Claude's discretion on opt-in flag vs always-include behavior -- base on backward compatibility requirements from success criteria

### Variance decomposition depth
- Curated essentials: FMI, lambda, and Barnard-Rubin adjusted df (not full V_w/V_b/V_t breakdown)
- Include relative efficiency (RE = 1 / (1 + FMI/M)) alongside the essentials -- immediately tells reviewers if M was sufficient
- Claude's discretion on classical vs adjusted FMI version -- follow MI literature best practice
- Barnard-Rubin small-sample adjustment for degrees of freedom (user-specified)

### Non-Rubin's-rules handling
- Omit diagnostic rows entirely for non-Rubin pooling methods -- no NA rows, cleaner ARD
- Emit informative message (cli::cli_inform) when diagnostics are omitted due to non-Rubin pooling
- Auto-detect pooling method from pool object structure -- zero user burden, no method= parameter
- Claude's discretion on which rbmi methods qualify as Rubin's-based (frequentist and/or Bayesian)

### Output presentation
- Claude's discretion on print layout (integrated vs separate section) -- follow existing print.ard patterns
- Claude's discretion on decimal precision for diagnostic stats -- follow regulatory reporting norms
- Degrees of freedom displayed as decimal (e.g., 42.7), not rounded integer -- Barnard-Rubin df is non-integer
- Claude's discretion on efficacy_table() compatibility boundary -- determine based on phase success criteria

### Claude's Discretion
- stat_name naming convention (short lowercase vs descriptive vs statistical notation)
- Grouping approach (inline with parameters vs separate diagnostic section)
- Whether stat_label column is included
- Per-visit vs summary diagnostics for multi-visit parameters
- Opt-in flag vs automatic inclusion when analysis_obj provided
- FMI version (classical vs adjusted)
- Print layout for enriched ARD
- Decimal precision for FMI/lambda/RE display
- Whether enriched ARD must pass through efficacy_table() without modification

</decisions>

<specifics>
## Specific Ideas

- Curated essentials (FMI, lambda, Barnard-Rubin df, RE) rather than full variance decomposition -- reviewers care about these, not raw V_w/V_b/V_t
- Clean omission for non-Rubin methods (no NA rows) with informative message -- don't pollute the ARD
- Auto-detection of pooling method from pool object -- friction-free for the user
- df as decimal to preserve Barnard-Rubin precision

</specifics>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 08-mi-diagnostic-statistics*
*Context gathered: 2026-02-10*
