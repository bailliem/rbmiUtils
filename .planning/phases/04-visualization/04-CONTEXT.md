# Phase 4: Visualization - Context

**Gathered:** 2026-02-08
**Status:** Ready for planning

<domain>
## Phase Boundary

Forest plot function for treatment effects across visits from rbmi pool objects. Returns a ggplot2 object. Includes an aligned table panel alongside the plot. Does not include other visualization types, interactive plots, or Shiny integration.

</domain>

<decisions>
## Implementation Decisions

### Plot structure
- Claude's discretion on orientation (horizontal vs vertical) — standard clinical trial forest plot conventions
- Claude's discretion on whether to show treatment difference only or LSM + difference — but user wants the option via argument
- Claude's discretion on visit ordering approach
- Claude's discretion on multi-endpoint handling (single vs faceted)
- Standard clinical trial forest plot style (regulatory-appropriate, CSR-ready)
- Include a table panel alongside the forest plot showing numeric values
- Claude's discretion on table panel columns (estimate, CI, p-value range)
- Claude's discretion on table panel implementation approach (annotation vs patchwork)

### Visual encoding
- Treatment-coded colors — different colors per arm, colorblind-friendly accessible palette
- Confidence intervals displayed as whisker lines without caps (clean, modern)
- Solid line at zero for reference line (not dashed)
- Uniform filled circles for all point estimates (no shape differentiation by type)

### Data display
- User's choice via argument for content: treatment difference vs LSM by arm
- Default to treatment difference when user doesn't specify
- Visually indicate significance — highlight visits where CI excludes zero
- Claude's discretion on significance highlighting method (filled vs open, bold, etc.)

### Customization API
- Function name: `plot_forest()`
- Primary input: pool object directly — `plot_forest(pool_obj)`, consistent with `efficacy_table()` pattern
- Returns ggplot2 object only — no file saving (user calls `ggsave()`)
- Built-in clinical theme (white background, minimal gridlines) — user can override with `+ theme()`
- Claude's discretion on API surface size (minimal vs moderate arguments)

### Claude's Discretion
- Plot orientation (horizontal vs vertical)
- Which estimates to show by default and panel arrangement
- Visit ordering approach (factor order vs auto-detect)
- Multi-endpoint support (single vs faceted)
- Table panel column selection and implementation approach
- Significance highlighting visual style
- API argument count and naming
- Point size and line weight defaults
- Exact color palette selection (must be colorblind-friendly)

</decisions>

<specifics>
## Specific Ideas

- Standard clinical trial forest plot — like what you'd see in a CSR or FDA submission
- Table panel alongside forest plot with aligned numeric values (common clinical convention)
- `plot_forest()` verb-first naming to group with potential future `plot_*` functions
- Consistent with `efficacy_table(pool_obj)` pattern — pool object as primary input

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 04-visualization*
*Context gathered: 2026-02-08*
