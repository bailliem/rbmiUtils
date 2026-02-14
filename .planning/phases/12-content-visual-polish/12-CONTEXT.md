# Phase 12: Content & Visual Polish - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Deliver a standalone binary responder vignette demonstrating imputed data storage and reanalysis workflow, and refine forest plot visual styling to publication/regulatory submission quality. No new exported functions or features — pure content and visual polish.

</domain>

<decisions>
## Implementation Decisions

### Vignette scope & workflow
- Post-imputation focus: assume rbmi already run, focus on extracting imputed data, defining binary endpoint, and reanalyzing
- Cross-references the pipeline vignette for setup steps, links to it, but stands fully on its own with a pre-built rbmi object
- Uses existing package example data (antidepressant trial data shipped with rbmi/rbmiUtils)
- Demonstrates BOTH endpoint definitions: threshold-based (e.g., >=50% improvement) AND clinical cutoff (e.g., HAMD <=10)
- Shows ARD flow: demonstrates storing binary responder results into an ARD object using rbmiUtils helpers
- Ends with a formatted responder table showing rates by treatment arm with CI — no forest plot or continuous comparison needed
- Uses rbmi's built-in analyse/pool for MI pooling — vignette shows WHAT to pass, not HOW pooling works internally
- Endpoint-only analysis (final visit) — no visit-level analysis
- No citations or academic references needed
- No comparison with continuous primary analysis
- Include brief caveats on assumptions/limitations of binary responder approach from imputed data

### Vignette naming & structure
- Article named "deriving-endpoints" with title "Deriving Endpoints from Imputed Data" — broader framing
- Explicit prerequisites/setup section with library() calls and data loading upfront
- Inline comments in code chunks explaining key steps
- Assumes reader familiarity with rbmi core concepts (draws, impute, analyse, pool) — no re-explanation

### Vignette tone & audience
- Serves both trial statisticians and R programmers equally
- Practical tutorial tone: "Next, we extract the imputed datasets..."
- Output table is self-explanatory — no interpretation walkthrough needed

### Forest plot typography
- Keep default sans-serif font family
- Slightly larger text sizes than current defaults (text_size=3, point_size=3, base_size=11 all need bumping up)
- Bold visit labels in table panel to distinguish from estimate values
- Bold column headers (slightly larger than data rows) in table panel
- Bold and larger plot title when present
- Keep current p-value format: 3 decimal places, threshold at <0.001
- No legend for significance indicators (filled vs open circles)
- Add descriptive x-axis title like "Treatment Difference (95% CI)"

### Forest plot spacing & layout
- Keep current panel width ratios (table:forest:pvalue = 3:4:1.5)
- Add subtle horizontal lines between visit rows for readability
- Tighter row spacing (compact layout to fit more visits)
- Dashed reference line (at 0 or null value) — replaces current solid line
- Document suggested default dimensions for A4/US Letter regulatory documents
- Precise panel alignment: ensure row baselines align perfectly across all three panels
- Add thin border around the forest panel
- Uniform point sizes (no weight-proportional scaling)

### Claude's Discretion
- Whether to include inline results or suppress output (balance readability vs CRAN timing)
- Vignette build strategy (pre-computed vs live) based on CRAN timing constraints
- Whether to include estimand framework context (include if flows naturally)
- Vignette structure and section headings
- Prose depth (balance explanation with conciseness)
- Vignette length (size based on content needed)
- Estimate decimal format adjustments
- CI whisker line weight adjustments for print quality
- LSM mode color palette (optimize for both screen and print)
- Exact text size bump values for regulatory document quality

</decisions>

<specifics>
## Specific Ideas

- Vignette should demonstrate BOTH threshold-based responder (e.g., >=50% reduction) AND clinical cutoff responder (e.g., HAMD <=10) to show flexibility of imputed data reuse
- Forest plot reference line should be dashed — common convention in regulatory forest plots
- Forest panel should have a visible thin border to clearly delineate from table panel
- Forest plot should suggest default dimensions suitable for regulatory document page sizes (A4/US Letter)

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 12-content-visual-polish*
*Context gathered: 2026-02-14*
