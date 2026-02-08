# Phase 3: Efficacy Tables - Context

**Gathered:** 2026-02-08
**Status:** Ready for planning

<domain>
## Phase Boundary

Generate regulatory-style efficacy summary tables from rbmi pool objects via gtsummary + gt. Users call a function with a pool object and receive a formatted gt table showing LS means by arm, treatment differences, CIs, and p-values organized by visit. Safety tables, demographic tables, and other table types are out of scope.

</domain>

<decisions>
## Implementation Decisions

### Table structure & layout
- Visits as row groups: each visit is a row group header with statistics beneath
- Within each visit group: show LS mean per arm, then treatment difference row
- Standard CDISC/ICH Table 14.2.x style formatting
- Full header: table title (endpoint name) + subtitle (model description)
- Visit labels auto-cleaned: replace underscores with spaces, title case (e.g., 'Week_4' -> 'Week 4')

### Content & statistics shown
- LS mean rows show: estimate, standard error, 95% CI
- Treatment difference statistics: Claude's discretion (standard regulatory set)
- CI level configurable with 95% as default

### Formatting conventions
- Decimal precision configurable with 2 decimal places as default
- CI bracket style: Claude's discretion (pick sensible default)
- P-value formatting: Claude's discretion (regulatory standard)
- Default footnotes include ALL of: model description, number of imputations, pooling method

### Function interface
- gt only output (no flextable/huxtable)
- Returns gt object with polished defaults, but user can pipe into further gt customization
- gt output suitable for HTML and PDF inclusion in clinical study reports

### Claude's Discretion
- Single function vs two-step (prepare + render) API design
- Whether to accept pool objects only or also ARD tibbles from pool_to_ard()
- Column header organization (unified vs separate spanner headers for LS means and differences)
- N count display approach (per-visit rows, column headers, or omitted)
- Multi-parameter handling (one table per parameter vs nested row groups)
- Treatment difference statistics selection (estimate + SE + CI + p-value set)
- CI source: display from pool object vs allow recalculation
- CI bracket style choice
- P-value decimal places and threshold

</decisions>

<specifics>
## Specific Ideas

- Follow standard CDISC/ICH Table 14.2.x conventions for regulatory submissions
- Visit labels should be human-readable (auto-clean underscores, title case)
- Footnotes should document the full analysis context: model, imputations, pooling method
- Table should be self-contained (title + subtitle + body + footnotes) but gt object returned for power-user customization

</specifics>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 03-efficacy-tables*
*Context gathered: 2026-02-08*
