# Phase 6: Documentation - Context

**Gathered:** 2026-02-08
**Status:** Ready for planning

<domain>
## Phase Boundary

End-to-end vignette walking through raw data to rbmi to rbmiUtils to regulatory tables and forest plots. README with rendered visual teasers. Rendered examples in function help pages. NEWS.md with version history. Cross-references to rbmi and beeca documentation. Site structure and branding belong in Phase 7.

</domain>

<decisions>
## Implementation Decisions

### Vignette narrative & scope
- Use package datasets (ADEFF/ADMI) — no external dependencies, always reproducible
- Show full rbmi pipeline (draws/impute/analyse/pool) with explanation so reader understands the complete workflow
- Continuous analysis as primary walkthrough, binary/responder as a shorter appendix section showing how it differs
- Pipeline-focused framing (e.g., "From rbmi Analysis to Regulatory Tables") rather than rbmiUtils-centric title
- Keep both the existing analyse2 vignette and the new end-to-end vignette — link between them
- README.Rmd → README.md (standard R package pattern, executable code examples)

### Rendered output strategy
- README.Rmd renders output so images stay in sync with code
- Rendered table and plot images appear in README as visual teasers

### Cross-referencing approach
- Inline hyperlinks woven into prose (not callout boxes) in vignettes
- Add cross-references to existing vignettes (analyse2) where rbmi/beeca links add value
- Add @seealso sections in roxygen for functions that wrap or depend on rbmi/beeca

### NEWS.md structure
- Organize by version number: v1 milestone = 0.1.0, v2 milestone = 0.2.0
- Group entries with sub-bullets (New features, Bug fixes, Improvements, Breaking changes)
- NEWS.md already exists — reorganize existing content into proper versioned structure

### Claude's Discretion
- Writing tone for vignette (tutorial vs reference — pick what fits the audience)
- Whether to include data prep steps (validate_data, prepare_data_ice) in the pipeline narrative
- Customization depth in vignette (defaults only vs showing key customizations)
- README depth (minimal teaser vs quick-start snippet)
- Rendered output approach for function help pages (static images vs live rendering)
- Cross-reference link targets (pkgdown sites vs CRAN pages)

</decisions>

<specifics>
## Specific Ideas

- Pipeline narrative: raw data → data prep → rbmi (draws/impute/analyse/pool) → rbmiUtils (tidy/format) → efficacy table + forest plot
- Continuous analysis is the main walkthrough; binary/responder gets a concise "how it differs" appendix
- Existing analyse2 vignette stays — new vignette is the getting-started guide, analyse2 is the focused analysis reference
- NEWS.md already has content (development version + older entries 0.1.4–0.1.8) — restructure into 0.1.0 and 0.2.0

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 06-documentation*
*Context gathered: 2026-02-08*
