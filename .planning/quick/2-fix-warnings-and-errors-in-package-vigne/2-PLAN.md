---
phase: quick
plan: 2
type: execute
wave: 1
depends_on: []
files_modified:
  - vignettes/analyse2.Rmd
  - inst/WORDLIST
autonomous: true
must_haves:
  truths:
    - "All vignettes build without warnings from vignette content issues"
    - "VignetteIndexEntry titles match YAML titles in all vignettes"
    - "Spelling check passes with no new flagged words"
  artifacts:
    - path: "vignettes/analyse2.Rmd"
      provides: "Corrected VignetteIndexEntry title"
      contains: "VignetteIndexEntry{Storing and Analyzing Imputed Data with rbmiUtils}"
    - path: "inst/WORDLIST"
      provides: "Updated spelling wordlist"
      contains: "gtsummary"
  key_links: []
---

<objective>
Fix vignette warnings: correct title mismatch in analyse2.Rmd and update spelling wordlist with missing words.

Purpose: Eliminate warnings that appear during vignette building and R CMD check spelling tests, ensuring clean package checks for CRAN submission.
Output: Warning-free vignette builds and passing spelling checks.
</objective>

<execution_context>
@/Users/bailliem/.claude/get-shit-done/workflows/execute-plan.md
@/Users/bailliem/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@vignettes/analyse2.Rmd (VignetteIndexEntry title mismatch)
@inst/WORDLIST (missing spelling words)
@vignettes/deriving-endpoints.Rmd (source of spelling warnings)
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix VignetteIndexEntry title mismatch and update spelling wordlist</name>
  <files>vignettes/analyse2.Rmd, inst/WORDLIST</files>
  <action>
    1. In vignettes/analyse2.Rmd, change line 10 from:
       `%\VignetteIndexEntry{Storing and Analyzing with rbmiUtils}`
       to:
       `%\VignetteIndexEntry{Storing and Analyzing Imputed Data with rbmiUtils}`
       This makes the VignetteIndexEntry match the YAML title exactly.

    2. Add three words to inst/WORDLIST (maintaining alphabetical order):
       - gtsummary (R package name, used in deriving-endpoints.Rmd:244)
       - standardised (British English spelling, used in deriving-endpoints.Rmd:243)
       - unblinding (clinical trial term, used in deriving-endpoints.Rmd:270)

    Note: The pandoc `--highlight-style` deprecation warnings are an upstream rmarkdown/pandoc compatibility issue and cannot be fixed in the vignette source files. These do not appear in R CMD check results.
  </action>
  <verify>
    Run: `Rscript -e 'rmarkdown::render("vignettes/analyse2.Rmd", output_dir=tempdir(), quiet=FALSE)' 2>&1 | grep -i "WARNING.*VignetteIndexEntry"` -- should produce no output.
    Run: `Rscript -e 'spelling::spell_check_test(vignettes=TRUE, error=FALSE)' 2>&1` -- should show "All Done!" with no potential spelling errors.
  </verify>
  <done>
    VignetteIndexEntry in analyse2.Rmd matches YAML title. Spelling wordlist includes gtsummary, standardised, and unblinding. No content-level warnings remain during vignette building.
  </done>
</task>

</tasks>

<verification>
- `Rscript -e 'devtools::check(vignettes=TRUE, args="--no-tests --no-examples")' 2>&1 | tail -5` shows 0 errors, 0 warnings
- All six vignettes render without content-level warnings
</verification>

<success_criteria>
- analyse2.Rmd VignetteIndexEntry title matches YAML title
- inst/WORDLIST contains gtsummary, standardised, unblinding in alphabetical order
- Spelling check reports no potential spelling errors
</success_criteria>

<output>
After completion, create `.planning/quick/2-fix-warnings-and-errors-in-package-vigne/2-SUMMARY.md`
</output>
