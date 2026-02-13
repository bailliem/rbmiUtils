---
phase: quick-1
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - vignettes/pipeline.Rmd
  - R/analysis_utils.R
autonomous: true
must_haves:
  truths:
    - "Vignette builds without Stan compilation warnings in output"
    - "No deprecation warnings triggered during internal function calls"
    - "Deprecated functions remain exported for backwards compatibility"
  artifacts:
    - path: "vignettes/pipeline.Rmd"
      provides: "Suppressed Stan compilation warnings"
      contains: "warning = FALSE"
    - path: "R/analysis_utils.R"
      provides: "Inlined logic replacing deprecated internal calls"
  key_links:
    - from: "R/analysis_utils.R:gcomp_responder"
      to: "inlined extract_covariates logic"
      via: "direct strsplit/unique/trimws calls"
    - from: "R/analysis_utils.R:gcomp_responder"
      to: "inlined as_simple_formula logic"
      via: "direct as.formula/paste0 calls"
---

<objective>
Fix two sources of warnings during vignette build and R CMD check: (1) Stan compilation noise leaking through as warnings in the pipeline vignette, and (2) deprecation warnings from internal calls to extract_covariates2() and as_simple_formula2().

Purpose: Clean warning-free vignette builds and R CMD check output.
Output: Updated vignettes/pipeline.Rmd and R/analysis_utils.R.
</objective>

<execution_context>
@/Users/bailliem/.claude/get-shit-done/workflows/execute-plan.md
@/Users/bailliem/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@vignettes/pipeline.Rmd
@R/analysis_utils.R
</context>

<tasks>

<task type="auto">
  <name>Task 1: Suppress Stan compilation warnings in vignette</name>
  <files>vignettes/pipeline.Rmd</files>
  <action>
    In vignettes/pipeline.Rmd, line 158, change the chunk header from:
    ```{r draws, message = FALSE}
    to:
    ```{r draws, message = FALSE, warning = FALSE}

    Stan compilation emits messages that come through as R warnings, so both message=FALSE and warning=FALSE are needed.
  </action>
  <verify>Grep the file to confirm the chunk has both `message = FALSE` and `warning = FALSE`.</verify>
  <done>The `draws` chunk in pipeline.Rmd has `warning = FALSE` added to its options.</done>
</task>

<task type="auto">
  <name>Task 2: Inline deprecated helper logic in gcomp_responder</name>
  <files>R/analysis_utils.R</files>
  <action>
    In R/analysis_utils.R, replace the two internal calls to deprecated functions inside `gcomp_responder()` with inlined equivalents:

    1. Replace `extract_covariates2(covariates)` (used at lines 72 and 110) with a local helper defined at the top of `gcomp_responder()`:

       ```r
       .extract_covars <- function(x) {
         if (is.null(x)) return(x)
         unique(trimws(unlist(strsplit(x, ":|\\*"))))
       }
       ```

       Then change line 72 from:
         `required_cols <- unique(c(group, outcome, extract_covariates2(covariates)))`
       to:
         `required_cols <- unique(c(group, outcome, .extract_covars(covariates)))`

       And change the extract_covariates2 call at line 110 similarly to use `.extract_covars(covariates)`.

    2. Replace `as_simple_formula2(outcome, ...)` (line 108-111) with inline code:

       ```r
       covars_for_model <- setdiff(unique(c(group, .extract_covars(covariates))), visit)
       frm <- stats::as.formula(
         paste0(outcome, " ~ 1 + ", paste0(covars_for_model, collapse = " + "))
       )
       environment(frm) <- globalenv()
       ```

    3. Do NOT remove or modify the exported `extract_covariates2()` and `as_simple_formula2()` functions -- they remain for backwards compatibility of any external callers.

    4. Do NOT change any other logic in gcomp_responder().
  </action>
  <verify>
    Run `Rscript -e "devtools::load_all(); suppressMessages(library(testthat)); devtools::test()"` to confirm all existing tests pass.
    Also run `Rscript -e "devtools::check(args = '--no-examples')"` or at minimum `Rscript -e "devtools::load_all()"` to confirm no load errors.
  </verify>
  <done>
    gcomp_responder() no longer calls extract_covariates2() or as_simple_formula2() internally.
    The deprecated functions still exist as exported functions.
    All existing tests pass.
  </done>
</task>

</tasks>

<verification>
- `grep -n "warning = FALSE" vignettes/pipeline.Rmd` shows the draws chunk
- `grep -n "extract_covariates2\|as_simple_formula2" R/analysis_utils.R` shows only the function definitions, NOT calls from gcomp_responder
- `devtools::test()` passes
- `devtools::load_all()` succeeds without deprecation warnings
</verification>

<success_criteria>
- The pipeline.Rmd draws chunk suppresses both messages and warnings
- gcomp_responder() uses inlined logic instead of deprecated helpers
- All tests pass, package loads cleanly
</success_criteria>

<output>
After completion, create `.planning/quick/1-fix-vignette-warnings-suppress-stan-comp/1-SUMMARY.md`
</output>
