# Generate documentation images for README and help pages
# Run: source("data-raw/generate-doc-images.R")
#
# This script creates all pre-rendered PNG images used in:
# - README.Rmd (forest plot and efficacy table visual teasers)
# - plot_forest() help page (example output image)
# - efficacy_table() help page (example output image)
#
# Requires: rbmiUtils, rbmi, dplyr, ggplot2, patchwork, gt, webshot2, chromote
#
# Note: gt::gtsave() for PNG output requires a Chromium-based browser.
# If Google Chrome is not installed, set CHROMOTE_CHROME to an alternative:
#   Sys.setenv(CHROMOTE_CHROME = "/path/to/chromium-browser")

library(rbmiUtils)
library(rbmi)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gt)

# --- Data preparation (same as README.Rmd) ---
data("ADMI", package = "rbmiUtils")

ADMI <- ADMI %>%
  mutate(
    TRT = factor(TRT, levels = c("Placebo", "Drug A")),
    USUBJID = factor(USUBJID),
    AVISIT = factor(AVISIT)
  )

# --- Analysis pipeline ---
vars <- set_vars(
  subjid = "USUBJID",
  visit = "AVISIT",
  group = "TRT",
  outcome = "CHG",
  covariates = c("BASE", "STRATA", "REGION")
)

method <- rbmi::method_bayes(
  n_samples = 100,
  control = rbmi::control_bayes(
    warmup = 200,
    thin = 5
  )
)

ana_obj <- analyse_mi_data(
  data = ADMI,
  vars = vars,
  method = method,
  fun = ancova
)

pool_obj <- pool(ana_obj)

# --- Generate images ---

# Ensure output directory exists
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

# 1. README forest plot
p_forest <- plot_forest(
  pool_obj,
  title = "Treatment Effect: Change from Baseline",
  arm_labels = c(ref = "Placebo", alt = "Drug A")
)
ggsave(
  "man/figures/README-forest-plot-1.png",
  plot = p_forest,
  width = 10, height = 4.5, dpi = 150
)

# 2. README efficacy table
tbl <- efficacy_table(
  pool_obj,
  title = "Table 14.2.1: ANCOVA of Change from Baseline",
  subtitle = "Mixed Model for Repeated Measures",
  arm_labels = c(ref = "Placebo", alt = "Drug A")
)
gt::gtsave(tbl, "man/figures/README-efficacy-table-1.png", vwidth = 800)

# 3. Help page forest plot (same as README version)
ggsave(
  "man/figures/plot_forest-trt.png",
  plot = p_forest,
  width = 10, height = 4.5, dpi = 150
)

# 4. Help page efficacy table (same as README version)
gt::gtsave(tbl, "man/figures/efficacy_table-example.png", vwidth = 800)

message("All documentation images generated successfully.")
message("Files created:")
message("  man/figures/README-forest-plot-1.png")
message("  man/figures/README-efficacy-table-1.png")
message("  man/figures/plot_forest-trt.png")
message("  man/figures/efficacy_table-example.png")
