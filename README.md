# Plasma Specificity Scanning (PSS)

Plasma Specificity Scanning (PSS) is a neutralization fingerprinting method
developed to delineate multi-specific broadly neutralizing antibody (bnAb)
activity in plasma. It differentiates plasma responses driven by a single bnAb
specificity from those that develop multiple bnAb specificities simultaneously.

PSS uses the same data as the Spearman-based MSBP method: plasma neutralization
data against the 41-virus multiclade panel, compared with a reference bnAb set.
First, the reference bnAb data is used to identify sub-panels of the reference
pseudovirus panel (of varying sizes) that best classify reference bnAbs into
their known epitope clusters. The best-performing virus sub-panels (n = 100)
are then used to predict plasma specificity by the MSBP method; for each plasma,
predictions across all 100 sub-panels are retrieved and combined.

## Method overview

The method has four main steps:

1. **Draw candidate sub-panels.** Randomly draw 63,000 smaller virus sub-panels
   (15–35 viruses out of the 41 reference viruses; 3,000 sub-panels per size).
2. **Score each sub-panel.** For each sub-panel, cluster the bnAbs by
   complete-linkage clustering on the correlation distance and compare this
   clustering to the bnAbs' known epitope classes, using *homogeneity* — an
   entropy-based measure of similarity between two clusterings.
3. **Select the best sub-panels.** Keep the 100 sub-panels with the highest
   homogeneity, i.e. those best able to classify bnAbs into their epitope classes.
4. **Predict plasma specificity.** For each plasma, predict the neutralization
   specificity on each of the 100 sub-panels using the classical Spearman
   fingerprinting method, then combine the results with two complementary
   approaches:
   - **(a) Individual bnAb-based predictions.** The distribution of plasma-vs-bnAb
     correlations across the 100 sub-panels, highlighting bnAbs with a mean
     Spearman correlation above 0.40.
   - **(b) bnAb class-based predictions.** The percentage of epitope predictions
     per bnAb class across all sub-panels. Correlations below 0.40 are grouped
     into an "unknown" class. Each percentage estimates the contribution of a
     bnAb class to the plasma's neutralizing capacity.

## Files

| File | Role |
|------|------|
| `run_PSS.R` | **Entry point.** Reads the supplementary tables, formats the data, runs the whole workflow, and writes all (date-stamped) outputs. |
| `create_virus_panels.R` | Defines `build_optimal_panels()` — steps 1–3. |
| `prediction_panels.R` | Defines `prediction_panels()` — step 4. |

The two `*_panels.R` files only **define** functions; they are sourced by
`run_PSS.R` and are not run directly.

Approximate runtime (MacBook Pro M1, 8 cores): ~5–10 min for steps 1–3 and
~1-2 min for step 4.

## Requirements

R (≥ 4.0) with:

```r
install.packages(c("readxl", "openxlsx", "dplyr", "tidyr",
                   "clevr", "foreach", "doParallel"))
```

## Input data

Two of the paper's supplementary Excel files are used as input:

| File | Sheet used | Provides |
|------|------------|----------|
| `Table S3_compiled.xlsx` | `Table S3_01` | Reference bnAb neutralization (75 bnAbs × 41 viruses, IC50) |
| `Table S1.xlsx` | `Table S1_02` | Plasma/serum neutralization (donors × 41 viruses, NT50) |

`run_PSS.R` formats these sheets automatically — you do **not** need to edit the
spreadsheets by hand.

## How to run

1. Place these four files in one folder: `run_PSS.R`, `create_virus_panels.R`,
   `prediction_panels.R`, plus `Table S3.xlsx` and `Table S1.xlsx`.
2. Set that folder as the working directory (in RStudio, see the commented
   `setwd(...)` line at the top of `run_PSS.R`).
3. Run the whole pipeline:

```r
source("run_PSS.R")
```

To change the panel-building settings or the number of cores, edit the
**CONFIGURATION** block near the top of `run_PSS.R`.

## Outputs

All outputs are written to a `results/` folder, and every file name is
automatically stamped with the execution date (`YYYYMMDD`). The formatted inputs
are saved alongside them for transparency.

| Output file | Content |
|-------------|---------|
| `data_mAb_panel.xlsx` | Formatted reference bnAb panel (input to the method) |
| `data_donor.xlsx` | Formatted donor data (input to the method) |
| `allpanels_<date>.txt` | Score of every candidate sub-panel |
| `panelopt_<date>.txt` | The 100 best-scoring sub-panels (`panel_id`, `score`) |
| `panelopt_<date>.RData` | `optimal_panels`: the virus names in each top sub-panel |
| `prediction_<date>.xlsx` | Epitope-class counts per donor (bnAb class-based prediction) |
| `averagecor_<date>.xlsx` | Mean correlation per donor and bnAb |
| `cormatrix_<date>.RData` | Full bnAb × donor × sub-panel correlation array |
| `cordata_<date>.RData` | Per-donor correlation distributions (individual bnAb-based prediction) |
| `ID_index_<date>.xlsx` | Donor ID ↔ row index lookup |
