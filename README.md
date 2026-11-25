# Longitudinal HDL-C and Interstitial Lung Abnormalities

R analysis specifications for **Pre-Screening Longitudinal HDL Cholesterol Measures and Interstitial Lung Abnormalities in Korean National Lung Cancer Screening**.

This repository documents how longitudinal high-density lipoprotein cholesterol (HDL-C) exposures are constructed and how their associations with prevalent interstitial lung abnormalities (ILAs) are modeled. It provides explicit covariate coding, functional forms, spline specifications, and confidence interval calculations for marginally standardized prevalences and prevalence differences.

## Analysis code

The main script is [`hdl_ila_analysis_specification.R`](hdl_ila_analysis_specification.R).

It includes:

- Construction of time-weighted mean HDL-C, the proportion of low HDL-C measurements, and the lowest HDL-C value.
- Final multivariable logistic regression Models 1 and 2.
- Marginal standardization and delta-method 95% confidence intervals.
- Sensitivity analyses using continuous recent HDL-C, CT-reported emphysema, inverse probability of selection weighting, and time-weighted mean HDL-C excluding the most recent measurement.
- Joint categories of time-weighted mean and recent HDL-C status.
- Functions for the final latent class mixed models under first-measurement and screening-date anchoring, followed by hard-assignment and posterior probability–weighted logistic regression.
- Pooled and participant-level summaries of intermeasurement intervals.

## Key specifications

**Exposure construction.** Valid HDL-C measurements are restricted to 10–150 mg/dL from January 1, 2009, through the baseline lung cancer screening (LCS) date. At least three valid measurements and a positive observation span are required. Low HDL-C is defined as <40 mg/dL in men and <50 mg/dL in women.

Time-weighted mean HDL-C is the trapezoidal area under the HDL-C–time curve divided by the interval between the first and last measurements. The interval from the last measurement to LCS is excluded. The proportion of low measurements is visit-based and expressed as a percentage. For tied lowest values, the most recent measurement determines the interval to LCS. The historical sensitivity analysis requires at least three measurements remaining after excluding the latest measurement.

**Model definitions.** Model 1 adjusts for age, sex, body mass index category, smoking status, pack-years, drinking status, family history of lung cancer, history of other cancer, COPD, diabetes, hypertension, ischemic heart disease, stroke, region, low-income status, statin use, LDL-C and triglyceride categories, time from the recent HDL-C measurement to LCS, number of measurements, and observation duration. Model 2 additionally includes recent HDL-C status. The lowest HDL-C analysis uses Model 1 only and substitutes time from the lowest measurement to LCS. Trajectory association models omit measurement count and observation duration.

Continuous adjustment covariates enter linearly unless otherwise specified; categorical covariates use treatment coding. Category boundaries and reference levels are defined in the script. Odds ratios for the proportion of low measurements are expressed per 10-percentage-point increase; those for lowest HDL-C are expressed per 5 mg/dL decrease.

**Splines.** The recent HDL-C sensitivity analysis uses four knots at 33, 44, 54, and 75 mg/dL, corresponding to the reported 5th, 35th, 65th, and 95th percentiles. The natural cubic spline basis has linear tails beyond the outer knots. The final trajectory models use a separate three-degree-of-freedom natural cubic spline of time, with internal knots at the one-third and two-thirds quantiles of observed measurement times and boundaries at their minimum and maximum.

**Standardization and confidence intervals.** Predicted probabilities are averaged after assigning each exposure value to the standardization population while retaining observed covariates. Confidence intervals use the delta method, retaining the covariance between standardized estimates when calculating prevalence differences. Unweighted models use model-based covariance matrices. Selection-weighted models use HC0 robust covariance, and posterior probability–weighted models use participant-clustered robust covariance. These calculations condition on the empirical covariate distribution and treat estimated selection weights and posterior probabilities as fixed. Prevalences are expressed as percentages and differences as percentage points.

## Requirements

The study-reported software environment was **R 4.0.5** and **lcmm 2.2.2**.

| Component | R packages |
|---|---|
| Exposure construction, logistic regression, and standardization | Base R and the bundled `splines` package |
| Selection-weighted and posterior probability–weighted regression | `sandwich` |
| Optional latent class mixed model fitting | `lcmm` |

The script does not install packages automatically. Each analysis run records its executing environment in `sessionInfo.txt`.

## Input data

Participant-level study data are not included. Analyses require separately obtained, appropriately authorized inputs, with source variables mapped to the names and categories documented at the beginning of the R script.

| Input object | Required structure |
|---|---|
| `lcs` | One row per participant, with the outcome, screening date, covariates, and preselected recent lipid panel |
| `hdl` | Longitudinal measurements with columns `id`, `date`, and `hdl` |
| `elig` — optional | Participants reaching the ILA-assessment step, including assessment status and selection-model covariates |
| `posterior_fu` — optional | First-measurement–anchored posterior table: `id`, `class`, `prob1`, `prob2`, `prob3` |
| `posterior_lcs` — optional | LCS-anchored posterior table with the same columns |

For the CT-emphysema sensitivity analysis, supply `emphysema` in `lcs` as 0/1, with missing assessments coded `NA`. Posterior class labels must correspond to the manuscript's Class 1, 2, and 3.

Study eligibility exclusions based on claims, questionnaire processing, and recent lipid-panel selection are performed upstream. The current recent-panel date validator permits dates from January 1 of the screening year minus four through the screening date. The selected panel must contain HDL-C, LDL-C, and triglycerides. The script documents this calendar-year rule explicitly.

## Running the analyses

After preparing the input data frames, run the following in R from the repository directory:

```r
source("hdl_ila_analysis_specification.R")

input <- list(lcs = lcs, hdl = hdl)

# Add available optional inputs:
# input$elig <- elig
# input$posterior_fu <- posterior_fu
# input$posterior_lcs <- posterior_lcs

results <- run_hdl_ila(input, output_dir = "results")
```

Use a new or empty output directory. The runner writes aggregate CSV results and analysis specifications; fitted models are returned in the R object `results`.

Alternatively, save the input list in R:

```r
saveRDS(input, "analysis_inputs.rds")
```

Then run from a terminal:

```bash
Rscript hdl_ila_analysis_specification.R analysis_inputs.rds results
```

The output includes adjusted odds ratios, standardized prevalence contrasts, spline specifications and likelihood-ratio tests, model formulas, sample counts, and applicable weighting diagnostics. Optional analyses run when their corresponding inputs are supplied.

The common standardization function can also be applied to individual returned models. For example, after fitting the emphysema sensitivity analysis:

```r
standardized_prevalence(
  results$emphysema$Emphysema_Model2,
  exposure = "low_twm",
  values = c("Normal", "Low")
)
```

CT-emphysema standardized prevalences are available through this function but are not automatically exported by the current runner.

## Trajectory fitting and scope

`fit_final_trajectory()` implements the selected three-class, three-degree-of-freedom trajectory specification. It accepts a saved converged subsample model for initialization; otherwise, it fits that selected specification in a seeded 50% subsample before the full-cohort fit. The main runner does not automatically launch LCMM fitting. Usage examples and explicit class-label mapping are provided in the R script.

This repository focuses on exposure construction and final analysis specifications. The full exploratory comparison of class counts and spline degrees of freedom, subgroup and interaction analyses, unadjusted analyses, and manuscript table/figure production are outside the automated workflow. Numerical reproduction additionally requires the original prepared data, trajectory initialization, and software environment.

## Implementation checks

A self-test using synthetic data is included:

```bash
Rscript hdl_ila_analysis_specification.R --self-test
```

It checks irregular-interval TWM calculation, eligibility after removing the latest measurement, unit rescaling, complete-case alignment, delta-method gradients, and spline comparison degrees of freedom. It does not use study records.

This version has undergone a static methods review. R execution, including the synthetic self-test, and numerical reproduction of study estimates have not yet been verified for this version.

## Citation

If you use or adapt this code, please cite the associated study, **Pre-Screening Longitudinal HDL Cholesterol Measures and Interstitial Lung Abnormalities in Korean National Lung Cancer Screening**, and identify the repository version or commit used. Publication details will be added when available.

## License

The code is distributed under the [MIT License](LICENSE). This license applies to the code and does not grant access to participant-level study data.
