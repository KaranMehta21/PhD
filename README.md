# Development and Validation of All-Cause Mortality Risk Prediction Models for Individuals with Type 2 Diabetes

PhD thesis code — Karan Mehta, UCL Institute of Health Informatics (Health Data Science).

This repository contains the analysis code for the corrected thesis: the cohort comparison with the National Diabetes Audit (Chapter 5), the mortality analysis (Chapter 6) and the development and validation of the prediction models (Chapter 7).

## Data

The analyses use **IQVIA Medical Research Data (IMRD)**. IQVIA Medical Research Data (IMRD) incorporates data from THIN, A Cegedim Database. Reference made to THIN is intended to be descriptive of the data asset licensed by IQVIA. Use of the data was approved by the IQVIA Scientific Review Committee (reference 20SRC056).

**No data are included in this repository.** IMRD data are available from IQVIA under licence and cannot be shared. The scripts read data from the UCL secure research drive; file paths are set as constants at the top of each script and must be changed to run elsewhere. Scripts write only aggregate outputs (counts, rates, model estimates, figures).

## Repository structure

```
chapter5/
  c51_nda_benchmark_aggregates.py         descriptive aggregates for the comparison with the National Diabetes Audit
chapter6/
  01_build_2009_matched_cohort.py         T2D and frequency-matched comparison cohort at 1 Jan 2009, follow-up, covariates, summary table
  02_mortality_rate_ratio_models.py       Poisson models (A-D), age-sex specific rates and MRRs, interaction tests, E-values, stepwise comparison
  03_poisson_assumption_check.py          overdispersion statistics and negative binomial sensitivity analysis
  04_figures_5_to_8.py                    Figures 5-8 from the aggregate outputs of 02
chapter7/
  original_2023_2024/                     original analysis code (Stata do-files; notebooks exported as code only)
    00a_kaplan_meier_logrank.py           Kaplan-Meier curves and log-rank tests
    00b_univariable_cox_correlation.py    univariable Cox regression and correlation matrix
    01_mvrs_nonlinear_terms.do            multivariable regression splines (mvrs, select(0.05))
    02_interaction_terms.do               interactions with sex and age (likelihood ratio tests)
    03_compare_cox_models_AIC_BIC_LR.do   comparison of Cox models with/without splines and interactions
    04a_mice_with_splines_interactions_150224.do   multiple imputation including spline and interaction terms
    04b_mice_final_180224.do              final multiple imputation do-file (mi impute chained, 20 imputations)
  01_predictor_selection_funnel.py        reconstruction of the predictor selection steps (log-rank, univariable Cox, correlation)
  02_followup_exits_table.py              reasons for leaving follow-up and person-time
  03_refit_four_models_evaluation.py      four models on one seeded 70/30 split; C-index and IBS with bootstrap CIs; paired
                                          comparisons; calibration; proportional hazards tests; final model on all data
  04_calibration_figures_10_to_13.py      calibration figures from the saved test-set predictions
  05_survivaleval_crosscheck.py           D-calibration with SurvivalEVAL
  06_prediction_diagnostics.py            diagnostics of predicted survival distributions (quality control)
  07_predict_individual.py                applies the saved models to one individual's raw predictor values
```

## Thesis tables and figures

| Thesis item | Script |
|---|---|
| Table 7 (comparison with the National Diabetes Audit) | `chapter5/c51_nda_benchmark_aggregates.py` |
| Table 8 (study population, deaths, person-years) | `chapter6/01_build_2009_matched_cohort.py` |
| Table 9 (mortality rates and MRRs by age and sex) | `chapter6/02_mortality_rate_ratio_models.py` |
| Table 10 (extended adjustment and sensitivity analyses) | `chapter6/02_mortality_rate_ratio_models.py` |
| Overdispersion checks (Section 6.4) | `chapter6/03_poisson_assumption_check.py` |
| Figures 5–8 | `chapter6/04_figures_5_to_8.py` |
| Table 12 (reasons for leaving follow-up) | `chapter7/02_followup_exits_table.py` |
| Table 13 (univariable Cox regression) | `chapter7/original_2023_2024/00b_univariable_cox_correlation.py` (reproduced by `chapter7/01_predictor_selection_funnel.py`) |
| Predictor selection steps (Section 7.4) | `chapter7/01_predictor_selection_funnel.py` |
| Table 14 (C-index and IBS with 95% CIs) | `chapter7/03_refit_four_models_evaluation.py` |
| Table 15 (paired differences between models) | `chapter7/03_refit_four_models_evaluation.py` |
| Table 16 (proportional hazards tests) | `chapter7/03_refit_four_models_evaluation.py` |
| Figures 10–13 (calibration) | `chapter7/04_calibration_figures_10_to_13.py` |
| D-calibration (Section 7.4) | `chapter7/05_survivaleval_crosscheck.py` |
| Supplementary Methods 7A1–7A2 (Kaplan–Meier plots, log-rank p-values) | `chapter7/original_2023_2024/00a_kaplan_meier_logrank.py` |
| Supplementary 7A3 (final Cox model: coefficients, baseline survival) | `chapter7/03_refit_four_models_evaluation.py` |
| Non-linear terms, interactions, imputation (Sections 7.3.5.1–7.3.5.3) | `chapter7/original_2023_2024/01_*.do` to `04b_*.do` |

Table 11 (baseline characteristics) reproduces the Chapter 5 cohort description. Figure 9 is a diagram.

Note: the final imputed dataset used for modelling (`impv2_phddata_230224`) was produced by the multiple imputation do-files in `original_2023_2024/`; the exact version of the do-file that generated that file was not located, and `04b_mice_final_180224.do` is the closest version.

## Software

- Stata — `mvrs`, `mi impute chained`
- Python 3.11.7 — pandas 2.1.4, NumPy 1.26.4, statsmodels, lifelines 0.29.0, scikit-survival 0.23.0, pycox 0.2.3, PyTorch 2.4.0, torchtuples 0.2.2, scikit-learn 1.5.1, SurvivalEVAL, SciPy, matplotlib

## Licence

Code released under the MIT licence (see `LICENSE`). The data are not covered by this licence.
