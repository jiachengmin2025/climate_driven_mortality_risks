## Function

`Code/Function/` contains the shared R functions used by the main paper and
supplementary R Markdown files.

### Stochastic Mortality Models

- `LC_model.R` fits the Lee--Carter model using singular value decomposition.
- `LL_model.R` fits the Li--Lee model by combining a common mortality trend with
  region-specific deviations.

### DLNM--LC and DLNM--LL Models

- `DLNM_LC.R` fits the proposed DLNM--LC model using the backfitting algorithm.
- `DLNM_LC.forecast.R` forecasts mortality rates from a fitted DLNM--LC model.
- `DLNM_LL.R` fits the proposed DLNM--LL model using the backfitting algorithm.
- `DLNM_LL.forecast.R` forecasts mortality rates from a fitted DLNM--LL model.
- `dlnm_proc.R` fits the age-specific Gaussian DLNM component on log mortality
  rates.

### Cross-Basis Matrices and Cross-Validation

- `crossbasis_mixed_frequency.R` builds historical weekly cross-basis matrices
  from daily UTCI data.
- `create_train_test_sets.R` constructs expanding-window training and testing
  sets for DLNM--LC and DLNM--LL model comparison.

### Comparator Models

- `Madaniyazi.fit.R` and `Madaniyazi.forecast.R` fit and forecast the
  DLNM--Madaniyazi comparator.
- `Guibert.fit.R`, `Guibert.forecast.R`, and
  `Guibert.create_train_test_sets.R` fit, forecast, and construct
  expanding-window inputs for the Guibert et al. comparator.

### Time-Series Models and Diagnostics

- `kt_model.R` fits time-series models for time-varying factors.
- `kappa_fit.R` forecasts common and region-specific time-varying
  factors.
- `sarima_table_helpers.R` formats optimal SARIMA models and coefficients for
  Supplementary Section D.
- `DLNM_residual_diagnostics.R` runs the residual seasonality diagnostics for
  Supplementary Section D.

### RCP Mortality Projection

- `rcp_inputs.R` prepares historical mortality inputs and loads future RCP
  cross-basis and wave-indicator inputs.
- `rcp_dlnm_simulation.R` simulates future mortality under RCP scenarios using
  fitted DLNM--LC and DLNM--LL models.
- `rcp_annualization.R` converts weekly simulated mortality rates into
  annualized mortality rates.
- `rcp_visualization.R` prepares weekly and annualized RCP simulation results for
  plotting.
- `rcp_legacy_sarima.R` stores the SARIMA specifications used in the mortality
  projection workflow.
- `rcp_region_names.R` standardizes region labels in RCP simulation outputs.
