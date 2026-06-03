# Mortality Forecasting under Climate Risk: A Stochastic Approach with Distributed Lag Non-Linear Models
This repository supports the paper ''Mortality Forecasting under Climate Risk: A Stochastic Approach with Distributed Lag Non-Linear Models (available on <a href="https://arxiv.org/abs/2506.00561">arXiv</a>)'' written by Jiacheng Min, Han Li, Thomas Nagler, and Shuanming Li. Assessing climate-driven mortality risk has become an emerging area of research in recent decades. In this paper, we propose a novel approach to explicitly incorporate climate-driven effects into both single- and multi-population stochastic mortality models. The new model consists of two components: a stochastic mortality model, and a distributed lag non-linear model (DLNM). The stochastic component captures the non-climate long-term trend, volatility, and seasonal patterns in mortality rates. The DLNM component captures non-linear and lagged effects of climate variables on mortality, as well as the impact of heat waves and cold waves across different age groups. For model calibration, we propose a novel backfitting algorithm that allows us to disentangle the climate-driven mortality risk from the non-climate-driven stochastic mortality risk. We illustrate the effectiveness and improved short-term (1--18 month) forecasting performance of our model against four alternative models, using data from three European regions: Athens, Lisbon, and Rome. Furthermore, as an application of the proposed modeling framework, we utilize future UTCI data generated from climate models to provide total mortality forecasts into 2045 across these regions under two Representative Concentration Pathway (RCP) scenarios, taking both stochastic mortality improvement trend and climate risk into account. The projections show a noticeable decrease in winter mortality alongside a rise in summer mortality, driven by a general increase in UTCI over time. Although we expect slightly lower overall mortality in the short term under RCP8.5 compared to RCP2.6, a long-term increase in total mortality is anticipated under the RCP8.5 scenario.

This repository contains the mortality dataset, UTCI dataset and necessary code to reproduce tables and figures in the paper.

## Overview
- `environment.Rproj` initializes the project.
- `Code` folder contains all code to reproduce the results. **To facilitate reproducibility checks, we also provide knitted PDF files in the** `Code` **folder showing that the code runs and reproduces the results.**
- `Data` folder contains the following datasets:
  - The death count data and population data from Eurostat can be found in `Data/Combined_data`,
  - Historical UTCI data from ERA5 can be found in `Data/UTCI_data/Daily_data`,
  - Scenario-based UTCI data can be found in `Data/Simulation_data/UTCI`.

## Package Requirements
The following packages are required to finish the experiments.
- `readxl`, `writexl`, `dplyr`, `tidyr`, `reshape2`, `zoo`, `ISOweek` for data imputation.
- `dlnm`, `splines`, `forecast`, `astsa`, `MTS`, `mgcv`, `demography` for modelling.
- `seastests`, `tseries`, `nortest`, `hwwntest`, `Metrics`, for model diagnostics and statistical tests.
- `ggplot2`, `patchwork`, `gridExtra`, `ggpubr`, `scales` for visualizations.

All packages are available on CRAN. To check whether all packages are installed or find out which are missing run:
```shell
Rscript package_requirements.R
```

## Function
`Code/Function/` contains main functions used in the project.

- `LC_model.R` and `LL_model.R` fit the Lee--Carter and Li--Lee models.
- `DLNM_LC.R`, `DLNM_LC.forecast.R`, `DLNM_LL.R`, and `DLNM_LL.forecast.R` fit and forecast DLNM--LC and DLNM--LL models.
- `dlnm_proc.R` fits the age-specific Gaussian DLNM component on log mortality rates.
- `create_train_test_sets.R` constructs expanding-window training and testing sets.
- `crossbasis_mixed_frequency.R` builds historical weekly cross-basis matrices from daily UTCI data.
- `Madaniyazi.fit.R`, `Madaniyazi.forecast.R`, `Guibert.fit.R`, `Guibert.forecast.R`, and `Guibert.create_train_test_sets.R` implement the Madaniyazi et al. and Guibert et al. models.
- `kt_model.R`, `kappa_fit.R`, and `sarima_table_helpers.R` fit and summarize time-series models for the time-varying mortality factors.
- `DLNM_residual_diagnostics.R` runs the residual seasonality diagnostics in Supplementary Section D.
- `rcp_inputs.R`, `rcp_dlnm_simulation.R`, `rcp_annualization.R`, `rcp_visualization.R`, `rcp_legacy_sarima.R`, and `rcp_region_names.R` support the RCP mortality projection and visualization workflow.

## Manuscript
The code related to  **Section 3, 4,** and **5** of the manuscript can be found in `Code/Main_paper`.
### Section 3: Data
- The code for **Section 3** of the manuscript can be found `Code/Main_paper/Section_3_Data`.
- We visualize the weekly historical mortality rates (2015-2019) in **Figure 1**, and UTCI data in **Figure 2**.
- Use the following code to reproduce **Figure 1** and **Figure 2** of the manuscript.
```shell
Rscript -e "render('Code/Main_paper/Section_3_Data/01_Historical_data_aggregation.Rmd')"
``` 

### Section 4: Empirical results
- The code for **Section 4** of the manuscript can be found `Code/Main_paper/Section_4_Empirical_results`.
- We calibrate stochastic mortality component (**Section 4.1**) by visualizing the fitted time-varying factors for LC, LL, DLNM--LC, and DLNM--LL models in **Figure 3** and **Figure 4**. Then we calibrate DLNM climate-driven mortality component (**Section 4.2**) in DLNM--LC and DLNM--LC model by visualizing overall cumulative effects of UTCI (**Section 4.3**) in **Figure 5** and present the coefficients of $\text{HWD}_t$ and $\text{CWD}_t$ in **Table 1**. 
- Use the following code to reproduce **Figure 3, 4**, and **5**, and **Table 1** of the manuscript.
```shell
Rscript -e "render('Code/Main_paper/Section_4_Empirical_results/01_DLNM_LC_LL_calibration.Rmd')"
```
- We perform the expanding-window cross-validation in **Section 4.4** for six models: **LC, LL, DLNM--LC, DLNM--LL, Madaniyazi et al., Guibert et al.**. The forecast mean absolute error (MAE) under $\times 100$ scale is reported in **Table 2**.
- Use the following code to reproduce **Table 2** of the manuscript.
```shell
Rscript -e "render('Code/Main_paper/Section_4_Empirical_results/02_model_comparison_MAE_table.Rmd')"
```

### Section 5: Mortality projections under RCP scenarios
- The code for **Section 5** of the manuscript can be found `Code/Main_paper/Section_5_Mortality_projection`
- We first process the future UTCI data under RCP2.6 and RCP8.5.
- Use the following code to initialize and process the future UTCI data.
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/00_RCP_future_input_processing.Rmd')"
```
- We present the results of weekly mortality projections (**Section 5.3**) under RCP2.6 and RCP8.5 in **Figure 7, 8**, and **9**. 
- Use the following code to run the weekly mortality projections for DLNM--LC and DLNM--LL model, respectively.    
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/01_DLNM_LC_RCP_simulation.Rmd')"
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/02_DLNM_LL_RCP_simulation.Rmd')"
```
- Use the following code to reproduce **Figure 7, 8**, and **9** of the manuscript.
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/03_RCP_weekly_visualization.Rmd')"
```
- We illustrate the results of annual mortality projections (**Section 5.4**) in **Figure 10** of the manuscript, and **Figure F.1** and **F.2** in the Supplementary Material.
- Use the following code to run the annual mortality projections for DLNM--LC and DLNM--LL model, respectively.
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/04_RCP_annualization.Rmd')"
```
- Use the following code to reproduce **Figure 10**, **Figure F.1** of the manuscript, and **F.2** of Supplementary Material.
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/05_RCP_annualized_visualization.Rmd')"
```
## Supplementary Material
- The code related to **Section C, D** and **E** of Supplementary Material can be found in `Code/Supplementary`.
### **Section C**
- The code for **Section C** of Supplementary Material can be found `Code/Supplementary/Section_C_Figures_of_specific_lagged_effects_via_DLNM.Rmd`.
- We present the lagged effects via DLNM for DLNM--LC and DLNM--LL models in **Figure C.1, C.2, C.3, C.4, C.5** and **C.6**.
- Use the following code to reproduce **Figure C.1, C.2, C.3, C.4, C.5** and **C.6** of Supplementary Material.
```shell
Rscript -e "render('Code/Supplementary/Section_C_Figures_of_specific_lagged_effects_via_DLNM.Rmd')"
```
### **Section D**
- The code for **Section D** of Supplementary Material can be found `Code/Supplementary/Section_D_Model_specifications_and_diagnostics.Rmd`
- We present the optimal time series model for time-varying factors (**Section D.1**) and seasonality test on model residuals (**Section D.2**).
- Use the following code to reproduce **Table D.1, D.2, D.3** and **D.4** of Supplementary Material.
```shell
Rscript -e "render('Code/Supplementary/Section_D_Model_specifications_and_diagnostics.Rmd')"
```
### **Section E**
- The code for **Section E** of Supplementary Material can be found `Code/Supplementary/Section_E_Alternative_heat_and_cold_wave_definitions.Rmd`
- We report the forecast performance under alternative heat and cold wave definitions, reporting forecast MAE from expanding-window cross-validation in **Table E.5** and **E.6**.
- Use the following code to reproduce **Table E.5, E.6** of Supplementary Material.
```shell
Rscript -e "render('Code/Supplementary/Section_E_Alternative_heat_and_cold_wave_definitions.Rmd')"
```

## References
Min, J., Li, H., Nagler, T., & Li, S. (2025). *Assessing Climate-Driven Mortality Risk: A Stochastic Approach with Distributed Lag Non-Linear Models*. arXiv preprint https://arxiv.org/abs/2506.00561.
