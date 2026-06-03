  # Mortality Forecasting under Climate Risk: A Stochastic Approach with Distributed Lag Non-Linear Models
This repository supports the paper ''Mortality Forecasting under Climate Risk: A Stochastic Approach with Distributed Lag Non-Linear Models (available on <a href="https://arxiv.org/abs/2506.00561">arXiv</a>)'' written by Jiacheng Min, Han Li, Thomas Nagler, and Shuanming Li. Assessing climate-driven mortality risk has become an emerging area of research in recent decades. In this paper, we propose a novel approach to explicitly incorporate climate-driven effects into both single- and multi-population stochastic mortality models. The new model consists of two components: a stochastic mortality model, and a distributed lag non-linear model (DLNM). The stochastic component captures the non-climate long-term trend, volatility, and seasonal patterns in mortality rates. The DLNM component captures non-linear and lagged effects of climate variables on mortality, as well as the impact of heat waves and cold waves across different age groups. For model calibration, we propose a novel backfitting algorithm that allows us to disentangle the climate-driven mortality risk from the non-climate-driven stochastic mortality risk. We illustrate the effectiveness and improved short-term (1--18 month) forecasting performance of our model against four alternative models, using data from three European regions: Athens, Lisbon, and Rome. Furthermore, as an application of the proposed modeling framework, we utilize future UTCI data generated from climate models to provide total mortality forecasts into 2045 across these regions under two Representative Concentration Pathway (RCP) scenarios, taking both stochastic mortality improvement trend and climate risk into account. The projections show a noticeable decrease in winter mortality alongside a rise in summer mortality, driven by a general increase in UTCI over time. Although we expect slightly lower overall mortality in the short term under RCP8.5 compared to RCP2.6, a long-term increase in total mortality is anticipated under the RCP8.5 scenario.

This repository contains the mortality dataset, UTCI dataset and necessary code to produce tables and figures in the paper.

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

## Overview
- `environment.Rproj` initializes the project.
- `Code` folder contains all the code to reproduce the results.
- `Function` subfolder within `Code` contains the core model functions.
- `Data` folder contains the raw death count data, population data, historical UTCI data, and scenario-based UTCI data (RCP2.6 and RCP8.5).

## Section 3: Data
- We visualize the weekly historical mortality rates (2015-2019) in **Figure 1**, and UTCI data in **Figure 2**.
- Use the following code to reproduce **Figure 1** and **Figure 2**.
```shell
Rscript -e "render('Code/Main_paper/Section_3_Data/01_Historical_data_aggregation.Rmd')"
``` 

## Section 4: Empirical results
- We calibrate stochastic mortality component (**Section 4.1**) by visualizing the fitted time-varying factors for LC, LL, DLNM--LC, and DLNM--LL models in **Figure 3** and **Figure 4**. Then we calibrate DLNM climate-driven mortality component (**Section 4.2**) in DLNM--LC and DLNM--LC model by visualizing overall cumulative effects of UTCI (**Section 4.3**) in Figure 5 and present the coefficients of $\text{HWD}_t$ and $\text{CWD}_t$ in Table 1. 
- Use the following code to reproduce **Figure 3, 4**, and **5**, and **Table 1**.
```shell
Rscript -e "render('Code/Main_paper/Section_4_Empirical_results/01_DLNM_LC_LL_calibration.Rmd')"
```
- We perform the expanding-window cross-validation in **Section 4.4** for 6 models: **LC, LL, DLNM--LC, DLNM--LL, Madaniyazi et al., Guibert et al.**. The forecast mean absolute error (MAE) under $\times 100$ scale is reported in **Table 2**.
- Use the following code to reproduce **Table 2**.
```shell
Rscript -e "render('Code/Main_paper/Section_4_Empirical_results/02_model_comparison_MAE_table.Rmd')"
```

## Section 5: Mortality projections under RCP scenarios
- We first process the future UTCI data under RCP2.6 and RCP8.5.
- Use the following code to process the future UTCI data.
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/00_RCP_future_input_processing.Rmd')"
```
- We present the results of weekly mortality projections (**Section 5.3**) under RCP2.6 and RCP8.5 in **Figure 7, 8**, and **9**. 
- Use the following code to run the weekly mortality projections for DLNM--LC and DLNM--LL model, respectively.    
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/01_DLNM_LC_RCP_simulation.Rmd')"
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/02_DLNM_LL_RCP_simulation.Rmd')"
```
- Use the following code to reproduce **Figure 7, 8**, and **9**.
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/03_RCP_weekly_visualization.Rmd')"
```
- We illustrate the results of annual mortality projections (**Section 5.4**) in **Figure 10** of the paper, and **Figure F.1** and **F.2** in the Supplementary Material.
- Use the following code to run the annual mortality projections for DLNM--LC and DLNM--LL model, respectively.
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/04_RCP_annualization.Rmd')"
```
- Use the following code to reproduce **Figure 10**, **Figure F.1** and **F.2**.
```shell
Rscript -e "render('Code/Main_paper/Section_5_Mortality_projection/05_RCP_annualized_visualization.Rmd')"
```
## Supplementary Material
- In **Section C**, we present the lagged effects via DLNM across all regions for DLNM--LC and DLNM--LL models in **Figure C.1, C.2, C.3, C.4, C.5** and **C.6**.
- Use the following code to reproduce **Figure C.1, C.2, C.3, C.4, C.5** and **C.6**.
```shell
Rscript -e "render('Code/Supplementary/Section_C_Figures_of_specific_lagged_effects_via_DLNM.Rmd')"
```
- In Section D, we present the optimal time series model for time-varying factors (**Section D.1**) and seasonality test on model residuals (**Section D.2**).
- Use the following code to reproduce **Table D.1, D.2, D.3** and **D.4**.
```shell
Rscript -e "render('Code/Supplementary/Section_D_Model_specifications_and_diagnostics.Rmd')"
```
- In **Section E**, we report the forecast performance under alternative heat and cold wave definitions, reporting forecast MAE from expanding-window cross-validation in **Table E.5** and **E.6**.
- Use the following code to reproduce **Table E.5, E.6**.
```shell
Rscript -e "render('Code/Supplementary/Section_E_Alternative_heat_and_cold_wave_definitions.Rmd')"
```

## References
Min, J., Li, H., Nagler, T., & Li, S. (2025). *Assessing Climate-Driven Mortality Risk: A Stochastic Approach with Distributed Lag Non-Linear Models*. arXiv preprint https://arxiv.org/abs/2506.00561.
