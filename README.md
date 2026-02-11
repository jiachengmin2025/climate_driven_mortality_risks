# Mortality Forecasting under Climate Risk: A Stochastic Approach with Distributed Lag Non-Linear Models
This repository supports the paper ''Mortality Forecasting under Climate Risk: A Stochastic Approach with Distributed Lag Non-Linear Models (available on <a href="https://arxiv.org/abs/2506.00561">arXiv</a>)'' written by Jiacheng Min, Han Li, Thomas Nagler, and Shuanming Li. In this paper, we propose a novel approach to explicitly incorporate climate-driven effects into both single- and multi-population stochastic mortality models. The new model con-
sists of two components: a stochastic mortality model, and a distributed lag non-linear model (DLNM). The stochastic component captures the non-climate long-term trend, volatility, and seasonal patterns in mortality rates. The DLNM component captures non-linear and lagged effects of climate variables on mortality, as well as the impact of heat waves and cold waves across different age groups. For model calibration, we propose a novel backfitting algorithm that allows us to disentangle the climate-driven mortality risk from the non-climate-driven stochastic mortality risk. We illustrate the effectiveness and superior forecasting performance of our model against four alternative models, using data from three European regions: Athens, Lisbon, and Rome. Furthermore, we utilize future UTCI data generated from climate models to provide mortality forecasts into 2045 across these regions under two Representative Concentration Pathway (RCP) scenarios, taking both stochastic mortality improvement trend and climate risk into account. The projections show a noticeable decrease in winter mortality alongside a rise in summer mortality, driven by a general increase in UTCI over time. Although we expect slightly lower overall mortality in the short term under RCP 8.5 compared to RCP 2.6, a long-term increase in total mortality is anticipated under the RCP8.5 scenario.

This repository contains the mortality dataset, UTCI dataset and necessary code to produce tables and figures in the paper.

## Package Requirements
The following packages are required to finish the experiments.
- `readxl`, `writexl`, `dplyr`, `tidyr`, `reshape2`, `zoo`, `ISOweek` for data imputation.
- `dlnm`, `splines`, `forecast`, `astsa`, `MTS`, `mgcv`, `demography` for modelling.
- `seastests`, `tseries`, `nortest`, `hwwntest` for model diagnostics and statistical tests.
- `ggplot2`, `patchwork`, `gridExtra`, `ggpubr`, `scales` for visualizations.

All packages are available on CRAN.

## Overview
- `LC_DLNM.Rproj` initializes the project.
- `Code` folder contains all the code to reproduce the results.
  - The Function subfolder within Code contains the core model functions.
- `Data` folder contains the raw death count data, population data, historical UTCI data, and scenario-based UTCI data (RCP2.6 and RCP8.5).
- `Figures` folder includes all the produced figures.
- `Results` folder is not shown on Github due to the large size. The folder is available on <a href="https://drive.google.com/drive/folders/1MTxlj_AZO0nfjLoOq5icurK3x4g1YJEW?usp=sharing">Google Drive</a>. The folder contains necessary `.RData` files for mortality projections under RCP scenarios. These `.RData` files can also be reproduced via R scripts in `Code/RCP_simulation_study`.

## Empirical results
- Use the following code to visualize mortality rates and UTCI across regions and age groups.

```shell
Rscript Code/Historical_data_aggregation/Historical_data.R
```

- Use the following code to generate mixed-frequency cross-basis matrices across regions.

```shell
Rscript Code/Crossbasis_matrix/cb_Attiki.R
Rscript Code/Crossbasis_matrix/cb_Lisbon.R
Rscript Code/Crossbasis_matrix/cb_Roma.R
```

- Use the following code to run DLNM--LC and DLNM--LL model. In the expanding window cross-validation, the mean absolute error (MAE) for LC, LL, DLNM--LC, and DLNM--LL can be computed and compared.

```shell
Rscript Code/Model/DLNM-LC-Attica.R
Rscript Code/Model/DLNM-LC-Lisbon.R
Rscript Code/Model/DLNM-LC-Rome.R
Rscript Code/Model/DLNM-LL.R
```
- Use the following code to calibrate DLNM--LC and DLNM--LL model. The plot of overall cumulative effects, lagged effects, and time-varying factors are generated.

```shell
Rscript Code/Calibration/DLNM-LC-Calibration.R
Rscript Code/Calibration/DLNM-LL-Calibration.R
Rscript Code/Calibration/time-varying-factors.R
```

- Use the following code to run DLNM. In the expanding window cross-validation, the mean absolute error (MAE) can be computed.

```shell
Rscript Code/Model-comparison/time-spline-DLNM-Attiki.R
Rscript Code/Model-comparison/time-spline-DLNM-Lisbon.R
Rscript Code/Model-comparison/time-spline-DLNM-Roma.R
```
## Mortlaity projections under RCP scenarios
- Use the following code to prepare future generated UTCI data under RCP2.6 and RCP8.5.

```shell
Code/RCP_simulation_study/data_processing_rcp26_rcp85.R
```

- Use the following code to simulate weekly mortality rates from 2031 to 2045 under RCP2.6 and RCP8.5 for DLNM--LC model, and visualize the results.

```shell
Rscript Code/RCP_simulation_study/DLNM_LC_sim_rcp26.R
Rscript Code/RCP_simulation_study/DLNM_LC_sim_rcp85.R
Rscript Code/RCP_simulation_study/DLNM_LC_RCP_sim_comparison.R
```

- Use the following code to simulate weekly mortality rates from 2031 to 2045 under RCP2.6 and RCP8.5 for DLNM--LL model, and visualize the results.

```shell
Rscript Code/RCP_simulation_study/DLNM_LL_sim_rcp26.R
Rscript Code/RCP_simulation_study/DLNM_LL_sim_rcp85.R
Rscript Code/RCP_simulation_study/DLNM_LL_RCP_sim_comparison.R
```
- The simulation results are saved as `.RData` file in `Results`. If your computer cannot save the simulation outputs due to limited memory or computational resources, please download the `Results` folder on <a href="https://drive.google.com/drive/folders/1MTxlj_AZO0nfjLoOq5icurK3x4g1YJEW?usp=sharing">Google Drive</a> to conitinue reproducing the simulations.

- Use the following code to aggregate annualized mortality data, simulate annualized mortality rates from 2031 to 2045 under RCP2.6 and RCP8.5, and visualize the results.

```shell
Rscript Code/RCP_simulation_study/Annualization_data_processing.R
Rscript Code/RCP_simulation_study/Annualized_DLNM_LC_plotting.R
Rscript Code/RCP_simulation_study/Annualized_DLNM_LL_plotting.R
```

## References
Min, J., Li, H., Nagler, T., & Li, S. (2025). *Assessing Climate-Driven Mortality Risk: A Stochastic Approach with Distributed Lag Non-Linear Models*. arXiv preprint https://arxiv.org/abs/2506.00561.
