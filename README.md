# Mortality Forecasting under Climate Risk: A Stochastic Approach with Distributed Lag Non-Linear Models
This repository supports the paper ''Mortality Forecasting under Climate Risk: A Stochastic Approach with Distributed Lag Non-Linear Models (available on <a href="https://arxiv.org/abs/2506.00561">arXiv</a>)'' written by Jiacheng Min, Han Li, Thomas Nagler, and Shuanming Li. In this paper, we propose a novel approach to explicitly incorporate climate-driven effects into both single- and multi-population stochastic mortality models. The new model con-
sists of two components: a stochastic mortality model, and a distributed lag non-linear model (DLNM). The stochastic component captures the non-climate long-term trend, volatility, and seasonal patterns in mortality rates. The DLNM component captures non-linear and lagged effects of climate variables on mortality, as well as the impact of heat waves and cold waves across different age groups. For model calibration, we propose a novel backfitting algorithm that allows us to disentangle the climate-driven mortality risk from the non-climate-driven stochastic mortality risk. We illustrate the effectiveness and superior forecasting performance of our model against four alternative models, using data from three European regions: Athens, Lisbon, and Rome. Furthermore, we utilize future UTCI data generated from climate models to provide mortality forecasts into 2045 across these regions under two Representative Concentration Pathway (RCP) scenarios, taking both stochastic mortality improvement trend and climate risk into account. The projections show a noticeable decrease in winter mortality alongside a rise in summer mortality, driven by a general increase in UTCI over time. Although we expect slightly lower overall mortality in the short term under RCP 8.5 compared to RCP 2.6, a long-term increase in total mortality is anticipated under the RCP8.5 scenario.

This repository contains the mortality dataset, UTCI dataset and necessary code to produce figures in the paper.

## Package Requirements
The following packages are required to finish the experiments.
- `dlnm`
- `abcd`

All packages are available on CRAN.

## Overview
- `Data` folder contains the raw death count data, population data, historical UTCI data, and simulated UTCI data (RCP 2.6 and RCP 8.5).
- `Function` folder contains the necessary functions for LC, LL, DLNM--LC, DLNM--LL, and Guibert et al's model. 
- `B-DLNM-LC` and `B-DLNM-LL` folders contain DLNM--LC and DLNM--LL model validation, including expanding cross-validation results and MAE computations.
- `Calibration` contains DLNM--LC and DLNM--LL model calibration, associating with Section 4 in the paper.
