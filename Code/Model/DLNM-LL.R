## Load packages
packages <- c("readxl", "dlnm", "splines", "mgcv", "forecast", "demography", 
              "Metrics", "tseries", "nortest", "seastests")
invisible(lapply(packages, library, character.only = TRUE))

## Load data
region_name = c("Attiki", "Lisbon", "Roma")
age = c('Y20_64', 'Y65_74', 'Y75_84', 'Y_GE85')
Y20_64 = read_xlsx("Data/Combined_data/Y20_64_combined.xlsx")
Y65_74 = read_xlsx("Data/Combined_data/Y65_74_combined.xlsx")
Y75_84 = read_xlsx("Data/Combined_data/Y75_84_combined.xlsx")
Y_GE85 = read_xlsx("Data/Combined_data/Y_GE85_combined.xlsx")

## LC matrix (mortality rate)
col_name = c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")

Attiki = cbind(Y20_64$Attiki/Y20_64$Attiki_pop,Y65_74$Attiki/Y65_74$Attiki_pop, 
               Y75_84$Attiki/Y75_84$Attiki_pop, Y_GE85$Attiki/Y_GE85$Attiki_pop)
colnames(Attiki) = col_name

Lisbon = cbind(Y20_64$Lisbon/Y20_64$Lisbon_pop,Y65_74$Lisbon/Y65_74$Lisbon_pop, 
               Y75_84$Lisbon/Y75_84$Lisbon_pop, Y_GE85$Lisbon/Y_GE85$Lisbon_pop)
colnames(Lisbon) = col_name

Roma = cbind(Y20_64$Roma/Y20_64$Roma_pop,Y65_74$Roma/Y65_74$Roma_pop, 
             Y75_84$Roma/Y75_84$Roma_pop, Y_GE85$Roma/Y_GE85$Roma_pop)
colnames(Roma) = col_name

## Store heat/cold wave data
# Attiki
Attiki_wave = data.frame(cbind(Y20_64$Attiki_hot_wave3, Y20_64$Attiki_cold_wave3))
colnames(Attiki_wave) = c("Attiki_hot_wave3", "Attiki_cold_wave3")
# Lisbon
Lisbon_wave = data.frame(cbind(Y20_64$Lisbon_hot_wave3, Y20_64$Lisbon_cold_wave3))
colnames(Lisbon_wave) = c("Lisbon_hot_wave3", "Lisbon_cold_wave3")
# Roma
Roma_wave = data.frame(cbind(Y20_64$Roma_hot_wave3, Y20_64$Roma_cold_wave3))
colnames(Roma_wave) = c("Roma_hot_wave3", "Roma_cold_wave3")

# DLNM 
## Crossbasis matrix
load("Code/Crossbasis_matrix/cb_Attiki.RData")
load("Code/Crossbasis_matrix/cb_Lisbon.RData")
load("Code/Crossbasis_matrix/cb_Roma.RData")


## Convert to standard weekly mortality (*52)
Attiki = Attiki * 52
Lisbon = Lisbon * 52
Roma = Roma * 52
dat_list <- list(Attiki = Attiki, Lisbon = Lisbon, Roma = Roma)
wave_list <- list(Attiki = Attiki_wave, Lisbon = Lisbon_wave, Roma = Roma_wave)
cb_list <- list(Attiki = UTCI_ext.Attiki, Lisbon = UTCI_ext.Lisbon, Roma = UTCI_ext.Roma) 
### cb_list: for forecasting, we use ORIGINAL data! the parameter can be chosen in the forecasting function!


# Load functions
files <- list.files("Code/Function/", pattern = "\\.R$", full.names = TRUE)
invisible(sapply(files, source))


## Backfitting DLNM-LL model with expanding window cross validation
expanding_window_sets <- lapply(0:9, function(i) {
  create_train_test_sets(dat_list, wave_list, cb_list, train_size = 102 + i * 8, test_size = 180 + i * 8, lag = 21, var.df = 3, lag.df = 3, degree = 3)
})
names(expanding_window_sets) = paste0("set_", 1:10)


## Fitting + Forecasting  (Cross Validation)
# Iterate over expanding_window_sets and train the model, forecast, and extract test values
expanding_window_result <- lapply(expanding_window_sets, function(set) {
  # Extract the number of rows for the test set
  h_step <- nrow(set$test_set$dat_list[[1]])
  
  # Train model
  train_model <- B_DLNM_LL(dat_list = set$train_set$dat_list,
                           cb_list = set$train_set$cb_list,
                           wave_list = set$train_set$wave_list,
                           tol = 1e-2, max_iter = 20)
  
  # Forecast values using the dynamically determined number of periods
  model_forecast <- B_DLNM_LL.forecast(train_model,
                                       cb_matrix_test = set$test_set$cb_list,
                                       wave_data_test = set$test_set$wave_list,
                                       forecast_step = 78)
  
  # Log-transform test values
  test_values <- lapply(set$test_set$dat_list, function(region) {
    lapply(region, function(age_group) {
      log(age_group)
    })
  })
  
  # Assign names to test values to match the structure of train model outputs
  names(test_values) <- names(train_model$final_LL$log_fitted)
  
  # Return a list containing the train model, model forecast, and test values
  list(train_model = train_model,
       model_forecast = model_forecast,
       test_values = test_values)
})



# Forecasting mae function (calculate mae for each set by B-DLNM-LL vs baseline LL, and plot mae if necessary)
# Forecasting mae for each region each group (Backfitting DLNM-LL)
BDLNMLL_mae_list = lapply(expanding_window_result, function(set){
  ## Backfitting DLNM-LL
  forecast_BDLNMLL = lapply(set$model_forecast$B_DLNM_LL_final.forecast$model.forecast, function(region){
    as.data.frame(lapply(region, exp))
  })
  ## Test values
  test.values = lapply(set$test_values, function(region){
    as.data.frame(lapply(region, exp))
  })
  
  BDLNMLL_mae <- lapply(names(forecast_BDLNMLL), function(region) {
    sapply(names(forecast_BDLNMLL[[region]]), function(age_group) {
      mae(test.values[[region]][[age_group]], forecast_BDLNMLL[[region]][[age_group]])
    })
  })
  names(BDLNMLL_mae) <- names(forecast_BDLNMLL)
  BDLNMLL_mae <- lapply(BDLNMLL_mae, function(x) setNames(x, names(forecast_BDLNMLL[[1]])))
  return(BDLNMLL_mae)
})

# Forecasting mae for each region each group (Baseline LL)
baseline_LL_mae_list = lapply(expanding_window_result, function(set){
  ## Baseline LL
  forecast_baseline_LL = lapply(set$model_forecast$baseline_LL.forecast, function(region){
    as.data.frame(lapply(region, exp))
  })
  ## Test values
  test.values = lapply(set$test_values, function(region){
    as.data.frame(lapply(region, exp))
  })
  
  baseline_LL_mae <- lapply(names(forecast_baseline_LL), function(region) {
    sapply(names(forecast_baseline_LL[[region]]), function(age_group) {
      mae(test.values[[region]][[age_group]], forecast_baseline_LL[[region]][[age_group]])
    })
  })
  names(baseline_LL_mae) <- names(forecast_baseline_LL)
  baseline_LL_mae <- lapply(baseline_LL_mae, function(x) setNames(x, names(forecast_baseline_LL[[1]])))
  return(baseline_LL_mae)
})

extract_mae_line = function(mae_list){
  # Extract name from the list
  first_set <- names(mae_list)[1]
  first_region <- names(mae_list[[first_set]])[1]
  regions <- names(mae_list[[first_set]])
  age_groups <- names(mae_list[[first_set]][[first_region]])
  
  # Reconstruct the list
  BDLNMLL_mae_line  <- lapply(regions, function(region) {
    lapply(age_groups, function(age_group) {
      sapply(1:10, function(set_index) { # 10: for original cv, 12: by month
        mae_list[[paste0("set_", set_index)]][[region]][[age_group]]
      })
    })
  })
  names(BDLNMLL_mae_line) <- regions
  BDLNMLL_mae_line <- lapply(BDLNMLL_mae_line, function(age_group_list) {
    setNames(age_group_list, age_groups)
  })
}

BDLNMLL_mae_line = extract_mae_line(BDLNMLL_mae_list)
baseline_LL_mae_line = extract_mae_line(baseline_LL_mae_list)

## Compute average MAE
baseline_ave_mae = lapply(baseline_LL_mae_line, function(region){
  lapply(region, function(age_group){
    return(mean(age_group))
  })
})

B_DLNM_LL_ave_mae = lapply(BDLNMLL_mae_line, function(region){
  lapply(region, function(age_group){
    return(mean(age_group))
  })
})

# Attiki
round(100*unlist(baseline_ave_mae$Attiki),4)
round(100*unlist(B_DLNM_LL_ave_mae$Attiki),4)
# Lisbon
round(100*unlist(baseline_ave_mae$Lisbon),4)
round(100*unlist(B_DLNM_LL_ave_mae$Lisbon),4)
# Roma
round(100*unlist(baseline_ave_mae$Roma),4)
round(100*unlist(B_DLNM_LL_ave_mae$Roma),4)


# Model residual analysis (consider the whole dataset)
# Cross-basis matrix (smoothing)
# Attiki
varknots <- equalknots(UTCI_ext.Attiki, fun="ns",df = 3, degree = 3)
lagknots <- logknots(21,3)
Attiki_cb <- crossbasis(UTCI_ext.Attiki, lag=21, 
                        argvar=list(fun="ns", knots = varknots), 
                        arglag=list(knots=lagknots,df = 3))

# Lisbon
varknots <- equalknots(UTCI_ext.Lisbon,fun="ns",df = 3, degree = 3)
lagknots <- logknots(21, 3)
Lisbon_cb <- crossbasis(UTCI_ext.Lisbon, lag=21, 
                        argvar=list(fun="ns", knots = varknots), 
                        arglag=list(knots=lagknots,df = 3))

# Roma
varknots <- equalknots(UTCI_ext.Roma,fun="ns",df = 3, degree = 3)
lagknots <- logknots(21, 3)
Roma_cb <- crossbasis(UTCI_ext.Roma, lag=21, 
                      argvar=list(fun="ns", knots = varknots), 
                      arglag=list(knots=lagknots,df = 3))

cb_list <- list(Attiki = Attiki_cb, Lisbon = Lisbon_cb, Roma = Roma_cb) 

# run the model
DLNMLL = B_DLNM_LL(dat_list = dat_list, cb_list = cb_list, wave_list = wave_list, 
                   tol = 1e-2, max_iter = 20)

# Model residual
Attiki.residual = Attiki - exp(DLNMLL$log_fitted$Attiki)
Lisbon.residual = Lisbon - exp(DLNMLL$log_fitted$Lisbon)
Roma.residual = Roma - exp(DLNMLL$log_fitted$Roma)

## Seasonality test on model residual
## Model residual - Attiki
## QS test
Attiki.residual.QS.p <- sapply(Attiki.residual[age], function(residual) {
  round(qs(residual, freq = 52)$Pval, 3)
})

## Friedman test
Attiki.residual.Friedman.p <- sapply(Attiki.residual[age], function(residual) {
  round(fried(residual, freq = 52)$Pval, 3)
})

## Kruskall-Wallis test
Attiki.residual.KW.p <- sapply(Attiki.residual[age], function(residual) {
  round(kw(residual, freq = 52)$Pval, 3)
})


## Model residual - Lisbon
## QS test
Lisbon.residual.QS.p <- sapply(Lisbon.residual[age], function(residual) {
  round(qs(residual, freq = 52)$Pval, 3)
})

## Friedman test
Lisbon.residual.Friedman.p <- sapply(Lisbon.residual[age], function(residual) {
  round(fried(residual, freq = 52)$Pval, 3)
})

## Kruskall-Wallis test
Lisbon.residual.KW.p <- sapply(Lisbon.residual[age], function(residual) {
  round(kw(residual, freq = 52)$Pval, 3)
})


## Model residual - Roma
## QS test
Roma.residual.QS.p <- sapply(Roma.residual[age], function(residual) {
  round(qs(residual, freq = 52)$Pval, 3)
})

## Friedman test
Roma.residual.Friedman.p <- sapply(Roma.residual[age], function(residual) {
  round(fried(residual, freq = 52)$Pval, 3)
})

## Kruskall-Wallis test
Roma.residual.KW.p <- sapply(Roma.residual[age], function(residual) {
  round(kw(residual, freq = 52)$Pval, 3)
})


