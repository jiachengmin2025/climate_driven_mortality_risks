## Used Packages
packages <- c("readxl", "dlnm", "splines", "mgcv", "forecast", "demography", 
              "Metrics", "tseries", "nortest", "seastests")
invisible(lapply(packages, library, character.only = TRUE))


## Load data
Y20_64 = read_xlsx("Data/Combined_data/Y20_64_combined.xlsx")
Y65_74 = read_xlsx("Data/Combined_data/Y65_74_combined.xlsx")
Y75_84 = read_xlsx("Data/Combined_data/Y75_84_combined.xlsx")
Y_GE85 = read_xlsx("Data/Combined_data/Y_GE85_combined.xlsx")

## LC matrix (mortality rate), convert to weekly mortality rate
Lisbon = cbind(Y20_64$Lisbon/Y20_64$Lisbon_pop,Y65_74$Lisbon/Y65_74$Lisbon_pop, 
               Y75_84$Lisbon/Y75_84$Lisbon_pop, Y_GE85$Lisbon/Y_GE85$Lisbon_pop)
colnames(Lisbon) = c('Y20_64', 'Y65_74', 'Y75_84', 'Y_GE85')
rownames(Lisbon) = Y20_64$Week
Lisbon = Lisbon * 52

## Store heat/cold wave data and load crossbasis matrix data
Lisbon_wave = data.frame(cbind(Y20_64$Lisbon_hot_wave3, Y20_64$Lisbon_cold_wave3))
colnames(Lisbon_wave) = c("Lisbon_hot_wave3", "Lisbon_cold_wave3")
# crossbasis matrix (before smoothing)
load("Code/Crossbasis_matrix/cb_Lisbon.RData")

# Load functions
files <- list.files("Code/Function/", pattern = "\\.R$", full.names = TRUE)
sapply(files, source)

#  Expanding window cross-validation
sets = 10
size = 8
forecast_step = 78
start_ind = 102

expanding_window_sets <- lapply(0:(sets-1), function(i) {
  # Training set
  Lisbon_train <- Lisbon[1:(start_ind + i * size), ]
  Lisbon_wave_train <- Lisbon_wave[1:(start_ind + i * size), ]
  UTCI_ext.Lisbon_train <- UTCI_ext.Lisbon[1:(start_ind + i * size), ]
  
  # Smoothing cross-basis matrix (train)
  varknots <- equalknots(UTCI_ext.Lisbon_train, fun = "ns", df = 3, degree = 3)
  lagknots <- logknots(21, 3)
  Lisbon_cb_train <- crossbasis(UTCI_ext.Lisbon_train, lag = 21, 
                                argvar = list(fun = "ns", knots = varknots), 
                                arglag = list(knots = lagknots, df = 3))
  
  # Test set
  Lisbon_test <- Lisbon[(start_ind + i * size + 1):(start_ind + i * size + forecast_step), ]
  Lisbon_wave_test <- Lisbon_wave[(start_ind + i * size + 1):(start_ind + i * size + forecast_step), ]
  UTCI_ext.Lisbon_test <- UTCI_ext.Lisbon[(start_ind + i * size + 1):(start_ind + i * size + forecast_step), ]
  
  # Smoothing cross-basis matrix (test) *** using training set setting ***
  varknots <- equalknots(UTCI_ext.Lisbon_train, fun = "ns", df = 3, degree = 3)
  lagknots <- logknots(21, 3)
  Lisbon_cb_test <- crossbasis(UTCI_ext.Lisbon_test, lag = 21, 
                               argvar = list(fun = "ns", knots = varknots), 
                               arglag = list(knots = lagknots, df = 3)) 
  
  # Return a list of training and test sets for this iteration
  list(
    train = list(Lisbon_train = Lisbon_train, Lisbon_wave_train = Lisbon_wave_train,
                 Lisbon_cb_train = Lisbon_cb_train),
    test = list(Lisbon_test = Lisbon_test, Lisbon_wave_test = Lisbon_wave_test,
                Lisbon_cb_test = Lisbon_cb_test)
  )
})
names(expanding_window_sets) = paste0("set_", 1:sets)

# Apply the B_DLNM_LC fitting and forecasting across all sets
expanding_window_result <- lapply(expanding_window_sets, function(set) {
  # Train the model on the training set
  train_model <- B_DLNM_LC(set$train$Lisbon_train, 
                           set$train$Lisbon_cb_train, 
                           set$train$Lisbon_wave_train, 
                           tol = 1e-2, max_iter = 20, region = "Lisbon") ## change 
  
  # Forecast using the trained model on the test set (cross-basis matrix and hot/cold wave data)
  model_forecast <- B_DLNM_LC.forecast(train_model, 
                                       set$test$Lisbon_cb_test, 
                                       set$test$Lisbon_wave_test, 
                                       kt_model_type = "SARIMA", 
                                       forecast_step = forecast_step)
  
  # Return both the train model and the forecast for this iteration
  list(train_model = train_model, model_forecast = model_forecast,
       Lisbon_test = set$test$Lisbon_test)
})



# Forecasting mae function (calculate mae for each set by B-DLNM-LC vs baseline LC, and plot mae if necessary)
expanding_window_forecast_mae <- function(set) {
  # Test values
  test_value <- exp(set$Lisbon_test)
  # B-DLNM-LC model forecast
  forecast_BDLNMLC <- exp(set$model_forecast$B_DLNM_LC_final.forecast$model.forecast)
  B_DLNM_LC_mae <- c(100*mae(test_value[,1], forecast_BDLNMLC[,1]),
                     100*mae(test_value[,2], forecast_BDLNMLC[,2]),
                     100*mae(test_value[,3], forecast_BDLNMLC[,3]),
                     100*mae(test_value[,4], forecast_BDLNMLC[,4]))
  
  # Baseline LC model forecast
  baseline_LC <- exp(set$model_forecast$baseline_LC.forecast)
  baseline_LC_mae <- c(100*mae(test_value[,1], baseline_LC[,1]),
                       100*mae(test_value[,2], baseline_LC[,2]),
                       100*mae(test_value[,3], baseline_LC[,3]),
                       100*mae(test_value[,4], baseline_LC[,4]))
  
  # Converged LC model forecast
  converged_LC <- exp(set$model_forecast$B_DLNM_LC_final.forecast$converged_LC.forecast)
  converged_LC_mae <- c(100*mae(test_value[,1], converged_LC[,1]),
                        100*mae(test_value[,2], converged_LC[,2]),
                        100*mae(test_value[,3], converged_LC[,3]),
                        100*mae(test_value[,4], converged_LC[,4]))
  # Return mae values and the plot
  return(list(B_DLNM_LC_mae = B_DLNM_LC_mae, baseline_LC_mae = baseline_LC_mae, converged_LC_mae = converged_LC_mae))
}

# Apply the function across all sets
forecast_mae <- lapply(expanding_window_result, function(set){
  expanding_window_forecast_mae(set)
  
})

## Compute average Forecasting mse
B_DLNM_LC_mae_values <- do.call(rbind, lapply(forecast_mae, function(set) set$B_DLNM_LC_mae))
baseline_LC_mae_values <- do.call(rbind, lapply(forecast_mae, function(set) set$baseline_LC_mae))
converged_LC_mae_values <- do.call(rbind, lapply(forecast_mae, function(set) set$converged_LC_mae))
age_groups <- c("20-64", "65-74", "75-84", "85+")
round((baseline_LC_ave_mae = colMeans(baseline_LC_mae_values)), 4)
round((B_DLNM_LC_ave_mae = colMeans(B_DLNM_LC_mae_values)), 4)

# Model residual analysis (consider the whole dataset)

## Smoothing cross-basis matrix on the whole dataset
varknots <- equalknots(UTCI_ext.Lisbon, fun = "ns", df = 3, degree = 3)
lagknots <- logknots(21, 3)
Lisbon_cb <- crossbasis(UTCI_ext.Lisbon, lag = 21, 
                        argvar = list(fun = "ns", knots = varknots), 
                        arglag = list(knots = lagknots, df = 3))

Lisbon.DLNMLC <- B_DLNM_LC(Lisbon, Lisbon_cb, Lisbon_wave, 
                           tol = 1e-2, max_iter = 20, region = "Lisbon")
Lisbon.residual = Lisbon - exp(Lisbon.DLNMLC$log_fitted)

## Normality test (Shapiro–Wilk test, Anderson–Darling test, and Jarque–Bera test)
### In Shapiro–Wilk test, Anderson–Darling test, and Jarque–Bera test, H0: residual is Gaussian vs H1: residual is non-Gaussian. 
# Shapiro–Wilk test
shapiro.p.value <- sapply(Lisbon.residual, function(x) shapiro.test(x[is.finite(x)])$p.value)
ifelse(shapiro.p.value < 0.05, "Residual is non-Gaussian", "Residual is Guassian")
# Anderson–Darling test
ad.p.value <- sapply(Lisbon.residual, function(x) ad.test(x[is.finite(x)])$p.value)
ifelse(ad.p.value < 0.05, "Residual is non-Gaussian", "Residual is Guassian")
# Jarque–Bera test
jarque.bera.p.value <- sapply(Lisbon.residual, function(x) jarque.bera.test(x[is.finite(x)])$p.value)
ifelse(jarque.bera.p.value < 0.05, "Residual is non-Gaussian", "Residual is Guassian")

## Seasonality test on model residual
# Lisbon 
## QS test
Lisbon.residual.QS.p <- sapply(Lisbon.residual, function(residual) {
  round(qs(residual, freq = 52)$Pval, 3)
})

## Friedman test
Lisbon.residual.Friedman.p <- sapply(Lisbon.residual, function(residual) {
  round(fried(residual, freq = 52)$Pval, 3)
})

## Kruskall-Wallis test
Lisbon.residual.KW.p <- sapply(Lisbon.residual, function(residual) {
  round(kw(residual, freq = 52)$Pval, 3)
})




