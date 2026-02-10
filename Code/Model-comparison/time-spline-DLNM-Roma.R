# Used Packages
packages <- c("readxl", "dlnm", "splines")
invisible(lapply(packages, library, character.only = TRUE))

# Load data
Y20_64 = read_xlsx("Data/Combined_data/Y20_64_combined.xlsx")
Y65_74 = read_xlsx("Data/Combined_data/Y65_74_combined.xlsx")
Y75_84 = read_xlsx("Data/Combined_data/Y75_84_combined.xlsx")
Y_GE85 = read_xlsx("Data/Combined_data/Y_GE85_combined.xlsx")


# Weekly mortality rate
Roma = cbind(Y20_64[,6]/Y20_64[,12],Y65_74[,6]/Y65_74[,12], 
               Y75_84[,6]/Y75_84[,12], Y_GE85[,6]/Y_GE85[,12])
colnames(Roma) = c('Y20_64', 'Y65_74', 'Y75_84', 'Y_GE85')
rownames(Roma) = Y20_64$Week
Roma = Roma * 52

# Load wave data
Roma_wave = data.frame(cbind(Y20_64$Roma_hot_wave3, Y20_64$Roma_cold_wave3))
colnames(Roma_wave) = c("Roma_hot_wave3", "Roma_cold_wave3")


# Load functions
files <- list.files("Code/Function/", pattern = "\\.R$", full.names = TRUE)
invisible(sapply(files, source))


## Load cross-basis matrix
load("Code/Crossbasis_matrix/cb_Roma.RData")


# Expanding window sets
sets = 10
size = 8
forecast_step = 78
start_ind = 102

expanding_window_sets <- lapply(0:(sets-1), function(i) {
  # Training set
  Roma_train <- Roma[1:(start_ind + i * size), ]
  # add time index
  Roma_train <- cbind(time = 1:nrow(Roma_train), Roma_train)
  Roma_wave_train <- Roma_wave[1:(start_ind + i * size), ]
  UTCI_ext.Roma_train <- UTCI_ext.Roma[1:(start_ind + i * size), ]
  
  # Smoothing cross-basis matrix (train)
  varknots <- equalknots(UTCI_ext.Roma_train, fun = "ns", df = 3, degree = 3)
  lagknots <- logknots(21, 3)
  Roma_cb_train <- crossbasis(UTCI_ext.Roma_train, lag = 21, 
                                argvar = list(fun = "ns", knots = varknots), 
                                arglag = list(knots = lagknots, df = 3))
  
  # Test set
  Roma_test <- Roma[(start_ind + i * size + 1):(start_ind + i * size + forecast_step), ]
  # add time index
  Roma_test <- cbind(time = 1:nrow(Roma_test), Roma_test)
  Roma_wave_test <- Roma_wave[(start_ind + i * size + 1):(start_ind + i * size + forecast_step), ]
  UTCI_ext.Roma_test <- UTCI_ext.Roma[(start_ind + i * size + 1):(start_ind + i * size + forecast_step), ]
  
  # Smoothing cross-basis matrix (test) *** using training set setting ***
  varknots <- equalknots(UTCI_ext.Roma_train, fun = "ns", df = 3, degree = 3)
  lagknots <- logknots(21, 3)
  Roma_cb_test <- crossbasis(UTCI_ext.Roma_test, lag = 21, 
                               argvar = list(fun = "ns", knots = varknots), 
                               arglag = list(knots = lagknots, df = 3)) 
  
  # Return a list of training and test sets for this iteration
  list(
    train = list(Roma_train = Roma_train, Roma_wave_train = Roma_wave_train,
                 Roma_cb_train = Roma_cb_train),
    test = list(Roma_test = Roma_test, Roma_wave_test = Roma_wave_test,
                Roma_cb_test = Roma_cb_test)
  )
})
names(expanding_window_sets) = paste0("set_", 1:sets)


# DLNM MAE

dlnm_mae <- function(sets, set_index, age_col, df_time){
  sname <- paste0("set_", set_index)
  tr <- sets[[sname]]$train
  te <- sets[[sname]]$test
  
  # train
  y_tr <- tr$Roma_train[[age_col]]
  time_tr <- tr$Roma_train$time
  hot_tr <- tr$Roma_wave_train$Roma_hot_wave3
  cold_tr <- tr$Roma_wave_train$Roma_cold_wave3
  cb_tr   <- tr$Roma_cb_train
  
  train_df <- data.frame(
    y = y_tr, time = time_tr, hot = hot_tr, cold = cold_tr
  )
  train_df$cb <- cb_tr
  fit <- glm(y ~ cb + hot + cold + ns(time, df = df_time), 
             family = quasipoisson(), data = train_df)
  fitted_tr <- fitted(fit)
  
  # test
  y_te <- te$Roma_test[[age_col]]
  time_te <- te$Roma_test$time
  hot_te  <- te$Roma_wave_test$Roma_hot_wave3
  cold_te <- te$Roma_wave_test$Roma_cold_wave3
  cb_te <- te$Roma_cb_test
  
  test_df <- data.frame(
    time = time_te, hot = hot_te, cold = cold_te
  )
  test_df$cb <- cb_te
  
  pred <- predict(fit, newdata = test_df, type = "response", se.fit = TRUE)
  mae  <- mean(abs(y_te - pred$fit), na.rm = TRUE)
  
  # output
  list(model = fit, fitted_train = fitted_tr, pred_test = pred$fit, mae = mae)
}



res_list.20 <- lapply(1:10, function(i) dlnm_mae(expanding_window_sets, i, age_col = "Y20_64", df_time = 12))
mean_mae.20 <- 100*mean(sapply(res_list.20, `[[`, "mae"))

res_list.65 <- lapply(1:10, function(i) dlnm_mae(expanding_window_sets, i, age_col = "Y65_74", df_time = 12))
mean_mae.65 <- 100*mean(sapply(res_list.65, `[[`, "mae"))

res_list.75 <- lapply(1:10, function(i) dlnm_mae(expanding_window_sets, i, age_col = "Y75_84", df_time = 12))
mean_mae.75 <- 100*mean(sapply(res_list.75, `[[`, "mae"))

res_list.85 <- lapply(1:10, function(i) dlnm_mae(expanding_window_sets, i, age_col = "Y_GE85", df_time = 12))
mean_mae.85 <- 100*mean(sapply(res_list.85, `[[`, "mae"))

round(c(mean_mae.20, mean_mae.65, mean_mae.75, mean_mae.85), 4)





