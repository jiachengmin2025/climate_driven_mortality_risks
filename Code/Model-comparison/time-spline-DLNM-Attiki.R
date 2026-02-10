# Used Packages
packages <- c("readxl", "dlnm", "splines")
invisible(lapply(packages, library, character.only = TRUE))

# Load data
Y20_64 = read_xlsx("Data/Combined_data/Y20_64_combined.xlsx")
Y65_74 = read_xlsx("Data/Combined_data/Y65_74_combined.xlsx")
Y75_84 = read_xlsx("Data/Combined_data/Y75_84_combined.xlsx")
Y_GE85 = read_xlsx("Data/Combined_data/Y_GE85_combined.xlsx")


# Weekly mortality rate
Attiki = cbind(Y20_64[,3]/Y20_64[,9],Y65_74[,3]/Y65_74[,9], 
               Y75_84[,3]/Y75_84[,9], Y_GE85[,3]/Y_GE85[,9])
colnames(Attiki) = c('Y20_64', 'Y65_74', 'Y75_84', 'Y_GE85')
rownames(Attiki) = Y20_64$Week
Attiki = Attiki * 52

# Load wave data
Attiki_wave = data.frame(cbind(Y20_64$Attiki_hot_wave3, Y20_64$Attiki_cold_wave3))
colnames(Attiki_wave) = c("Attiki_hot_wave3", "Attiki_cold_wave3")


# Load functions
files <- list.files("Code/Function/", pattern = "\\.R$", full.names = TRUE)
invisible(sapply(files, source))


## Load cross-basis matrix
load("Code/Crossbasis_matrix/cb_Attiki.RData")


# Expanding window sets
sets = 10
size = 8
forecast_step = 78
start_ind = 102

expanding_window_sets <- lapply(0:(sets-1), function(i) {
  # Training set
  Attiki_train <- Attiki[1:(start_ind + i * size), ]
  # add time index
  Attiki_train <- cbind(time = 1:nrow(Attiki_train), Attiki_train)
  Attiki_wave_train <- Attiki_wave[1:(start_ind + i * size), ]
  UTCI_ext.Attiki_train <- UTCI_ext.Attiki[1:(start_ind + i * size), ]
  
  # Smoothing cross-basis matrix (train)
  varknots <- equalknots(UTCI_ext.Attiki_train, fun = "ns", df = 3, degree = 3)
  lagknots <- logknots(21, 3)
  Attiki_cb_train <- crossbasis(UTCI_ext.Attiki_train, lag = 21, 
                                argvar = list(fun = "ns", knots = varknots), 
                                arglag = list(knots = lagknots, df = 3))
  
  # Test set
  Attiki_test <- Attiki[(start_ind + i * size + 1):(start_ind + i * size + forecast_step), ]
  # add time index
  Attiki_test <- cbind(time = 1:nrow(Attiki_test), Attiki_test)
  Attiki_wave_test <- Attiki_wave[(start_ind + i * size + 1):(start_ind + i * size + forecast_step), ]
  UTCI_ext.Attiki_test <- UTCI_ext.Attiki[(start_ind + i * size + 1):(start_ind + i * size + forecast_step), ]
  
  # Smoothing cross-basis matrix (test) *** using training set setting ***
  varknots <- equalknots(UTCI_ext.Attiki_train, fun = "ns", df = 3, degree = 3)
  lagknots <- logknots(21, 3)
  Attiki_cb_test <- crossbasis(UTCI_ext.Attiki_test, lag = 21, 
                               argvar = list(fun = "ns", knots = varknots), 
                               arglag = list(knots = lagknots, df = 3)) 
  
  # Return a list of training and test sets for this iteration
  list(
    train = list(Attiki_train = Attiki_train, Attiki_wave_train = Attiki_wave_train,
                 Attiki_cb_train = Attiki_cb_train),
    test = list(Attiki_test = Attiki_test, Attiki_wave_test = Attiki_wave_test,
                Attiki_cb_test = Attiki_cb_test)
  )
})
names(expanding_window_sets) = paste0("set_", 1:sets)


# DLNM MAE

dlnm_mae <- function(sets, set_index, age_col, df_time){
  sname <- paste0("set_", set_index)
  tr <- sets[[sname]]$train
  te <- sets[[sname]]$test
  
  # train
  y_tr <- tr$Attiki_train[[age_col]]
  time_tr <- tr$Attiki_train$time
  hot_tr <- tr$Attiki_wave_train$Attiki_hot_wave3
  cold_tr <- tr$Attiki_wave_train$Attiki_cold_wave3
  cb_tr   <- tr$Attiki_cb_train
  
  train_df <- data.frame(
    y = y_tr, time = time_tr, hot = hot_tr, cold = cold_tr
  )
  train_df$cb <- cb_tr
  fit <- glm(y ~ cb + hot + cold + ns(time, df = df_time), 
             family = quasipoisson(), data = train_df)
  fitted_tr <- fitted(fit)
  
  # test
  y_te <- te$Attiki_test[[age_col]]
  time_te <- te$Attiki_test$time
  hot_te  <- te$Attiki_wave_test$Attiki_hot_wave3
  cold_te <- te$Attiki_wave_test$Attiki_cold_wave3
  cb_te <- te$Attiki_cb_test
  
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





