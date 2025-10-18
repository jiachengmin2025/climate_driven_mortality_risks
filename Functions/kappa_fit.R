kappa_fit = function(B_DLNM_LL_model, forecast_step){
  ## Final model
  Kt.df = data.frame(B_DLNM_LL_model$final_LL$Kt, B_DLNM_LL_model$final_LL$kt[[1]], 
                     B_DLNM_LL_model$final_LL$kt[[2]], B_DLNM_LL_model$final_LL$kt[[3]])
  names(Kt.df) = c("Common", "Attiki", "Lisbon", "Roma")
  
  # Common
  # Common_fit = arima(ts(Kt.df$Common, frequency = forecast_step), order = c(0,1,0)) # RWD
  Common_fit = auto.arima(ts(Kt.df$Common, frequency = 52), seasonal = T, max.d = 0, allowmean = F) # ARIMA
  # summary(Common_fit)
  Common_forecast = forecast(Common_fit, h = forecast_step)
  Common_pred = Common_forecast$mean
  # ts.plot(ts(Kt.df$Common, frequency = forecast_step), xlim = c(1,6))
  # lines(Common_pred, col = "red")
  
  # Attiki
  # Attiki_fit = arima(ts(Kt.df$Attiki, frequency = forecast_step), order = c(0,1,0))
  Attiki_fit = auto.arima(ts(Kt.df$Attiki, frequency = 52), seasonal = T, max.d = 0, allowmean = F)
  # summary(Attiki_fit)
  Attiki_forecast = forecast(Attiki_fit, h = forecast_step)
  Attiki_pred = Attiki_forecast$mean
  # ts.plot(ts(Kt.df$Attiki, frequency = forecast_step), xlim = c(1,6))
  # lines(Attiki_pred, col = "red")
  
  # Lisbon 
  # Lisbon_fit = arima(ts(Kt.df$Lisbon, frequency = forecast_step), order = c(0,1,0))
  Lisbon_fit = auto.arima(ts(Kt.df$Lisbon, frequency = 52), seasonal = T, max.d = 0, allowmean = F)
  # summary(Lisbon_fit)
  Lisbon_forecast = forecast(Lisbon_fit, h = forecast_step)
  Lisbon_pred = Lisbon_forecast$mean
  # ts.plot(ts(Kt.df$Lisbon, frequency = forecast_step), xlim = c(1,6))
  # lines(Lisbon_pred, col = "red")
  
  # Roma
  # Roma_fit = arima(ts(Kt.df$Roma, frequency = forecast_step), order = c(0,1,0))
  Roma_fit = auto.arima(ts(Kt.df$Roma, frequency = 52), seasonal = T, max.d = 0, allowmean = F)
  # summary(Roma_fit)
  Roma_forecast = forecast(Roma_fit, h = forecast_step)
  Roma_pred = Roma_forecast$mean
  # ts.plot(ts(Kt.df$Roma, frequency = forecast_step), xlim = c(1,6))
  # lines(Roma_pred, col = "red")
  
  kt.fit = list(Common_fit = Common_fit, Attiki_fit = Attiki_fit,
                Lisbon_fit = Lisbon_fit, Roma_fit = Roma_fit)
  kt.forecast = list(Common_forecast = Common_forecast, Attiki_forecast = Attiki_forecast, 
                     Lisbon_forecast = Lisbon_forecast, Roma_forecast = Roma_forecast)
  kt.pred = list(Common_pred = Common_pred, Attiki_pred = Attiki_pred,
                 Lisbon_pred = Lisbon_pred, Roma_pred = Roma_pred)
  
  ## Baseline model
  Kt.df_ref = data.frame(B_DLNM_LL_model$baseline_LL$Kt, B_DLNM_LL_model$baseline_LL$kt[[1]], 
                         B_DLNM_LL_model$baseline_LL$kt[[2]], B_DLNM_LL_model$baseline_LL$kt[[3]])
  names(Kt.df_ref) = c("Common", "Attiki", "Lisbon", "Roma")
  
  # Common
  Common_fit_ref = auto.arima(ts(Kt.df_ref$Common, frequency = 52), seasonal = T) # SARIMA
  # summary(Common_fit_ref)
  Common_forecast_ref = forecast(Common_fit_ref, h = forecast_step)
  Common_pred_ref = Common_forecast_ref$mean
  # ts.plot(ts(Kt.df_ref$Common, frequency = forecast_step), xlim = c(1,6))
  # lines(Common_pred_ref, col = "red")
  
  # Attiki
  # Attiki_fit_ref = arima(ts(Kt.df_ref$Attiki, frequency = 52), order = c(0,1,0))
  Attiki_fit_ref = auto.arima(ts(Kt.df_ref$Attiki, frequency = 52), seasonal = T) # SARIMA
  # summary(Attiki_fit_ref)
  Attiki_forecast_ref = forecast(Attiki_fit_ref, h = forecast_step)
  Attiki_pred_ref = Attiki_forecast_ref$mean

  # Lisbon
  # Lisbon_fit_ref = arima(ts(Kt.df_ref$Lisbon, frequency = 52),  order = c(0,1,0))
  Lisbon_fit_ref = auto.arima(ts(Kt.df_ref$Attiki, frequency = 52), seasonal = T) # SARIMA
  # summary(Lisbon_fit_ref)
  Lisbon_forecast_ref = forecast(Lisbon_fit_ref, h = forecast_step)
  Lisbon_pred_ref = Lisbon_forecast_ref$mean
  
  # Roma
  # Roma_fit_ref = arima(ts(Kt.df_ref$Roma, frequency = 52),  order = c(0,1,0))
  Roma_fit_ref = auto.arima(ts(Kt.df_ref$Attiki, frequency = 52), seasonal = T) # SARIMA
  # summary(Roma_fit_ref)
  Roma_forecast_ref = forecast(Roma_fit_ref, h = forecast_step)
  Roma_pred_ref = Roma_forecast_ref$mean
  
  kt.fit_ref = list(Common_fit_ref = Common_fit_ref, Attiki_fit_ref = Attiki_fit_ref,
                    Lisbon_fit_ref = Lisbon_fit_ref, Roma_fit_ref = Roma_fit_ref)
  kt.forecast_ref = list(Common_forecast_ref = Common_forecast_ref, Attiki_forecast_ref = Attiki_forecast_ref, 
                         Lisbon_forecast_ref = Lisbon_forecast_ref, Roma_forecast_ref = Roma_forecast_ref)
  kt.pred_ref = list(Common_pred_ref = Common_pred_ref, Attiki_pred_ref = Attiki_pred_ref,
                     Lisbon_pred_ref = Lisbon_pred_ref, Roma_pred_ref = Roma_pred_ref)
  
  kt_set = list(final_kt = list(kt.fit = kt.fit, kt.forecast = kt.forecast, kt.pred = kt.pred), 
                base_kt = list(kt.fit_ref = kt.fit_ref, kt.forecast_ref = kt.forecast_ref, kt.pred_ref = kt.pred_ref))
  return(kt_set)
}