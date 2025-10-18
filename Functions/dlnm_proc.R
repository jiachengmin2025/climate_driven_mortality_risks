dlnm_proc = function(dat, cb_train, wave_data, region){
  hot_wave_var = paste0(region, "_hot_wave3")
  cold_wave_var = paste0(region, "_cold_wave3")
  
  # 20-64
  DLNM20 = glm(exp(dat[,1])  ~ cb_train + wave_data[[hot_wave_var]] + wave_data[[cold_wave_var]], 
               family = quasipoisson())
  DLNM20_fitted = log(fitted.values(DLNM20))
  
  # 65-74
  DLNM65 = glm(exp(dat[,2])  ~ cb_train + wave_data[[hot_wave_var]] + wave_data[[cold_wave_var]], 
               family = quasipoisson())
  DLNM65_fitted = log(fitted.values(DLNM65))
  
  # 75-84
  DLNM75 = glm(exp(dat[,3])  ~ cb_train + wave_data[[hot_wave_var]] + wave_data[[cold_wave_var]], 
               family = quasipoisson())
  DLNM75_fitted = log(fitted.values(DLNM75))
  
  # 85+
  DLNM85 = glm(exp(dat[,4])  ~ cb_train + wave_data[[hot_wave_var]] + wave_data[[cold_wave_var]], 
               family = quasipoisson())
  DLNM85_fitted = log(fitted.values(DLNM85))
  
  DLNM_result = list(DLNM_fitted = data.frame(DLNM20_fitted, DLNM65_fitted, DLNM75_fitted, DLNM85_fitted),
                     DLNM20 = DLNM20, DLNM65 = DLNM65, DLNM75 = DLNM75, DLNM85 = DLNM85,
                     DLNM_coef = list(Y20_64 = coef(DLNM20), Y65_74 = coef(DLNM65), 
                                      Y75_84 = coef(DLNM75), Y_GE85 = coef(DLNM85)),
                     DLNM_vcov = list(Y20_64 = vcov(DLNM20), Y65_74 = vcov(DLNM65), 
                                      Y75_84 = vcov(DLNM75), Y_GE85 = vcov(DLNM85)))
  
  return(DLNM_result)
}