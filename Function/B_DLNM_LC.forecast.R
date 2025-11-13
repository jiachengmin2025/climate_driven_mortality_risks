B_DLNM_LC.forecast = function(B_DLNM_LC_model, cb_matrix_test, wave_data_test, kt_model_type, forecast_step){
  if (nrow(cb_matrix_test) != forecast_step){
    stop("The dimension of new data does not match the forecasting step.")
  }
  # Forecasting parameters from Lee-Carter model
  axh = B_DLNM_LC_model$final_LC$a_x
  bxh = B_DLNM_LC_model$final_LC$b_x
  kt.fit = kt_model(B_DLNM_LC_model$final_LC$k_t, kt_model_type)
  kth = c(forecast(kt.fit, h = forecast_step)$mean)
  # Lee-Carter part forecasting
  # log(m(x, t+h|t)) = a(x|t) + b(x|t) %*% k(t+h|t) 
  bxh_kth = outer(kth, bxh, FUN = "*")
  log_LC.forecast = sweep(bxh_kth, 2, axh, FUN = "+")
  
  # DLNM part forecasting
  ## (1) dlnm_coef %*% smoothing crossbasis matrix (each iteration each age group)
  cb_matrix.forecast <- lapply(B_DLNM_LC_model$dlnm_coef, function(iteration) {
    lapply(iteration, function(age_group) {
      apply(cb_matrix_test, 1, function(row) {
        sum(row * age_group[2:(length(age_group)-2)]) + age_group[1] # intercept
      })
    })
  })
  ## sum up all forecast values from DLNMs.
  cb_matrix.sum.forecast = Reduce(function(x, y) {
    mapply(`+`, x, y, SIMPLIFY = FALSE)
  }, cb_matrix.forecast)
  cb_matrix.sum.forecast = do.call(cbind, cb_matrix.sum.forecast)
  
  
  ## (2) dlnm_coef %*% hot/cold wave term (each iteration each age group)
  wave.forecast = lapply(B_DLNM_LC_model$dlnm_coef, function(iteration){
    lapply(iteration, function(age_group){
      apply(wave_data_test, 1,
            function(row) sum(row * age_group[(length(age_group)-1):length(age_group)])) 
    })
  })
  wave.sum.forecast = Reduce(function(x, y) {
    # sum by each age group
    mapply(`+`, x, y, SIMPLIFY = FALSE)
  }, wave.forecast)
  wave.sum.forecast = do.call(cbind, wave.sum.forecast)
  
  # Combine Lee--Carter forecasting part, DLNM forecasting part and hot/cold wave part
  model.forecast = exp(log_LC.forecast + cb_matrix.sum.forecast + wave.sum.forecast)
  
  # Baseline LC part
  ## Baseline forecasting
  ax_ref = B_DLNM_LC_model$baseline_LC$a_x
  bx_ref = B_DLNM_LC_model$baseline_LC$b_x
  kt_ref.fit = kt_model(B_DLNM_LC_model$baseline_LC$k_t, kt_model_type)
  kt_ref = c(forecast(kt_ref.fit, h = forecast_step)$mean)
  bxkt_ref = outer(kt_ref, bx_ref, FUN = "*")
  log_LC_baseline.forecast = sweep(bxkt_ref, 2, ax_ref, FUN = "+")
  
  # Output
  output = list(baseline_LC.parameter = list(ax_ref = ax_ref, bx_ref = bx_ref, 
                                             kt_ref.fit = kt_ref.fit, kt_ref = kt_ref),
                baseline_LC.forecast = exp(log_LC_baseline.forecast),
                B_DLNM_LC_final.parameter = list(ax = axh, bx = bxh, kt.fit = kt.fit,
                                                 kt = kth),
                B_DLNM_LC_final.forecast = list(model.forecast = model.forecast,
                                                converged_LC.forecast = exp(log_LC.forecast),
                                                cb_matrix.forecast = cb_matrix.forecast,
                                                cb_matrix.sum.forecast = cb_matrix.sum.forecast,
                                                wave.forecast = wave.forecast,
                                                wave.sum.forecast = wave.sum.forecast))
  return(output)
}