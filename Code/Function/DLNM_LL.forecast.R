#' DLNM_LL.forecast: Forecast mortality using a fitted DLNM-LL model.
#'
#' @param DLNM_LL_model Output object from \code{DLNM_LL()}, containing \code{final_LL}, \code{baseline_LL},
#'   and region-wise iteration lists \code{dlnm_coef}.
#' @param cb_matrix_test Named list of cross-basis design matrices for the forecast horizon, one per region.
#' @param wave_data_test Optional named list of wave covariate matrices/data.frames
#'   for the forecast horizon, one per region. If a region entry is \code{NULL}
#'   or empty, the wave contribution is set to zero for that region.
#' @param forecast_step Integer forecast horizon (number of steps ahead).
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{baseline_LL.parameter}: baseline LL parameters and and fitted time-series model.
#'   \item \code{baseline_LL.forecast}: baseline LL forecasts (log scale in current implementation).
#'   \item \code{DLNM_LL_final.parameter}: final (converged in algorithm) LL parameters.
#'   \item \code{DLNM_LL_final.forecast}: forecasts from the DLNM-LL model.
#' }

DLNM_LL.forecast = function(DLNM_LL_model, cb_matrix_test, wave_data_test, forecast_step){
  region_labels = if (exists("region_name", inherits = TRUE)) {
    region_name
  } else {
    names(DLNM_LL_model$dlnm_coef)
  }
  cb_matrix_test = lapply(cb_matrix_test, as.matrix)
  if (is.null(wave_data_test)) {
    wave_data_test = vector("list", length(cb_matrix_test))
  }
#' DLNM-LL forecast internal: retrieve region-specific wave covariates
#'
#' @param i Integer region index.
#'
#' @return Wave covariate matrix/data.frame for region \code{i}, or \code{NULL}
#'   if no wave covariates are available.
  get_region_wave = function(i) {
    region_key = names(DLNM_LL_model$dlnm_coef)[i]
    if (!is.null(names(wave_data_test)) && region_key %in% names(wave_data_test)) {
      return(wave_data_test[[region_key]])
    }
    if (length(wave_data_test) >= i) {
      return(wave_data_test[[i]])
    }
    NULL
  }

  Axh = DLNM_LL_model$final_LL$Ax # list
  Bxh = DLNM_LL_model$final_LL$Bx 
  bxh = DLNM_LL_model$final_LL$bx # list 
  # Fitting and forecasting Kt sets
  kt_fitting = kappa_fit(DLNM_LL_model = DLNM_LL_model, forecast_step = forecast_step)
  kth = kt_fitting$final_kt$kt.pred[-1]
  Kth = kt_fitting$final_kt$kt.pred[[1]]
  
  # Lee--Li part forecasting
  # log(m(x, t+h|t)) = a(x|t) + B(x|t) %*% K(t+h|t) + b(x,i|t) %*% k(t+h,i|t) 
  ## B(x|t) %*% K(t+h|t)
  Bxh_Kth = outer(Kth, Bxh, FUN = "*")
  ## b(x,i|t) %*% k(t+h,i|t) (list of n data frame)
  bxh_kth = mapply(function(kth.i, bxh.i) outer(kth.i, bxh.i, FUN = "*"), kth, bxh, SIMPLIFY = FALSE)
  ## B(x|t) %*% K(t+h|t) + b(x,i|t) %*% k(t+h,i|t) for region i
  new.BKbk = Map(function(x, y) x + Bxh_Kth, bxh_kth)
  log_LL.forecast = Map(function(new.BKbk.i, Axh.i) sweep(new.BKbk.i, 2, Axh.i, FUN = "+"), new.BKbk, Axh)
  ## converged LL forecast
  log_LL.forecast = setNames(log_LL.forecast, region_labels)

  # DLNM part forecasting
  ## (1) dlnm_coef %*% smoothing crossbasis matrix (each iteration each age group)
  cb_matrix.forecast = lapply(seq_along(DLNM_LL_model$dlnm_coef), function(i) { 
    region_cb_test <- cb_matrix_test[[i]] #####
    n_cb = ncol(region_cb_test)
    lapply(DLNM_LL_model$dlnm_coef[[i]], function(iteration) {
      lapply(iteration, function(age_group) {
        as.numeric(age_group[1] + region_cb_test %*% age_group[seq_len(n_cb) + 1])
      })
    })
  })
  cb_matrix.forecast = setNames(cb_matrix.forecast, region_labels) # rename
  ## (1.1) Sum up all forecast values from DLNMs.
  cb_matrix.sum.forecast = lapply(cb_matrix.forecast, function(region) {
    Reduce(function(x, y) {
      mapply(`+`, x, y, SIMPLIFY = FALSE)
    }, region)
  })
  ## (1.2) Combine data to dataframe for each region and store in one list
  cb_matrix.region.sum.forecast = lapply(cb_matrix.sum.forecast, function(region){
    do.call(cbind, region)
  })
  
  ## (2) dlnm_coef %*% hot/cold wave term (each iteration each age group)
  wave.forecast = lapply(seq_along(DLNM_LL_model$dlnm_coef), function(i) { 
    region_wave_test <- get_region_wave(i)
    has_wave = !is.null(region_wave_test) && ncol(as.data.frame(region_wave_test)) > 0
    region_wave_test = if (has_wave) as.matrix(region_wave_test) else NULL
    n_cb = ncol(cb_matrix_test[[i]])

    lapply(DLNM_LL_model$dlnm_coef[[i]], function(iteration) {
      lapply(iteration, function(age_group) {
        if (!has_wave) {
          return(rep(0, forecast_step))
        }

        wave_coef = age_group[-seq_len(n_cb + 1)]
        if (length(wave_coef) == 0) {
          return(rep(0, forecast_step))
        }
        as.numeric(region_wave_test %*% wave_coef)
      })
    })
  })
  wave.forecast = setNames(wave.forecast, region_labels)  # rename by region
  ## (2.1) Sum up all forecast values from DLNMs.
  wave.sum.forecast = lapply(wave.forecast, function(region) {
    Reduce(function(x, y) {
      mapply(`+`, x, y, SIMPLIFY = FALSE)
    }, region)
  })
  ## (2.2) Combine data to dataframe for each region and store in one list
  wave.sum.forecast = lapply(wave.sum.forecast, function(region){
    do.call(cbind, region)
  })
  
  # Combine Lee--Li forecasting part, cross-basis matrix forecasting part and hot/cold wave part
  model.forecast = lapply(names(log_LL.forecast), function(region) {
    log_LL.forecast[[region]] + cb_matrix.region.sum.forecast[[region]] + wave.sum.forecast[[region]]
  })
  names(model.forecast) = names(log_LL.forecast)
  model.forecast = lapply(model.forecast, function(region) {
    region = as.data.frame(region)
    colnames(region) =  c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
    return(region)
  })
  
  # Baseline LC part
  ## Baseline forecasting
  Ax_ref = DLNM_LL_model$baseline_LL$Ax
  Bx_ref = DLNM_LL_model$baseline_LL$Bx
  bx_ref = DLNM_LL_model$baseline_LL$bx
  kth_ref = kt_fitting$base_kt$kt.pred_ref[-1]
  Kth_ref = kt_fitting$base_kt$kt.pred_ref[[1]]
  BxKt_ref = outer(Kth_ref, Bx_ref, FUN = "*")
  bxkt_ref = mapply(function(kth.i, bxh.i) outer(kth.i, bxh.i, FUN = "*"), kth_ref, bx_ref,  SIMPLIFY = FALSE)
  new.BKbk_ref = Map(function(x, y) x + BxKt_ref, bxkt_ref)
  names(new.BKbk_ref) = names(DLNM_LL_model$log_fitted)
  log_LL_baseline.forecast = Map(function(region_df, ax_values) {
    sweep(region_df, 2, ax_values, FUN = "+")
  }, new.BKbk_ref, Ax_ref)
  log_LL_baseline.forecast = lapply(log_LL_baseline.forecast, function(region) {
    region = as.data.frame(region)
    colnames(region) = c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
    return(region)
  })
  
  # Output
  output = list(baseline_LL.parameter = list(Ax = Ax_ref, Bx = Bx_ref, bx = bx_ref, 
                                             Kt = Kth_ref, kt = kth_ref),
                baseline_LL.forecast = log_LL_baseline.forecast,
                DLNM_LL_final.parameter = list(Ax = Axh, Bx = Bxh, bx = bxh, 
                                                 Kt = Kth, kt = kth),
                DLNM_LL_final.forecast = list(model.forecast = model.forecast,
                                                converged_LL.forecast = log_LL.forecast,
                                                # need to change
                                                cb_matrix.forecast = cb_matrix.forecast,
                                                cb_matrix.sum.forecast = cb_matrix.sum.forecast,
                                                cb_matrix.region.sum.forecast = cb_matrix.region.sum.forecast,
                                                wave.forecast = wave.forecast,
                                                wave.sum.forecast = wave.sum.forecast))
  return(output)
}
