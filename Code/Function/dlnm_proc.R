#' dlnm_proc: Fit Gaussian log-scale DLNM across different age groups.
#'
#' @param dat Log-mortality residual surface from the current Lee--Carter step.
#'   Columns correspond to age groups: 20--64, 65--74, 75--84, 85+.
#' @param cb_train Cross-basis matrix.
#' @param wave_data A list or data.frame containing wave covariates. If supplied
#'   with standard names, \code{<region>_hot_wave3} and
#'   \code{<region>_cold_wave3} are used. If \code{NULL} or empty, the DLNM is
#'   fitted without wave terms.
#' @param region Character region label used to locate wave variables in
#'   \code{wave_data}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{DLNM_fitted}: fitted DLNM components on the log-rate scale.
#'   \item \code{DLNM20}, \code{DLNM65}, \code{DLNM75}, \code{DLNM85}: fitted GLM objects.
#'   \item \code{DLNM_coef}: list of coefficient vectors by age group.
#'   \item \code{DLNM_vcov}: list of variance--covariance matrices by age group.
#' }

dlnm_proc = function(dat, cb_train, wave_data = NULL, region = NULL){
  wave_df = NULL
  has_wave = !is.null(wave_data) && ncol(as.data.frame(wave_data)) > 0

  if (has_wave) {
    wave_df = as.data.frame(wave_data)
    if (!is.null(region)) {
      hot_wave_var = paste0(region, "_hot_wave3")
      cold_wave_var = paste0(region, "_cold_wave3")
      standard_wave_vars = c(hot_wave_var, cold_wave_var)

      if (all(standard_wave_vars %in% names(wave_df))) {
        wave_df = wave_df[, standard_wave_vars, drop = FALSE]
      }
    }
  }

#' dlnm_proc internal: fit one age-specific Gaussian DLNM
#'
#' @param y Numeric vector of log-rate residuals for one age group.
#'
#' @return A fitted Gaussian \code{glm} object with cross-basis terms and,
#'   when available, hot/cold wave covariates.
  fit_dlnm = function(y) {
    if (has_wave) {
      glm(y ~ cb_train + ., data = wave_df, family = gaussian())
    } else {
      glm(y ~ cb_train, family = gaussian())
    }
  }
  
  # 20-64
  DLNM20 = fit_dlnm(dat[,1])
  DLNM20_fitted = fitted.values(DLNM20)
  
  # 65-74
  DLNM65 = fit_dlnm(dat[,2])
  DLNM65_fitted = fitted.values(DLNM65)
  
  # 75-84
  DLNM75 = fit_dlnm(dat[,3])
  DLNM75_fitted = fitted.values(DLNM75)
  
  # 85+
  DLNM85 = fit_dlnm(dat[,4])
  DLNM85_fitted = fitted.values(DLNM85)
  
  DLNM_result = list(DLNM_fitted = data.frame(DLNM20_fitted, DLNM65_fitted, DLNM75_fitted, DLNM85_fitted),
                     DLNM20 = DLNM20, DLNM65 = DLNM65, DLNM75 = DLNM75, DLNM85 = DLNM85,
                     DLNM_coef = list(Y20_64 = coef(DLNM20), Y65_74 = coef(DLNM65), 
                                      Y75_84 = coef(DLNM75), Y_GE85 = coef(DLNM85)),
                     DLNM_vcov = list(Y20_64 = vcov(DLNM20), Y65_74 = vcov(DLNM65), 
                                      Y75_84 = vcov(DLNM75), Y_GE85 = vcov(DLNM85)))
  
  return(DLNM_result)
}
