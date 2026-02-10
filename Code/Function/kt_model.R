#' kt_model: Time-series model for time-varying factors
#'
#' @param kt Numeric vector of period index values (k_t).
#' @param method Character string specifying the model: \code{"RWD"} or \code{"SARIMA"}.
#'
#' @return A fitted ARIMA model object (from \code{stats::arima} or \code{forecast::auto.arima}),
#'   stored as \code{ts.fit}.

kt_model = function(kt, method){
  if (method == "RWD"){
    ts.fit = arima(ts(kt),  order = c(0,1,0))
  }
  if (method == "SARIMA"){
    kt = ts(kt, frequency = 52)
    ts.fit = auto.arima(kt, seasonal = T, max.d = 0, allowmean = F)
  }
  return(ts.fit)
}
