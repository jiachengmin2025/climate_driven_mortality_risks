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
