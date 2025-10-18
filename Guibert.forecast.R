Guibert.forecast <- function(LL_fit, dlnm_estimates_list, cb_matrix_test, forecast_step) {
  # Forecast common K_t with RWD
  K_fit <- auto.arima(ts(LL_fit$Kt, frequency = 52), seasonal = T, max.d = 0, allowmean = F)
  K_forecast  <- as.numeric(forecast(K_fit, h = forecast_step)$mean)
  
  # Forecast regional kappa_t with AR(1) 
  regions  <- names(LL_fit$kt)
  kappa_forecast <- setNames(lapply(regions, function(r){
    fit_r <- auto.arima(ts(LL_fit$kt[[r]], frequency = 52), seasonal = T, max.d = 0, allowmean = F)
    as.numeric(forecast(fit_r, h = forecast_step)$mean)
  }), regions)
  
  # Compute log baseline forecasts: A_x + B_x * K_t + b_x * kappa_{t,i}
  age_groups <- colnames(LL_fit$log_fitted[[regions[1]]])
  rown <- paste0("t+", seq_len(forecast_step))
  
  # Keep unchanged parameters
  Bxh <- setNames(as.numeric(LL_fit$Bx)[seq_along(age_groups)], age_groups)
  Axh <- lapply(LL_fit$Ax, function(v) setNames(as.numeric(v)[seq_along(age_groups)], age_groups))
  bxh <- lapply(LL_fit$bx, function(v) setNames(as.numeric(v)[seq_along(age_groups)], age_groups))
  
  ## B(x|t) %*% K(t+h|t)
  Bxh_Kth <- outer(K_forecast, Bxh, FUN = "*")
  
  ## b(x,i|t) %*% k(t+h,i|t)
  bxh_kth <- mapply(function(kth.i, bxh.i) {
    outer(kth.i, bxh.i, FUN = "*")
  }, kappa_forecast[regions], bxh[regions], SIMPLIFY = FALSE)
  
  ## B(x|t) %*% K(t+h|t) + b(x,i|t) %*% k(t+h,i|t) for region i
  new.BKbk <- lapply(bxh_kth, function(x) x + Bxh_Kth)
  
  ## Log LL forecast A(x|t) + B(x|t) %*% K(t+h|t) + b(x,i|t) %*% k(t+h,i|t)
  log_LL.forecast <- Map(function(new.BKbk.i, Axh.i) {
    out <- sweep(new.BKbk.i, 2, Axh.i, FUN = "+")
    colnames(out) <- age_groups
    rownames(out) <- rown
    out
  }, new.BKbk, Axh[regions])
  
  # Name by regions
  log_LL.forecast <- setNames(log_LL.forecast, regions)
  
  # DLNM forecast block (T and eta_cb on test cross-basis)
  regions_d <- names(dlnm_estimates_list)
  age_groups_d <- names(dlnm_estimates_list[[regions_d[1]]]$fits)
  
  T_list <- setNames(lapply(regions_d, function(r) {
    Xnew <- as.matrix(cb_matrix_test[[r]])
    rn <- rownames(Xnew)
    out <- sapply(age_groups_d, function(age) {
      fit <- dlnm_estimates_list[[r]]$fits[[age]]
      b <- coef(fit); b <- b[!is.na(b)]
      idx <- startsWith(names(b), "cb")
      b_cb <- b[idx]
      Xcb <- Xnew
      colnames(Xcb) <- names(b_cb)
      as.numeric(exp(-drop(Xcb %*% b_cb)))
    })
    colnames(out) <- age_groups_d
    rownames(out) <- if (is.null(rn)) seq_len(nrow(out)) else rn
    out
  }), regions_d)
  
  eta_list <- setNames(lapply(regions_d, function(r) {
    Xnew <- as.matrix(cb_matrix_test[[r]])
    rn <- rownames(Xnew)
    out <- sapply(age_groups_d, function(age) {
      fit <- dlnm_estimates_list[[r]]$fits[[age]]
      b <- coef(fit); b <- b[!is.na(b)]
      idx <- startsWith(names(b), "cb")
      b_cb <- b[idx]
      Xcb <- Xnew
      colnames(Xcb) <- names(b_cb)
      as.numeric(-drop(Xcb %*% b_cb))
    })
    colnames(out) <- age_groups_d
    rownames(out) <- if (is.null(rn)) seq_len(nrow(out)) else rn
    out
  }), regions_d)
  
  # m_hat = exp(log_base) * T
  regions_c <- intersect(names(log_LL.forecast), names(T_list))
  res <- setNames(lapply(regions_c, function(r) {
    log_base <- as.matrix(log_LL.forecast[[r]])[, age_groups, drop = FALSE]
    T_mat <- as.matrix(T_list[[r]])[, age_groups, drop = FALSE]
    out <- exp(log_base) / T_mat
    colnames(out) <- age_groups
    rownames(out) <- rownames(log_base)
    as.data.frame(out, check.names = FALSE)
  }), regions_c)
  
  # Output
  list(
    final_kt = list(kt.pred = c(list(K = K_forecast), kappa_forecast)),
    log_LL.forecast.list = log_LL.forecast,
    temp_adj_factor.forecast = T_list,
    cb.forecast = eta_list,
    Guibert.forecast.list = res
  )
}