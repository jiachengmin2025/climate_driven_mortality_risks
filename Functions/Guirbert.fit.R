Guibert.fit <- function(regions_data, cb_list,
                        age_groups = c("Y20_64","Y65_74","Y75_84","Y_GE85"),
                        ns_df = 35, time_var = "time") {
  
  ## Step 1 & 2: DLNM estimation and non-temperature attributable mortality rates
  dlnm_estimate <- function(dat, age_groups = c("Y20_64","Y65_74","Y75_84","Y_GE85"),
                            cb_matrix, ns_df = 35, time_var = "time") {
    # DLNM fitting
    fits <- setNames(lapply(age_groups, function(age) {
      cb <- cb_matrix
      glm(y ~ cb + ns(time, df = ns_df), family = quasipoisson(),
          data = list(y = dat[[age]], time = dat[[time_var]], cb = cb))
    }), age_groups)
    
    # Temperature Adjustment Factor (from fitted cross-basis matrix)
    temp_adjustment <- function(model) {
      X <- model.matrix(model)
      b <- coef(model)
      keep <- !is.na(b) & names(b) %in% colnames(X)
      b <- b[keep]; X <- X[, names(b), drop = FALSE]
      idx <- startsWith(names(b), "cb")        
      cb_fitted <- drop(X[, idx, drop = FALSE] %*% b[idx])
      exp(cb_fitted)
    }
    temp_adj.df <- do.call(cbind, lapply(fits, temp_adjustment))
    temp_adj.df <- data.frame(temp_adj.df)
    colnames(temp_adj.df) <- age_groups
    rownames(temp_adj.df) <- rownames(dat)
    
    # Non-temperature attributable mortality rates
    non_temp <- data.frame(dat[, age_groups, drop = FALSE]) / temp_adj.df
    rownames(non_temp) <- rownames(dat)
    
    # Output
    list(fits = fits, temp_adj_factor.df = temp_adj.df, non_temp_mort = non_temp)
  }
  
  ## Run example (internal call per region)
  dlnm_estimates_list <- setNames(lapply(names(regions_data), function(r) {
    dlnm_estimate(dat = regions_data[[r]], 
                  age_groups = age_groups, 
                  cb_matrix = cb_list[[r]], ns_df = ns_df, time_var = time_var)
  }), names(regions_data))
  
  ## Step 3: Fit Li--Lee model on non-temperature attributable mortality rates
  non.temp.data_list <- setNames(lapply(names(regions_data), function(r) {
    dlnm_estimates_list[[r]]$non_temp_mort
  }), names(regions_data))
  LL_fit <- LL_model(non.temp.data_list)
  
  ## Step 4: Combine fitted values from two sub-models
  Guibert_fitted <- function(LL_fit, dlnm_estimates_list, regions = NULL, age_groups = NULL) {
    if (is.null(regions)) {
      regions <- intersect(names(LL_fit$log_fitted), names(dlnm_estimates_list))
    }
    if (length(regions) == 0) stop("No overlapping regions between LL_fit and dlnm_estimates_list.")
    if (is.null(age_groups)) {
      age_groups <- colnames(LL_fit$log_fitted[[regions[1]]])
    }
    res <- setNames(lapply(regions, function(r) {
      LL.log.fitted <- LL_fit$log_fitted[[r]]                      
      temp_adj_factor <- dlnm_estimates_list[[r]]$temp_adj_factor.df  
      
      # combine fitted mortality rates from two sub-models
      mort.fitted <- sapply(age_groups, function(age) {
        m0_hat <- exp(LL.log.fitted[, age])  
        m0_hat * temp_adj_factor[, age] 
      })
      
      colnames(mort.fitted) <- age_groups
      rownames(mort.fitted) <- rownames(LL.log.fitted)
      data.frame(mort.fitted)
    }), regions)
    
    return(res)
  }
  
  dlnm_estimates_list <- dlnm_estimates_list
  LL_fit <- LL_fit
  Guibert_fitted.res <- Guibert_fitted(LL_fit, dlnm_estimates_list)
  
  list(
    dlnm_estimates_list = dlnm_estimates_list,
    non.temp.data_list  = non.temp.data_list,
    LL_fit = LL_fit,
    Guibert_fitted.res = Guibert_fitted.res
  )
}