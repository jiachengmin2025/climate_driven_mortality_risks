LL_model = function(dat_list){
  num_regions = length(dat_list)
  region_names = names(dat_list)
  # Product matrix for all regions
  product_matrix = Reduce(function(x, y) x * y, dat_list)^(1 / num_regions)
  # Fit Lee-Carter model on product matrix
  LC_product = LC_model(product_matrix)
  
  # Common factor B(x), K(t) (under log scale)
  Bx = LC_product$b_x
  Kt = LC_product$k_t
  
  # Common B(x) %*% K(t) + E(x,t)
  std_mxt = LC_product$std_mxt
  
  # Fit Lee-Carter model on ratio matrix
  results_list = lapply(dat_list, function(df) {
    ratio = df / product_matrix
    LC_ratio = LC_model(ratio)
    
    ## Residuals for each region (under log scale)
    ## res(t,x,i) = log(f(t,x,i)) - log(p(t,x)) - log(r(t,x,i))
    df_res = log(df) - LC_product$log_fitted - LC_ratio$log_fitted
    
    ## log_fitted = A(x) + B(x) %*% K(t) + b(t,i) %*% k(t,i)
    log_fitted = data.frame(LC_product$log_fitted + LC_ratio$log_fitted)
    colnames(log_fitted) <- colnames(df)
    
    ## Common-specific-mixed factor A(x), region-specific factor b(x,i), k(t,i) (under log scale)
    common_Ax = LC_product$a_x
    individual_ax = LC_ratio$a_x
    Ax = LC_product$a_x + LC_ratio$a_x
    bx = LC_ratio$b_x
    kt = LC_ratio$k_t
    
    
    ## std_mxti for each region
    ## std_mxti = b(x) %*$ k(t) + eps(x,t,i) (under log scale)
    std_mxti = LC_ratio$std_mxt
    
    return(list(df_res = df_res, log_fitted = log_fitted, Ax = Ax, bx = bx, kt = kt, 
                common_Ax = common_Ax, individual_ax = individual_ax, std_mxti = std_mxti))
  })
  
  # Obtain results
  df_res_list = lapply(results_list, function(x) x$df_res)
  log_fitted_list = lapply(results_list, function(x) x$log_fitted)
  common_Ax = lapply(results_list, function(x) x$common_Ax)
  individual_ax = lapply(results_list, function(x) x$individual_ax)
  Ax_list = lapply(results_list, function(x) x$Ax)
  bx_list = lapply(results_list, function(x) x$bx)
  kt_list = lapply(results_list, function(x) x$kt)
  std_mxti_list = lapply(results_list, function(x) x$std_mxti)
  final_std_mxt = lapply(std_mxti_list, function(df) df + std_mxt)
  
  names(df_res_list) <- region_names
  names(log_fitted_list) <- region_names
  names(Ax_list) <- region_names
  names(bx_list) <- region_names
  names(kt_list) <- region_names
  
  
  result = list(dat_list = dat_list, Ax = Ax_list, Bx = Bx, Kt = Kt, bx = bx_list, 
                kt = kt_list, common_Ax = common_Ax, individual_ax = individual_ax,
                std_mxt = std_mxt, std_mxti = std_mxti_list, final_std_mxt = final_std_mxt,
                log_fitted = log_fitted_list)
  return(result)
}