#' B-DLNM-LL: Function for fitting Back-fitting DLNM--Li--Lee model.
#'
#' @param dat_list Named list of mortality surfaces across regions.
#' @param cb_list Named list of cross-basis objects passed to \code{dlnm_proc()},
#'   aligned with \code{dat_list}.
#' @param wave_list Named list of wave covariates passed to \code{dlnm_proc()},
#'   aligned with \code{dat_list}.
#' @param tol Convergence tolerance for Li--Lee parameters (\code{Ax}, \code{Bx}, \code{Kt}, \code{bx}, \code{kt}).
#' @param max_iter Maximum number of back-fitting iterations.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{final_LL}: the final (converged) Li--Lee model, or the last fit if not converged.
#'   \item \code{baseline_LL}: baseline Li--Lee model.
#'   \item \code{log_fitted}: overall fitted log-mortality surfaces for each region (list).
#'   \item \code{dlnm_fitted}: fitted DLNM components.
#'   \item \code{dlnm_coef}: DLNM coefficient lists by region and iteration (age-group specific).
#'   \item \code{dlnm_vcov}: DLNM variance--covariance lists by region and iteration (age-group specific).
#'   \item \code{iteration}: number of iterations performed.
#' }


B_DLNM_LL = function(dat_list, cb_list, wave_list, tol, max_iter){
  # Li-Lee model
  result0 = LL_model(dat_list)
  # Baseline Li-Lee model
  baseline_LL = result0
  names(baseline_LL$log_fitted) = names(dat_list)
  baseline_LL$log_fitted = lapply(baseline_LL$log_fitted, function(region) {
    df = data.frame(region)
    colnames(df) = c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
    return(df)
  })
  
  # Initialize empty matrix for log_fitted, dlnm_coef, and dlnm_vcov
  iter = 0
  all_regions = names(dat_list) # record region names
  log_fitted_list = lapply(dat_list, function(region_data) {
    matrix(data = 0, nrow = nrow(region_data), ncol = ncol(region_data))
  })
  names(log_fitted_list) = names(dat_list)
  # DLNM coefficients
  dlnm_coef_list = lapply(all_regions, function(region) list())
  names(dlnm_coef_list) = all_regions
  # DLNM coefficient variance-covariance matrix
  dlnm_vcov_list = lapply(all_regions, function(region) list())
  names(dlnm_vcov_list) = all_regions
  # DLNM accumulation
  dlnm_fitted = lapply(all_regions, function(region) list())
  names(dlnm_fitted) = all_regions
  
  # DLNM iteration
  while (TRUE) {
    iter = iter + 1
    ## Fit DLNM on B(x) %*% K(t) + b(x) %*% k(t) + eps(x,t,i) (under log scale)
    # Fit DLNM for common components
    # Common BK
    dlnm_results_list = lapply(all_regions, function(region) {
      dlnm_proc(result0$final_std_mxt[[region]], cb_list[[region]], wave_list[[region]], 
                region = region)
    })
    names(dlnm_results_list) = names(dat_list)
    ## Extract coefficients and variance-covariance matrix
    for (i in seq_along(dlnm_results_list)) {
      region = all_regions[i]
      ## coefficients
      dlnm_coef_list[[region]] = append(dlnm_coef_list[[region]], list(dlnm_results_list[[i]]$DLNM_coef))
      names(dlnm_coef_list[[region]])[iter] = paste("DLNM_iteration_", iter, sep = "")
      ## variance-covariance matrix
      dlnm_vcov_list[[region]] = append(dlnm_vcov_list[[region]], list(dlnm_results_list[[i]]$DLNM_vcov))
      names(dlnm_vcov_list[[region]])[iter] = paste("DLNM_iteration_", iter, sep = "")
      dlnm_fitted[[region]] = append(dlnm_fitted[[region]], list(dlnm_results_list[[i]]$DLNM_fitted))
      names(dlnm_fitted[[region]])[iter] = paste("DLNM_iteration_", iter, sep = "")
    }
    ## Update LL_model
    ## Update data: log(m(x,t,i;j)) = log(m(x,t,i)) - B(x;j) %*% K(t;j) - b(x,i;j) %*% k(t,i;j)
    update_list = lapply(names(dat_list), function(region) {
      log(result0$dat_list[[region]]) - dlnm_results_list[[region]]$DLNM_fitted
    })
    names(update_list) = names(dat_list)
    result = LL_model(lapply(update_list, exp))
    
    # Converging condition
    if (all(abs(unlist(result$Ax) - unlist(result0$Ax)) < tol) &&
        all(abs(unlist(result$Bx) - unlist(result0$Bx)) < tol) &&
        all(abs(unlist(result$Kt) - unlist(result0$Kt)) < tol) &&
        all(abs(unlist(result$bx) - unlist(result0$bx)) < tol) &&
        all(abs(unlist(result$kt) - unlist(result0$kt)) < tol)) {
      
      log_fitted_list = mapply(function(log_fitted, region_fitted) {
        return(log_fitted + region_fitted)
      }, log_fitted_list, result0$log_fitted, SIMPLIFY = FALSE)
      log_fitted_list = lapply(log_fitted_list, function(region) {
        names(region) = c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
        return(region)
      })
      cat("converged\n")
      break
    }
    
    if (iter >= max_iter) {
      log_fitted_list = mapply(function(log_fitted, region_fitted) {
        return(log_fitted + region_fitted)
      }, log_fitted_list, result0$log_fitted, SIMPLIFY = FALSE)
      ## rename
      log_fitted_list = lapply(log_fitted_list, function(region) {
        names(region) = c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
        return(region)
      })
      cat("Reach maximum iterations/not converged\n")
      break
    }

    # Add DLNM fitted (r-1) part
    log_fitted_list = mapply(function(log_fitted, region) {
      # region DLNM fitted
      dlnm_region_fitted = dlnm_results_list[[region]]$DLNM_fitted
      final_log_fitted = log_fitted + dlnm_region_fitted
      # log_fitted + dlnm_region_fitted
      return(final_log_fitted)
    }, log_fitted_list, names(log_fitted_list), SIMPLIFY = FALSE)
    result0 = result
  }
  output = list(final_LL = result0, log_fitted = log_fitted_list, baseline_LL = baseline_LL,
                dlnm_fitted = dlnm_fitted, dlnm_coef = dlnm_coef_list, dlnm_vcov = dlnm_vcov_list,
                iteration = iter)
  return(output)
}