#' B-DLNM-LC: Backfitting DLNM--Lee--Carter model
#' Functions for fitting DLNM--LC model.
#'
#' @param dat Mortality input passed to \code{LC_model()}.
#' @param cross_matrix Cross-basis matrix used by \code{dlnm_proc()}.
#' @param wave_data Wave indicators used by \code{dlnm_proc()}.
#' @param tol Convergence tolerance for LC parameters (\code{a_x}, \code{b_x}, \code{k_t}).
#' @param max_iter Maximum number of back-fitting iterations.
#' @param region Region label passed to \code{dlnm_proc()}.
#'
#' @return A list:
#' \itemize{
#'   \item \code{final_LC}: the converged LC model.
#'   \item \code{log_fitted}: overall fitted log-mortality rates from the backfitting algorithm.
#'   \item \code{baseline_LC}: the initial Lee--Carter fit before adding the DLNM component.
#'   \item \code{dlnm_coef}: DLNM coefficient vectors for each iteration.
#'   \item \code{dlnm_vcov}: DLNM variance--covariance matrices for each iteration.
#'   \item \code{iteration}: number of iterations performed.
#' }

B_DLNM_LC = function(dat, cross_matrix, wave_data, tol, max_iter, region){
  # Initial set up
  iter = 0
  log_fitted = matrix(data = 0, nrow = nrow(dat), ncol = ncol(dat))
  dlnm_coef = list()
  dlnm_vcov = list()
  baseline_LC = LC_model(dat) # Original LC
  result0 = LC_model(dat)
  # DLNM iteration
  while (TRUE) {
    iter = iter + 1
    dlnm_result = dlnm_proc(result0$std_m, cross_matrix, wave_data, region)
    ## Store DLNM coefficients each time
    dlnm_coef = append(dlnm_coef, list(dlnm_result$DLNM_coef))
    dlnm_vcov = append(dlnm_vcov, list(dlnm_result$DLNM_vcov))
    names(dlnm_coef)[iter] = paste("DLNM_iteration_", iter, sep = "")
    names(dlnm_vcov)[iter] = paste("DLNM_iteration_", iter, sep = "")
    update.df = log(result0$df) - dlnm_result$DLNM_fitted
    result = LC_model(exp(update.df))
    if (all(abs(result$a_x - result0$a_x) < tol) && 
        all(abs(result$b_x - result0$b_x) < tol) &&
        all(abs(result$k_t - result0$k_t) < tol)) {
      log_fitted = log_fitted + result0$log_fitted
      cat("converged\n")
      break
    }
    if (iter >= max_iter) {  
      log_fitted = log_fitted + result0$log_fitted
      cat("reach maximum iterations/not converged")
      break
    }
    # Add DLNM fitted (r-1) part
    log_fitted = log_fitted + dlnm_result$DLNM_fitted
    result0 = result
  }
  return(list(final_LC = result0, log_fitted = log_fitted, 
              baseline_LC = baseline_LC, dlnm_coef = dlnm_coef,
              dlnm_vcov = dlnm_vcov, iteration = iter))
}
