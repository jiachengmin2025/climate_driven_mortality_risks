#' LC_model: function for fitting Lee--Carter model via SVD.
#'
#' @param df Mortality surface (matrix/data.frame). Rows/columns represent age/time
#'   depending on your data layout.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{df}: mortality data (age * time).
#'   \item \code{a_x}: mean log-mortality by the margin used for centring.
#'   \item \code{b_x}: age-specific coefficient.
#'   \item \code{k_t}: time-varying factors.
#'   \item \code{std_mxt}: standardised log-mortality \eqn{\log m_{x,t} - a_x} (includes residual).
#'   \item \code{log_fitted}: fitted \eqn{a_x + b_x k_t} on the log scale.
#'   \item \code{log_matrix}: observed \eqn{\log m_{x,t}}.
#'   \item \code{log_bk_fitted}: fitted \eqn{b_x k_t} component on the log scale.
#' }
#'
#' @details
#' The normalisation used here sets \code{b_x = v_1 / sum(v_1)} and
#' \code{k_t = d_1 u_1 * sum(v_1)}, so that \code{b_x} sums to 1.
#' This is one of several equivalent LC identifiability conventions.


LC_model = function(df){
  log_matrix = log(as.matrix(df))   
  # a_x
  a_x = apply(log_matrix, 2, mean, na.rm = T)
  log_sd_matrix = sweep(log_matrix, 2, a_x, FUN = "-")
  # SVD
  U = svd(log_sd_matrix)$u
  V = svd(log_sd_matrix)$v
  d = svd(log_sd_matrix)$d
  # k_t
  kappa_t = d[1]*U[,1]*sum(V[,1])
  # fitted standard matrix bx %*% kt (log scale)
  log_sd_matrix_fitted = as.matrix(kappa_t) %*% t(V[,1]/sum(V[,1])) 
  # fitted matrix ax + bx %*% kt (log scale)
  log_matrix_fitted = sweep(log_sd_matrix_fitted, 2, a_x, FUN = "+")
  # b_x
  b_x = as.numeric(t(V[,1]/sum(V[,1])))
  # residual result log(mxt) - ax - bx %*% kt (log scale)
  LC_residual = data.frame(log_matrix - log_matrix_fitted)
  # std_mxt = bx %*% kt + eps
  std_mxt = log_sd_matrix_fitted + LC_residual
  result = list(df=df, a_x=a_x, b_x=b_x, k_t=kappa_t, std_mxt=std_mxt,
                log_fitted = log_matrix_fitted, log_matrix = log_matrix,
                log_bk_fitted = log_sd_matrix_fitted)
  return(result)
}