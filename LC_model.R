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