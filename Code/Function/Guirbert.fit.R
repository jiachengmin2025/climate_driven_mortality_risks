#' Guibert comparator: fit DLNM adjustment and Li--Lee model
#' Estimates age-specific DLNM temperature effects, constructs
#' non-temperature-attributable mortality rates, and fits a Li--Lee model to
#' the adjusted rates.
#'
#' @param regions_data Named list of region-specific mortality-rate data.frames.
#' @param cb_list Named list of UTCI cross-basis matrices aligned with
#'   \code{regions_data}.
#' @param age_groups Character vector of age-group column names.
#' @param ns_df Degrees of freedom for the time spline used in the DLNM step.
#' @param time_var Column name containing the time index.
#'
#' @return A list:
#' \itemize{
#'   \item \code{dlnm_estimates_list}: fitted DLNM models and temperature factors.
#'   \item \code{non.temp.data_list}: non-temperature-attributable mortality rates.
#'   \item \code{LL_fit}: Li--Lee fit to the adjusted rates.
#'   \item \code{Guibert_fitted.res}: fitted mortality rates after recombining both components.
#' }
Guibert.fit <- function(regions_data, cb_list,
                        age_groups = c("Y20_64","Y65_74","Y75_84","Y_GE85"),
                        ns_df = 35, time_var = "time") {
  
  ## Step 1 & 2: DLNM estimation and non-temperature attributable mortality rates
#' Guibert internal: estimate age-specific DLNM effects
#'
#' @param dat Region-specific data.frame containing mortality-rate age columns
#'   and a time column.
#' @param age_groups Character vector of age-group column names.
#' @param cb_matrix UTCI cross-basis matrix for this region.
#' @param ns_df Degrees of freedom for the time spline.
#' @param time_var Column name containing the time index.
#'
#' @return A list containing fitted GLMs, temperature-adjustment factors, and
#'   non-temperature-attributable mortality rates.
  dlnm_estimate <- function(dat, age_groups = c("Y20_64","Y65_74","Y75_84","Y_GE85"),
                            cb_matrix, ns_df = 35, time_var = "time") {
    # DLNM fitting
    fits <- setNames(lapply(age_groups, function(age) {
      cb <- cb_matrix
      glm(y ~ cb + ns(time, df = ns_df), family = quasipoisson(),
          data = list(y = dat[[age]], time = dat[[time_var]], cb = cb))
    }), age_groups)
    
    # Temperature Adjustment Factor (from fitted cross-basis matrix)
#' Guibert internal: compute temperature-adjustment factor
#'
#' @param model Fitted age-specific quasipoisson GLM from the Guibert DLNM step.
#'
#' @return Numeric vector \code{exp(eta_cb)} containing the fitted
#'   temperature-attributable multiplicative factor.
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
#' Guibert internal: recombine Li--Lee and DLNM fitted values
#'
#' @param LL_fit Fitted Li--Lee model on non-temperature-attributable rates.
#' @param dlnm_estimates_list Region-specific DLNM estimates.
#' @param regions Optional character vector of regions to recombine.
#' @param age_groups Optional character vector of age-group column names.
#'
#' @return Named list of fitted mortality-rate data.frames by region.
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
