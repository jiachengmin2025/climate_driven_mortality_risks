#' DLNM residual diagnostics and seasonality tests
#' Fits DLNM--LC and DLNM--LL models and computes residual seasonality tests
#' for Supplementary Section D.
#'
#' @param base_dir Project root containing \code{Data}, \code{Code}, and
#'   \code{Results} folders.
#' @param residual_scale Character string. If \code{"log_rate"}, residuals are
#'   observed log rates minus fitted log rates; otherwise residuals are computed
#'   on the mortality-rate scale.
#' @param lag_days Maximum lag in days for the UTCI cross-basis.
#' @param var_df Degrees of freedom for the UTCI exposure dimension.
#' @param lag_df Degrees of freedom for the lag dimension.
#' @param degree Degree passed to \code{dlnm::equalknots()}.
#' @param tol Convergence tolerance passed to \code{DLNM_LC()} and
#'   \code{DLNM_LL()}.
#'
#' @return A list with \code{DLNM_LC} and \code{DLNM_LL} entries. Each contains
#'   fitted models, residuals, seasonality-test p-values, and display names.
DLNM_residual_diagnostics <- function(base_dir,
                                      residual_scale = "log_rate",
                                      lag_days = 21,
                                      var_df = 3,
                                      lag_df = 3,
                                      degree = 3,
                                      tol = 1e-2) {
  data_dir <- file.path(base_dir, "Data/Combined_data")

  Y20_64 <- readxl::read_xlsx(file.path(data_dir, "Y20_64_combined.xlsx"))
  Y65_74 <- readxl::read_xlsx(file.path(data_dir, "Y65_74_combined.xlsx"))
  Y75_84 <- readxl::read_xlsx(file.path(data_dir, "Y75_84_combined.xlsx"))
  Y_GE85 <- readxl::read_xlsx(file.path(data_dir, "Y_GE85_combined.xlsx"))

  age_groups <- c("20-64", "65-74", "75-84", "85+")
  rate_colnames <- c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
  region_config <- list(
    Attiki = list(display_region = "Athens", death_col = 3,
                  exposure_col = 9, max_iter = 20),
    Lisbon = list(display_region = "Lisbon", death_col = 4,
                  exposure_col = 10, max_iter = 20),
    Roma = list(display_region = "Rome", death_col = 6,
                exposure_col = 12, max_iter = 15)
  )

#' DLNM diagnostics internal: build annualized mortality rates
#'
#' @param cfg Region configuration list containing death/exposure columns.
#'
#' @return Data.frame of annualized mortality rates by age group.
  make_rate <- function(cfg) {
    rate <- cbind(
      Y20_64[[cfg$death_col]] / Y20_64[[cfg$exposure_col]],
      Y65_74[[cfg$death_col]] / Y65_74[[cfg$exposure_col]],
      Y75_84[[cfg$death_col]] / Y75_84[[cfg$exposure_col]],
      Y_GE85[[cfg$death_col]] / Y_GE85[[cfg$exposure_col]]
    )
    colnames(rate) <- rate_colnames
    rownames(rate) <- Y20_64$Week
    rate * 52
  }

#' DLNM diagnostics internal: build hot/cold wave covariates
#'
#' @param region Legacy region label used in the combined data files.
#'
#' @return Data.frame with region-specific hot-wave and cold-wave indicators.
  make_wave <- function(region) {
    wave <- data.frame(
      hot_wave = Y20_64[[paste0(region, "_hot_wave3")]],
      cold_wave = Y20_64[[paste0(region, "_cold_wave3")]]
    )
    colnames(wave) <- c(paste0(region, "_hot_wave3"),
                        paste0(region, "_cold_wave3"))
    wave
  }

#' DLNM diagnostics internal: load historical UTCI lag matrix
#'
#' @param region Legacy region label.
#'
#' @return Matrix of historical weekly UTCI lag values for the region.
  load_region_utci <- function(region) {
    load_crossbasis_history(
      base_dir,
      display_region = region_config[[region]]$display_region,
      legacy_region = region
    )
  }

#' DLNM diagnostics internal: create LC cross-basis object
#'
#' @param utci Historical UTCI lag matrix.
#'
#' @return A \code{dlnm::crossbasis} object for a single-region DLNM--LC fit.
  make_lc_crossbasis <- function(utci) {
    varknots <- dlnm::equalknots(utci, fun = "ns", df = var_df,
                                 degree = degree)
    lagknots <- dlnm::logknots(lag_days, lag_df)
    dlnm::crossbasis(utci, lag = lag_days,
                     argvar = list(fun = "ns", knots = varknots),
                     arglag = list(knots = lagknots, df = lag_df))
  }

#' DLNM diagnostics internal: create LL cross-basis object
#'
#' @param utci Historical UTCI lag matrix.
#'
#' @return A \code{dlnm::crossbasis} object for a multi-region DLNM--LL fit.
  make_ll_crossbasis <- function(utci) {
    varknots <- dlnm::equalknots(utci, fun = "ns", df = var_df,
                                 degree = degree)
    lagknots <- dlnm::logknots(lag_days, lag_df)
    dlnm::crossbasis(utci, lag = lag_days,
                     argvar = list(fun = "ns", knots = varknots,
                                   df = var_df),
                     arglag = list(knots = lagknots, df = lag_df))
  }

#' DLNM diagnostics internal: calculate model residuals
#'
#' @param rate Observed annualized mortality-rate matrix/data.frame.
#' @param log_fitted Fitted log-mortality matrix/data.frame from a DLNM model.
#'
#' @return Data.frame of residuals with standard age-group column names.
  calculate_residual <- function(rate, log_fitted) {
    rate_matrix <- as.matrix(rate)
    log_fitted_matrix <- as.matrix(log_fitted)

    residual <- if (residual_scale == "log_rate") {
      log(rate_matrix) - log_fitted_matrix
    } else {
      rate_matrix - exp(log_fitted_matrix)
    }
    residual[!is.finite(residual)] <- NA_real_

    residual <- as.data.frame(residual)
    colnames(residual) <- rate_colnames
    residual
  }

#' DLNM diagnostics internal: run residual seasonality tests
#'
#' @param residual Residual data.frame with one column per age group.
#' @param frequency Seasonal frequency used by the tests.
#' @param digits Number of digits used to round p-values.
#'
#' @return Matrix of p-values for QS, Friedman, and Kruskal--Wallis tests.
  test_residual_seasonality <- function(residual, frequency = 52,
                                        digits = 3) {
#' DLNM diagnostics internal: safely extract one seasonality p-value
#'
#' @param x Numeric residual vector.
#' @param test_fun Seasonality test function from \pkg{seastests}.
#'
#' @return Numeric p-value, or \code{NA_real_} when the test cannot be run.
    safe_p_value <- function(x, test_fun) {
      x <- x[is.finite(x)]
      if (length(x) < frequency || length(unique(x)) < 2) {
        return(NA_real_)
      }
      tryCatch(round(test_fun(x, freq = frequency)$Pval, digits),
#' DLNM diagnostics internal: ignore failed seasonality tests
#'
#' @param e Error object from \code{tryCatch()}.
#'
#' @return \code{NA_real_}.
               error = function(e) NA_real_)
    }

    test_table <- rbind(
      QS = sapply(residual, safe_p_value, test_fun = seastests::qs),
      Friedman = sapply(residual, safe_p_value, test_fun = seastests::fried),
      `Kruskal-Wallis` = sapply(residual, safe_p_value,
                                test_fun = seastests::kw)
    )
    colnames(test_table) <- age_groups
    test_table
  }

#' DLNM diagnostics internal: fit one regional DLNM--LC model
#'
#' @param region Legacy region label.
#'
#' @return A list containing the fitted model, residuals, and seasonality tests.
  fit_dlnm_lc <- function(region) {
    cfg <- region_config[[region]]
    rate <- make_rate(cfg)
    wave <- make_wave(region)
    utci <- load_region_utci(region)

    model <- DLNM_LC(rate,
                       make_lc_crossbasis(utci),
                       wave,
                       tol = tol,
                       max_iter = cfg$max_iter,
                       region = region)
    residual <- calculate_residual(rate, model$log_fitted)

    list(model = model,
         residual = residual,
         seasonality_tests = test_residual_seasonality(residual))
  }

#' DLNM diagnostics internal: fit the multi-region DLNM--LL model
#'
#' @return A list containing the fitted model, residuals, and seasonality tests.
  fit_dlnm_ll <- function() {
    rate_list <- lapply(region_config, make_rate)
    wave_list <- lapply(names(region_config), make_wave)
    names(wave_list) <- names(region_config)
    utci_list <- lapply(names(region_config), load_region_utci)
    names(utci_list) <- names(region_config)
    cb_list <- lapply(utci_list, make_ll_crossbasis)

    model <- DLNM_LL(dat_list = rate_list,
                       cb_list = cb_list,
                       wave_list = wave_list,
                       tol = tol,
                       max_iter = 20)
    residual <- lapply(names(rate_list), function(region) {
      calculate_residual(rate_list[[region]], model$log_fitted[[region]])
    })
    names(residual) <- names(rate_list)

    list(model = model,
         residual = residual,
         seasonality_tests = lapply(residual, test_residual_seasonality))
  }

  lc <- lapply(names(region_config), fit_dlnm_lc)
  names(lc) <- names(region_config)

  list(
    DLNM_LC = list(
      results = lc,
      seasonality_tests = lapply(lc, `[[`, "seasonality_tests"),
      display_names = vapply(region_config, `[[`, character(1),
                             "display_region")
    ),
    DLNM_LL = c(
      fit_dlnm_ll(),
      list(display_names = vapply(region_config, `[[`, character(1),
                                  "display_region"))
    )
  )
}
