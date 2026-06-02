#' Read historical combined mortality/exposure data for RCP projection
#'
#' @param base_dir Project root containing \code{Data/Combined_data}.
#'
#' @return A list of combined data tables by age group.
read_rcp_combined_data <- function(base_dir) {
  data_dir <- file.path(base_dir, "Data/Combined_data")
  list(
    Y20_64 = readxl::read_xlsx(file.path(data_dir, "Y20_64_combined.xlsx")),
    Y65_74 = readxl::read_xlsx(file.path(data_dir, "Y65_74_combined.xlsx")),
    Y75_84 = readxl::read_xlsx(file.path(data_dir, "Y75_84_combined.xlsx")),
    Y_GE85 = readxl::read_xlsx(file.path(data_dir, "Y_GE85_combined.xlsx"))
  )
}

#' Build historical mortality-rate and wave inputs for RCP projection
#'
#' @param base_dir Project root.
#'
#' @return A list:
#' \itemize{
#'   \item \code{dat_list}: annualized mortality-rate matrices by legacy region.
#'   \item \code{wave_list}: hot/cold wave covariate data.frames by legacy region.
#' }
build_historical_inputs <- function(base_dir) {
  combined <- read_rcp_combined_data(base_dir)
  death_cols <- c(Attiki = 3, Lisbon = 4, Roma = 6)
  exposure_cols <- c(Attiki = 9, Lisbon = 10, Roma = 12)

  dat_list <- lapply(rcp_internal_region_names, function(region) {
    death <- cbind(
      combined$Y20_64[[death_cols[[region]]]],
      combined$Y65_74[[death_cols[[region]]]],
      combined$Y75_84[[death_cols[[region]]]],
      combined$Y_GE85[[death_cols[[region]]]]
    )
    exposure <- cbind(
      combined$Y20_64[[exposure_cols[[region]]]],
      combined$Y65_74[[exposure_cols[[region]]]],
      combined$Y75_84[[exposure_cols[[region]]]],
      combined$Y_GE85[[exposure_cols[[region]]]]
    )
    out <- death / exposure * 52
    colnames(out) <- rcp_age_names
    out
  })
  names(dat_list) <- rcp_internal_region_names

  wave_list <- lapply(rcp_internal_region_names, function(region) {
    out <- data.frame(
      hot_wave = combined$Y20_64[[paste0(region, "_hot_wave3")]],
      cold_wave = combined$Y20_64[[paste0(region, "_cold_wave3")]]
    )
    colnames(out) <- c(paste0(region, "_hot_wave3"),
                       paste0(region, "_cold_wave3"))
    out
  })
  names(wave_list) <- rcp_internal_region_names

  list(dat_list = dat_list, wave_list = wave_list)
}

#' Build a UTCI cross-basis matrix for RCP inputs
#'
#' @param x UTCI lag matrix or exposure vector.
#' @param lag_days Maximum lag in days.
#' @param var_df Degrees of freedom for the exposure dimension.
#' @param lag_df Degrees of freedom for the lag dimension.
#' @param degree Degree passed to \code{dlnm::equalknots()}.
#' @param legacy_lag_knots Logical; if \code{TRUE}, use the legacy
#'   \code{dlnm::logknots(lag, df)} call.
#'
#' @return A \code{dlnm::crossbasis} object.
make_rcp_crossbasis_matrix <- function(x, lag_days = 21, var_df = 3,
                                       lag_df = 3, degree = 3,
                                       legacy_lag_knots = FALSE) {
  varknots <- dlnm::equalknots(x, fun = "ns", df = var_df, degree = degree)
  lagknots <- if (legacy_lag_knots) {
    dlnm::logknots(lag_days, lag_df)
  } else {
    dlnm::logknots(lag_days, fun = "ns", df = lag_df)
  }

  dlnm::crossbasis(
    x,
    lag = lag_days,
    argvar = list(fun = "ns", knots = varknots),
    arglag = list(knots = lagknots, df = lag_df)
  )
}

#' Load historical UTCI cross-basis matrices for all regions
#'
#' @param base_dir Project root.
#' @param lag_days Maximum lag in days.
#' @param var_df Degrees of freedom for the exposure dimension.
#' @param lag_df Degrees of freedom for the lag dimension.
#' @param degree Degree passed to \code{dlnm::equalknots()}.
#' @param legacy_lag_knots Logical; passed to
#'   \code{make_rcp_crossbasis_matrix()}.
#'
#' @return Named list of historical cross-basis objects by legacy region.
load_historical_cb_list <- function(base_dir, lag_days = 21, var_df = 3,
                                    lag_df = 3, degree = 3,
                                    legacy_lag_knots = FALSE) {
  if (!exists("ensure_crossbasis_matrix", mode = "function")) {
    source(file.path(base_dir, "Code/Function/crossbasis_mixed_frequency.R"))
  }
  ensure_crossbasis_matrix(base_dir)

  setNames(lapply(rcp_internal_region_names, function(region) {
    make_rcp_crossbasis_matrix(
      load_crossbasis_history(base_dir,
                              display_region = rcp_display_names[[region]],
                              legacy_region = region),
      lag_days = lag_days,
      var_df = var_df,
      lag_df = lag_df,
      degree = degree,
      legacy_lag_knots = legacy_lag_knots
    )
  }), rcp_internal_region_names)
}

#' Load future RCP UTCI cross-basis and wave inputs
#'
#' @param base_dir Project root.
#' @param scenario RCP scenario label, e.g. \code{"rcp26"} or \code{"rcp85"}.
#' @param lag_days Maximum lag in days.
#' @param var_df Degrees of freedom for the exposure dimension.
#' @param lag_df Degrees of freedom for the lag dimension.
#' @param degree Degree passed to \code{dlnm::equalknots()}.
#' @param legacy_lag_knots Logical; passed to
#'   \code{make_rcp_crossbasis_matrix()}.
#'
#' @return A list with \code{cb_list} and \code{wave_list}, both named by
#'   legacy region.
load_future_rcp_inputs <- function(base_dir, scenario, lag_days = 21,
                                   var_df = 3, lag_df = 3, degree = 3,
                                   legacy_lag_knots = FALSE) {
  scenario_dir <- file.path(base_dir, "Data/Simulation_data/UTCI", scenario)

  cb_list <- setNames(lapply(rcp_internal_region_names, function(region) {
    env <- new.env(parent = emptyenv())
    load(file.path(scenario_dir, paste0(region, "_cb_", scenario, ".RData")),
         envir = env)
    make_rcp_crossbasis_matrix(
      env[[paste0(region, "_prep_cb_", scenario)]],
      lag_days = lag_days,
      var_df = var_df,
      lag_df = lag_df,
      degree = degree,
      legacy_lag_knots = legacy_lag_knots
    )
  }), rcp_internal_region_names)

  wave_env <- new.env(parent = emptyenv())
  load(file.path(scenario_dir, paste0("wave_list_", scenario, ".RData")),
       envir = wave_env)
  wave_list <- to_internal_region_names(wave_env[[paste0("wave_list_", scenario)]])
  wave_list <- wave_list[rcp_internal_region_names]

  list(cb_list = cb_list, wave_list = wave_list)
}
