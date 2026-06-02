#' Legacy DLNM--LC SARIMA specifications for RCP simulation
#'
#' @return Named list of SARIMA/ARIMA parameter specifications for regional
#'   DLNM--LC \code{kappa(t)} simulations.
legacy_lc_sarima_specs <- function() {
  list(
    Attiki = list(engine = "arima",
                  ar = c(1.3039156, -0.3653115),
                  ma = -0.9733007),
    Lisbon = list(engine = "sarima",
                  ar = 0.52663241, d = 0, ma = 0,
                  sar = 0.08157594, D = 0, sma = 0, S = 52),
    Roma = list(engine = "sarima",
                ar = c(1.0683698, -0.1036257), d = 0,
                ma = -0.8438067, sar = 0.0964450, D = 0,
                sma = 0, S = 52)
  )
}

#' Legacy DLNM--LL SARIMA specifications for RCP simulation
#'
#' @return Named list of SARIMA/ARIMA parameter specifications for the common
#'   \code{K(t)} and regional \code{kappa(t, i)} simulations.
legacy_ll_sarima_specs <- function() {
  list(
    Common = list(engine = "sarima",
                  ar = 0.757553149, d = 0,
                  ma = c(-0.319981635, 0.009969013),
                  sar = 0, D = 0, sma = 0.106674692, S = 52),
    Attiki = list(engine = "sarima",
                  ar = c(0.40358338, 0.25462489), d = 0,
                  ma = c(-0.04433412, -0.30213904),
                  sar = 0, D = 0, sma = 0.14377324, S = 52),
    Lisbon = list(engine = "arima",
                  ar = 0.3361887),
    Roma = list(engine = "sarima",
                ar = c(1.53866985, -0.55864782), d = 0,
                ma = c(-1.38872694, 0.43411636),
                sar = -0.06684183, D = 0, sma = 0, S = 52)
  )
}

#' Simulate paths from a legacy ARIMA/SARIMA specification
#'
#' @param spec List defining the simulation engine and ARIMA/SARIMA
#'   coefficients.
#' @param n Number of time points to simulate in each path.
#' @param nsim Number of simulation paths.
#'
#' @return Numeric matrix with \code{n} rows and \code{nsim} columns.
simulate_legacy_sarima_paths <- function(spec, n, nsim) {
  replicate(n = nsim, {
    if (identical(spec$engine, "arima")) {
      model <- list()
      if (!is.null(spec$ar)) model$ar <- spec$ar
      if (!is.null(spec$ma)) model$ma <- spec$ma
      as.numeric(stats::arima.sim(n = n, model = model))
    } else {
      as.numeric(astsa::sarima.sim(
        n = n,
        ar = spec$ar,
        d = spec$d,
        ma = spec$ma,
        sar = spec$sar,
        D = spec$D,
        sma = spec$sma,
        S = spec$S
      ))
    }
  })
}

#' Simulate DLNM--LC regional kappa paths
#'
#' @param n Number of time points to simulate.
#' @param nsim Number of simulation paths.
#'
#' @return Named list of simulated regional \code{kappa(t)} matrices.
simulate_legacy_lc_kappa_paths <- function(n, nsim) {
  specs <- legacy_lc_sarima_specs()
  setNames(
    lapply(rcp_internal_region_names, function(region) {
      simulate_legacy_sarima_paths(specs[[region]], n = n, nsim = nsim)
    }),
    rcp_internal_region_names
  )
}

#' Simulate DLNM--LL common and regional kappa paths
#'
#' @param n Number of time points to simulate.
#' @param nsim Number of simulation paths.
#'
#' @return A list containing \code{Common} for \code{K(t)} paths and
#'   \code{kappa} for region-specific \code{kappa(t, i)} paths.
simulate_legacy_ll_kappa_paths <- function(n, nsim) {
  specs <- legacy_ll_sarima_specs()
  list(
    Common = simulate_legacy_sarima_paths(specs$Common, n = n, nsim = nsim),
    kappa = setNames(
      lapply(rcp_internal_region_names, function(region) {
        simulate_legacy_sarima_paths(specs[[region]], n = n, nsim = nsim)
      }),
      rcp_internal_region_names
    )
  )
}

#' Build a table of legacy SARIMA coefficients
#'
#' @param model_type Model label: \code{"DLNM_LC"}, \code{"DLNM_LL"},
#'   \code{"LC"}, or \code{"LL"}.
#'
#' @return Data.frame of SARIMA coefficient names and estimates.
legacy_sarima_coefficient_table <- function(model_type) {
  specs <- switch(
    model_type,
    DLNM_LC = legacy_lc_sarima_specs(),
    DLNM_LL = legacy_ll_sarima_specs(),
    LC = legacy_lc_sarima_specs(),
    LL = legacy_ll_sarima_specs(),
    stop("Unknown legacy SARIMA model type: ", model_type)
  )

  do.call(rbind, lapply(names(specs), function(series_name) {
    spec <- specs[[series_name]]
    parameters <- setdiff(names(spec), "engine")
    do.call(rbind, lapply(parameters, function(parameter) {
      values <- spec[[parameter]]
      data.frame(
        model = model_type,
        time_series = series_name,
        parameter = paste0(parameter, seq_along(values)),
        estimate = as.numeric(values),
        row.names = NULL
      )
    }))
  }))
}

#' Write legacy SARIMA inputs used by RCP simulation
#'
#' @param model_type Model label passed to
#'   \code{legacy_sarima_coefficient_table()}.
#' @param scenario RCP scenario label.
#' @param base_dir Project root.
#'
#' @return Data.frame of SARIMA coefficients. Also writes a CSV under
#'   \code{Results/Simulated_mortality_rates}.
write_legacy_sarima_inputs <- function(model_type, scenario, base_dir) {
  out_dir <- file.path(base_dir, "Results/Simulated_mortality_rates")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  table <- legacy_sarima_coefficient_table(model_type)
  write.csv(
    table,
    file.path(out_dir, paste0("sarima_inputs_", model_type, "_", scenario, ".csv")),
    row.names = FALSE
  )
  table
}
