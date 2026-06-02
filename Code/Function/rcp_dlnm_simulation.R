#' Bootstrap DLNM coefficient vectors
#'
#' @param dlnm_coef Nested list of DLNM coefficient vectors by region,
#'   iteration, and age group.
#' @param dlnm_vcov Nested list of DLNM variance--covariance matrices aligned
#'   with \code{dlnm_coef}.
#' @param nsim Number of bootstrap draws.
#'
#' @return Nested list of simulated coefficient matrices by region, iteration,
#'   and age group.
bootstrap_dlnm_coefficients <- function(dlnm_coef, dlnm_vcov, nsim) {
  boot_coef <- lapply(names(dlnm_coef), function(region) {
    lapply(names(dlnm_coef[[region]]), function(iteration) {
      lapply(names(dlnm_coef[[region]][[iteration]]), function(age_group) {
        MASS::mvrnorm(
          n = nsim,
          mu = dlnm_coef[[region]][[iteration]][[age_group]],
          Sigma = dlnm_vcov[[region]][[iteration]][[age_group]]
        )
      })
    })
  })
  names(boot_coef) <- names(dlnm_coef)

  for (region in names(boot_coef)) {
    names(boot_coef[[region]]) <- names(dlnm_coef[[region]])
    for (iteration in names(boot_coef[[region]])) {
      names(boot_coef[[region]][[iteration]]) <- names(dlnm_coef[[region]][[iteration]])
    }
  }

  boot_coef
}

#' Sum DLNM simulation terms over backfitting iterations
#'
#' @param iteration_terms Nested list of simulated DLNM linear predictors by
#'   iteration and age group.
#'
#' @return List of summed simulated DLNM terms by age group.
sum_iteration_simulations <- function(iteration_terms) {
  age_group_sums <- list()
  for (iteration in iteration_terms) {
    for (age_group_name in names(iteration)) {
      if (!age_group_name %in% names(age_group_sums)) {
        age_group_sums[[age_group_name]] <- iteration[[age_group_name]]
      } else {
        age_group_sums[[age_group_name]] <- mapply(
          `+`,
          age_group_sums[[age_group_name]],
          iteration[[age_group_name]],
          SIMPLIFY = FALSE
        )
      }
    }
  }
  age_group_sums
}

#' Simulate future DLNM linear predictor terms
#'
#' @param boot_coef Bootstrapped DLNM coefficients from
#'   \code{bootstrap_dlnm_coefficients()}.
#' @param cb_list_future Named list of future UTCI cross-basis matrices.
#' @param wave_list_future Named list of future hot/cold wave covariates.
#'
#' @return Nested list of simulated DLNM terms by region and age group.
simulate_dlnm_terms <- function(boot_coef, cb_list_future, wave_list_future) {
  out <- lapply(names(boot_coef), function(region) {
    region_cb <- as.matrix(cb_list_future[[region]])
    region_wave <- wave_list_future[[region]]
    has_wave <- !is.null(region_wave) && ncol(as.data.frame(region_wave)) > 0
    region_wave <- if (has_wave) as.matrix(region_wave) else NULL
    n_cb <- ncol(region_cb)

    iteration_terms <- lapply(boot_coef[[region]], function(iteration) {
      lapply(iteration, function(age_group_boot) {
        lapply(seq_len(nrow(age_group_boot)), function(bootstrap_index) {
          beta <- age_group_boot[bootstrap_index, ]
          eta <- beta[[1]] + as.numeric(region_cb %*% beta[seq_len(n_cb) + 1])
          wave_beta <- beta[-seq_len(n_cb + 1)]
          if (has_wave && length(wave_beta) > 0) {
            eta <- eta + as.numeric(region_wave %*% wave_beta)
          }
          eta
        })
      })
    })

    sum_iteration_simulations(iteration_terms)
  })
  names(out) <- names(boot_coef)
  out
}

#' Split simulated mortality matrices by age group
#'
#' @param region_sim_matrices List of simulation matrices, each with columns
#'   corresponding to age groups.
#'
#' @return List of simulations by age group.
split_sim_matrices_by_age <- function(region_sim_matrices) {
  n_age_groups <- ncol(region_sim_matrices[[1]])
  age_group_sims <- vector("list", n_age_groups)
  names(age_group_sims) <- rcp_age_names

  for (sim_index in seq_along(region_sim_matrices)) {
    for (age_group_index in seq_len(n_age_groups)) {
      age_group_sims[[age_group_index]][[sim_index]] <-
        region_sim_matrices[[sim_index]][, age_group_index]
    }
  }

  age_group_sims
}

#' Convert age-group simulation lists to data.frames
#'
#' @param age_list_by_region Nested list of simulation vectors by region and
#'   age group.
#'
#' @return Nested list of simulation data.frames by display region name.
age_lists_to_simulation_dfs <- function(age_list_by_region) {
  simulation_dfs <- lapply(age_list_by_region, function(region) {
    simulations <- lapply(seq_along(region[[1]]), function(sim_index) {
      sim_data <- sapply(rcp_age_names, function(age_group) {
        region[[age_group]][[sim_index]]
      })
      data.frame(sim_data, check.names = FALSE) |>
        stats::setNames(rcp_age_names)
    })
    names(simulations) <- paste0("simulation_", seq_along(simulations))
    simulations
  })
  to_output_region_names(simulation_dfs)
}

#' Combine converged LC/LL components with simulated DLNM terms
#'
#' @param converged_age_sims Simulated converged LC/LL log-rate components.
#' @param dlnm_age_sims Simulated DLNM log-rate components.
#'
#' @return Nested list of combined log-rate simulation paths.
combine_converged_and_dlnm <- function(converged_age_sims, dlnm_age_sims) {
  Map(function(region_converged, region_dlnm) {
    Map(function(age_group_converged, age_group_dlnm) {
      mapply(function(sim_converged, sim_dlnm) {
        n_converged <- length(sim_converged)
        n_dlnm <- length(sim_dlnm)

        if (n_dlnm >= n_converged) {
          return(head(sim_converged + tail(sim_dlnm, n_converged), n_converged))
        }

        unchanged_part <- sim_converged[seq_len(n_converged - n_dlnm)]
        adjusted_part <- tail(sim_converged, n_dlnm) + sim_dlnm
        c(unchanged_part, adjusted_part)
      }, age_group_converged, age_group_dlnm, SIMPLIFY = FALSE)
    }, region_converged, region_dlnm)
  }, converged_age_sims, dlnm_age_sims)
}

#' Add Gaussian simulation error to log-rate simulation paths
#'
#' @param simulation_dfs Nested list of simulation data.frames.
#' @param sd Standard deviation of the Gaussian error.
#' @param seed Optional random seed.
#'
#' @return Nested list of simulation data.frames with additive errors.
add_simulation_errors <- function(simulation_dfs, sd, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  lapply(simulation_dfs, function(region) {
    lapply(region, function(simulation) {
      error <- matrix(
        rnorm(prod(dim(simulation)), mean = 0, sd = sd),
        nrow = nrow(simulation),
        ncol = ncol(simulation),
        dimnames = list(NULL, colnames(simulation))
      )
      as.data.frame(as.matrix(simulation) + error)
    })
  })
}

#' Run one DLNM--LC RCP mortality simulation scenario
#'
#' @param scenario RCP scenario label.
#' @param DLNM_LC_models Named list of fitted regional \code{DLNM_LC()} models.
#' @param dlnm_coef Nested DLNM coefficient list.
#' @param dlnm_vcov Nested DLNM variance--covariance list.
#' @param base_dir Project root.
#' @param simulation_weeks Number of weeks to simulate.
#' @param nsim Number of simulation paths.
#' @param error_sd Standard deviation of additive simulation error.
#' @param seed Random seed.
#'
#' @return One-row data.frame identifying the saved simulation object and file.
run_lc_rcp_scenario <- function(scenario, DLNM_LC_models, dlnm_coef, dlnm_vcov,
                                base_dir, simulation_weeks, nsim,
                                error_sd = 0.01, seed = 111) {
  message("Running legacy-compatible DLNM-LC RCP simulation for ", scenario)

  set.seed(seed)
  kappa_sim <- simulate_legacy_lc_kappa_paths(
    n = simulation_weeks,
    nsim = nsim
  )

  converged_age_sims <- lapply(rcp_internal_region_names, function(region) {
    model <- DLNM_LC_models[[region]]
    region_matrices <- lapply(seq_len(ncol(kappa_sim[[region]])), function(sim_index) {
      lc_part <- outer(kappa_sim[[region]][, sim_index],
                       model$final_LC$b_x,
                       FUN = "*")
      sweep(lc_part, 2, model$final_LC$a_x, FUN = "+")
    })
    split_sim_matrices_by_age(region_matrices)
  })
  names(converged_age_sims) <- rcp_internal_region_names

  boot_coef <- bootstrap_dlnm_coefficients(
    dlnm_coef = dlnm_coef,
    dlnm_vcov = dlnm_vcov,
    nsim = nsim
  )

  future_inputs <- load_future_rcp_inputs(
    base_dir,
    scenario = scenario,
    lag_days = 21,
    var_df = 2,
    lag_df = 2,
    degree = 3,
    legacy_lag_knots = TRUE
  )

  dlnm_age_sims <- simulate_dlnm_terms(
    boot_coef = boot_coef,
    cb_list_future = future_inputs$cb_list,
    wave_list_future = future_inputs$wave_list
  )

  final_sim <- combine_converged_and_dlnm(
    converged_age_sims = converged_age_sims,
    dlnm_age_sims = dlnm_age_sims
  )
  final_sim <- age_lists_to_simulation_dfs(final_sim)
  final_sim <- add_simulation_errors(final_sim, sd = error_sd)

  out_dir <- file.path(base_dir, "Results/Simulated_mortality_rates")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  final_name <- paste0("final_LC.sim.", scenario)
  output_file <- file.path(out_dir, paste0("final_LC_", scenario, ".RData"))
  assign(final_name, final_sim)
  save(list = final_name, file = output_file, envir = environment())

  data.frame(
    scenario = scenario,
    object = final_name,
    output_file = output_file,
    row.names = NULL
  )
}

#' Run one DLNM--LL RCP mortality simulation scenario
#'
#' @param scenario RCP scenario label.
#' @param DLNM_LL_fitting Fitted \code{DLNM_LL()} object.
#' @param base_dir Project root.
#' @param simulation_weeks Number of weeks to simulate.
#' @param nsim Number of simulation paths.
#' @param error_sd Standard deviation of additive simulation error.
#' @param seed Random seed.
#'
#' @return One-row data.frame identifying the saved simulation object and file.
run_ll_rcp_scenario <- function(scenario, DLNM_LL_fitting, base_dir,
                                simulation_weeks, nsim,
                                error_sd = 0.0001, seed = 111) {
  message("Running legacy-compatible DLNM-LL RCP simulation for ", scenario)

  set.seed(seed)
  kappa_paths <- simulate_legacy_ll_kappa_paths(
    n = simulation_weeks,
    nsim = nsim
  )
  K_sim <- kappa_paths$Common
  kappa_sim <- kappa_paths$kappa

  Axs <- DLNM_LL_fitting$final_LL$Ax
  Bxs <- DLNM_LL_fitting$final_LL$Bx
  bxs <- DLNM_LL_fitting$final_LL$bx

  common_component <- lapply(seq_len(ncol(K_sim)), function(sim_index) {
    outer(K_sim[, sim_index], Bxs, FUN = "*")
  })

  regional_component <- lapply(rcp_internal_region_names, function(region) {
    lapply(seq_len(ncol(kappa_sim[[region]])), function(sim_index) {
      outer(kappa_sim[[region]][, sim_index], bxs[[region]], FUN = "*")
    })
  })
  names(regional_component) <- rcp_internal_region_names

  converged_age_sims <- lapply(rcp_internal_region_names, function(region) {
    region_matrices <- Map(
      function(common_part, regional_part) {
        sweep(common_part + regional_part, 2, Axs[[region]], FUN = "+")
      },
      common_component,
      regional_component[[region]]
    )
    split_sim_matrices_by_age(region_matrices)
  })
  names(converged_age_sims) <- rcp_internal_region_names

  boot_coef <- bootstrap_dlnm_coefficients(
    dlnm_coef = DLNM_LL_fitting$dlnm_coef,
    dlnm_vcov = DLNM_LL_fitting$dlnm_vcov,
    nsim = nsim
  )

  future_inputs <- load_future_rcp_inputs(
    base_dir,
    scenario = scenario,
    lag_days = 21,
    var_df = 3,
    lag_df = 3,
    degree = 3
  )

  dlnm_age_sims <- simulate_dlnm_terms(
    boot_coef = boot_coef,
    cb_list_future = future_inputs$cb_list,
    wave_list_future = future_inputs$wave_list
  )

  final_sim <- combine_converged_and_dlnm(
    converged_age_sims = converged_age_sims,
    dlnm_age_sims = dlnm_age_sims
  )
  final_sim <- age_lists_to_simulation_dfs(final_sim)
  final_sim <- add_simulation_errors(final_sim, sd = error_sd)

  out_dir <- file.path(base_dir, "Results/Simulated_mortality_rates")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  final_name <- paste0("final_LL.sim.", scenario)
  output_file <- file.path(out_dir, paste0("final_LL_", scenario, ".RData"))
  assign(final_name, final_sim)
  save(list = final_name, file = output_file, envir = environment())

  data.frame(
    scenario = scenario,
    object = final_name,
    output_file = output_file,
    row.names = NULL
  )
}
