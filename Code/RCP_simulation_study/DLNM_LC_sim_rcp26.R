# Load packages
packages <- c("dlnm", "splines", "forecast", "astsa", "MASS")
invisible(lapply(packages, library, character.only = TRUE))

# Load functions
files <- list.files("Code/Function/", pattern = "\\.R$", full.names = TRUE)
sapply(files, source)


# Load historical UTCi and mortality data
load("Code/Historical_data_aggregation/Historical_data.RData")


# Load historical crossbasis matrix data
# *** Before 2D smoothing ***
load("Code/Crossbasis_matrix/cb_Attiki.RData")
load("Code/Crossbasis_matrix/cb_Lisbon.RData")
load("Code/Crossbasis_matrix/cb_Roma.RData")

# Attiki
varknots <- equalknots(UTCI_ext.Attiki, fun="ns",df = 2, degree = 3)
lagknots <- logknots(21, 2)
Attiki_cb <- crossbasis(UTCI_ext.Attiki, lag=21, 
                        argvar=list(fun="ns", knots = varknots), 
                        arglag=list(knots=lagknots,df = 2))

# Lisbon
varknots <- equalknots(UTCI_ext.Lisbon,fun="ns",df = 2, degree = 3)
lagknots <- logknots(21, 2)
Lisbon_cb <- crossbasis(UTCI_ext.Lisbon, lag=21, 
                        argvar=list(fun="ns", knots = varknots), 
                        arglag=list(knots=lagknots,df = 2))

# Roma
varknots <- equalknots(UTCI_ext.Roma,fun="ns",df = 2, degree = 3)
lagknots <- logknots(21, 2)
Roma_cb <- crossbasis(UTCI_ext.Roma, lag=21, 
                      argvar=list(fun="ns", knots = varknots), 
                      arglag=list(knots=lagknots,df = 2))

cb_list <- list(Attiki = Attiki_cb, Lisbon = Lisbon_cb, Roma = Roma_cb)
rm(UTCI_ext.Attiki, UTCI_ext.Lisbon, UTCI_ext.Roma)
rm(Attiki_cb, Lisbon_cb, Roma_cb)


# Load future UTCI data and future corssbasis matrix data (rcp26)
# *** Before 2D smoothing ***
load("Data/Simulation_data/UTCI/rcp26/Attiki_cb_rcp26.RData")
load("Data/Simulation_data/UTCI/rcp26/Lisbon_cb_rcp26.RData")
load("Data/Simulation_data/UTCI/rcp26/Roma_cb_rcp26.RData")


# Attiki
varknots <- equalknots(Attiki_prep_cb_rcp26, fun="ns",df = 2, degree = 3)
lagknots <- logknots(21, 2)
Attiki_cb.rcp26 <- crossbasis(Attiki_prep_cb_rcp26, lag=21, 
                              argvar=list(fun="ns", knots = varknots), 
                              arglag=list(knots=lagknots,df = 2))

# Lisbon
varknots <- equalknots(Lisbon_prep_cb_rcp26,fun="ns",df = 2, degree = 3)
lagknots <- logknots(21, 2)
Lisbon_cb.rcp26 <- crossbasis(Lisbon_prep_cb_rcp26, lag=21, 
                              argvar=list(fun="ns", knots = varknots), 
                              arglag=list(knots=lagknots,df = 2))

# Roma
varknots <- equalknots(Roma_prep_cb_rcp26,fun="ns",df = 2, degree = 3)
lagknots <- logknots(21, 2)
Roma_cb.rcp26 <- crossbasis(Roma_prep_cb_rcp26, lag=21, 
                            argvar=list(fun="ns", knots = varknots), 
                            arglag=list(knots=lagknots,df = 2))
cb_list_rcp26 = list(Attiki = Attiki_cb.rcp26, Lisbon = Lisbon_cb.rcp26, Roma = Roma_cb.rcp26)
rm(Attiki_prep_cb_rcp26, Lisbon_prep_cb_rcp26, Roma_prep_cb_rcp26)
rm(Attiki_cb.rcp26, Lisbon_cb.rcp26, Roma_cb.rcp26)


# Load Future Heat/Cold Wave Data (rcp26)
load("Data/Simulation_data/UTCI/rcp26/wave_list_rcp26.RData")


# Fit Backfitting Lee-Carter Model
B_DLNM_LC_Attiki = B_DLNM_LC(dat_list$Attiki, cb_list$Attiki, wave_list$Attiki, tol = 1e-2, max_iter = 20, region = "Attiki")
B_DLNM_LC_Lisbon = B_DLNM_LC(dat_list$Lisbon, cb_list$Lisbon, wave_list$Lisbon, tol = 1e-2, max_iter = 20, region = "Lisbon")
B_DLNM_LC_Roma = B_DLNM_LC(dat_list$Roma, cb_list$Roma, wave_list$Roma, tol = 1e-2, max_iter = 20, region = "Roma")
nsim = 10000


# Fit SARIMA model on kappa(t)
kt_fit_Attiki = kt_model(B_DLNM_LC_Attiki$final_LC$k_t, "SARIMA")
kt_fit_Lisbon = kt_model(B_DLNM_LC_Lisbon$final_LC$k_t, "SARIMA")
kt_fit_Roma = kt_model(B_DLNM_LC_Roma$final_LC$k_t, "SARIMA")


# Simulation of K(t) and kappa(t,i) (Mannualy extract coefficeints and orders)
set.seed(111)
K_sim_list = list(
  Attiki_k.sim = replicate(n = nsim, arima.sim(n = 52*26, list(ar = c(1.3039156, -0.3653115), ma = -0.9733007))),
  Lisbon_k.sim = replicate(n = nsim, sarima.sim(n = 52*26, ar = 0.52663241, d = 0, ma = 0, sar = 0.08157594, D = 0, sma = 0, S = 52)),
  Roma_k.sim = replicate(n = nsim, sarima.sim(n = 52*26, ar = c(1.0683698, -0.1036257), d = 0, ma = -0.8438067, sar = 0.0964450, D = 0, S = 52))
)


# Converged LC
axs = list(Attiki = B_DLNM_LC_Attiki$final_LC$a_x, Lisbon = B_DLNM_LC_Lisbon$final_LC$a_x, Roma = B_DLNM_LC_Roma$final_LC$a_x)
bxs = list(Attiki = B_DLNM_LC_Attiki$final_LC$b_x, Lisbon = B_DLNM_LC_Lisbon$final_LC$b_x, Roma = B_DLNM_LC_Roma$final_LC$b_x)
kts = K_sim_list

bk.sim = list()
for (ind in 1:length(kts)){
  bk.sim[[ind]] = lapply(1:ncol(kts[[ind]]), function(col) {
    outer(kts[[ind]][, col], bxs[[ind]], FUN = "*")
  })
}
names(bk.sim) = c("Attiki", "Lisbon", "Roma")

converged_LC.sim.rcp26 = Map(function(sim_sublist, axs_values) {
  # sublist
  lapply(sim_sublist, function(df) {
    # add A(x, i)
    sweep(df, 2, axs_values, FUN = "+")
  })
}, bk.sim, axs)

converged_LC.sim.rcp26 <- lapply(converged_LC.sim.rcp26, function(region) {
  n_age_groups <- ncol(region[[1]])
  age_group_sims <- vector("list", n_age_groups)
  names(age_group_sims) <- c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
  for (sim_index in seq_along(region)) {
    for (age_group_index in seq_len(n_age_groups)) {
      if (is.null(age_group_sims[[age_group_index]])) {
        age_group_sims[[age_group_index]] <- list()
      }
      age_group_sims[[age_group_index]][[sim_index]] <- region[[sim_index]][, age_group_index]
    }
  }
  age_group_sims
})


# DLNM Coefficient Bootstrap Simulation
dlnm_coef = list(Attiki = B_DLNM_LC_Attiki$dlnm_coef, Lisbon = B_DLNM_LC_Lisbon$dlnm_coef, Roma = B_DLNM_LC_Roma$dlnm_coef)
dlnm_vcov = list(Attiki = B_DLNM_LC_Attiki$dlnm_vcov, Lisbon = B_DLNM_LC_Lisbon$dlnm_vcov, Roma = B_DLNM_LC_Roma$dlnm_vcov)
boot_coef <- lapply(names(dlnm_coef), function(region) {
  lapply(names(dlnm_coef[[region]]), function(iteration) {
    lapply(names(dlnm_coef[[region]][[iteration]]), function(age_group) {
      coef_vec <- dlnm_coef[[region]][[iteration]][[age_group]]
      vcov_mat <- dlnm_vcov[[region]][[iteration]][[age_group]]
      mvrnorm(n = nsim, mu = coef_vec, Sigma = vcov_mat)
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



# DLNM Part Simulation
## Crossbasis matrix simulation results
## (1) dlnm_coef %*% smoothing crossbasis matrix (each iteration each age group)
cb_matrix.sim = lapply(seq_along(boot_coef), function(i) { 
  region_cb <- cb_list_rcp26[[i]]  # crossbasis matrix
  
  lapply(boot_coef[[i]], function(iteration) {
    lapply(iteration, function(age_group) {
      # 100 simulations
      lapply(seq_len(nrow(age_group)), function(bootstrap_index) {
        coef <- age_group[bootstrap_index, ]
        apply(region_cb, 1, function(row) {
          sum(row * coef[2:(length(coef) - 2)]) + coef[1]
        })
      })
    })
  })
})
names(cb_matrix.sim) = names(boot_coef)

## (1.1) Sum up all forecast values from DLNMs.
cb_matrix.sum.sim <- lapply(cb_matrix.sim, function(region) {
  age_group_sums <- list()
  for (iteration in region) {
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
})


## Heat/Cold Wave Simulation Results
## (2) dlnm_coef %*% hot/cold wave term (each iteration each age group)
wave.sim = lapply(seq_along(boot_coef), function(i) { 
  region_wave <- wave_list_rcp26[[i]] # regional wave data
  lapply(boot_coef[[i]], function(iteration) {
    lapply(iteration, function(age_group) {
      # 100 simulations
      lapply(seq_len(nrow(age_group)), function(bootstrap_index) {
        coef <- age_group[bootstrap_index, ]
        apply(region_wave, 1, function(row) {
          sum(row * coef[(length(coef)-1):length(coef)])
        })
      })
    })
  })
})
names(wave.sim) = names(boot_coef)
## (2.1) Sum up all forecast values from DLNMs.
wave.sum.sim <- lapply(wave.sim, function(region) {
  age_group_sums <- list()
  for (iteration in region) {
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
})


## Summation of Crossbasis Matrix And Wave Results
dlnm.sum.sim <- Map(function(cb_region, wave_region) {
  age_group_sums <- lapply(names(cb_region), function(age_group_name) {
    cb_age_group <- cb_region[[age_group_name]]
    wave_age_group <- wave_region[[age_group_name]]
    
    mapply(function(cb_row, wave_row) {
      cb_row + wave_row
    }, cb_age_group, wave_age_group, SIMPLIFY = FALSE)
  })
  
  names(age_group_sums) <- names(cb_region)
  age_group_sums
}, cb_matrix.sum.sim, wave.sum.sim)


## Summation of Converged LC Part And DLNM Part
final_LC.sim.rcp26<- Map(function(region_converged, region_dlnm) {
  Map(function(age_group_converged, age_group_dlnm) {
    mapply(function(sim_converged, sim_dlnm) {
      # UTCI is unknown
      unchanged_part <- sim_converged[1:(length(sim_converged) - length(sim_dlnm))]
      # UTCI is known
      adjusted_part <- sim_converged[(length(sim_converged) - length(sim_dlnm) + 1):length(sim_converged)] + sim_dlnm
      c(unchanged_part, adjusted_part)
    }, age_group_converged, age_group_dlnm, SIMPLIFY = FALSE)
  }, region_converged, region_dlnm)
}, converged_LC.sim.rcp26, dlnm.sum.sim)

age_group_names <- c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
final_LC.sim.rcp26 <- lapply(final_LC.sim.rcp26, function(region) {
  lapply(seq_along(region[[1]]), function(sim_index) {
    sim_data <- sapply(region, function(age_group) age_group[[sim_index]])
    data.frame(sim_data, check.names = FALSE) |> 
      setNames(age_group_names)
  })
})


# Add Error Terms, Store Converged LC And B-DLNM-LC Simulation Results 
# Generate error terms
generate_error_terms <- function(original_data, sd = 0.01) {
  # original_data: Reference structure for generating error terms
  # sd: Standard deviation of the error terms
  lapply(original_data, function(region) {
    lapply(region, function(simulation) {
      # Generate random error terms matching the dimensions of each simulation
      as.data.frame(matrix(rnorm(n = prod(dim(simulation)), mean = 0, sd = sd), 
                           nrow = nrow(simulation), ncol = ncol(simulation),
                           dimnames = list(NULL, colnames(simulation))))
    })
  })
}

# Generate error terms with the same structure as `final_LC.sim.rcp26` and `converged_LC.sim.rcp26`
error_terms <- generate_error_terms(final_LC.sim.rcp26, sd = 0.01)

# Add error terms to `final_LC.sim.rcp26`
final_LC.sim.rcp26 <- Map(function(data, error) {
  lapply(seq_along(data), function(i) {
    Map(`+`, data[[i]], error[[i]])  # Add error terms to corresponding data
  })
}, final_LC.sim.rcp26, error_terms)

# Rename the simulations within each region
final_LC.sim.rcp26 <- lapply(final_LC.sim.rcp26, function(region) {
  names(region) <- paste0("simulation_", seq_along(region))
  region
})

# Data Structure
final_LC.sim.rcp26 <- lapply(final_LC.sim.rcp26, function(region) {
  lapply(region, function(simulations) {
    as.data.frame(do.call(cbind, simulations))
  })
})

# Add error terms to `converged_LC.sim.rcp26`
converged_LC.sim.rcp26 <- Map(function(data, error) {
  lapply(seq_along(data), function(i) {
    Map(`+`, data[[i]], error[[i]])  # Add error terms to corresponding data
  })
}, converged_LC.sim.rcp26, error_terms)

# Rename the simulations within each region and rename age groups
converged_LC.sim.rcp26 <- lapply(converged_LC.sim.rcp26, function(region) {
  # Rename the age groups inside each region
  age_group_names <- c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")  # Ensure correct names for age groups
  names(region) <- age_group_names
  
  # Rename the simulations inside each age group
  lapply(region, function(age_group) {
    names(age_group) <- paste0("simulation_", seq_along(age_group))
    age_group
  })
})

converged_LC.sim.rcp26 <- lapply(converged_LC.sim.rcp26, function(region) {
  age_group_names <- names(region)
  num_simulations <- length(region[[1]])
  simulations <- lapply(seq_len(num_simulations), function(sim_index) {
    sim_data <- sapply(age_group_names, function(age_group) {
      region[[age_group]][[sim_index]] 
    })
    data.frame(sim_data) |> setNames(age_group_names)  
  })
  names(simulations) <- paste0("simulation_", seq_along(simulations))
  simulations
})

save(final_LC.sim.rcp26, file = "Results/Simulated_mortality_rates/final_LC_rcp26.RData")
