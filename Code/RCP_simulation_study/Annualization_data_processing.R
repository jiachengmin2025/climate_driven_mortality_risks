# Load Converged LL and B-DLNM-LL Results
load("Results/Simulated_mortality_rates/final_LC_rcp26.RData")
load("Results/Simulated_mortality_rates/final_LC_rcp85.RData")
load("Results/Simulated_mortality_rates/final_LL_rcp26.RData")
load("Results/Simulated_mortality_rates/final_LL_rcp85.RData")


# Convert all to mortality rates
# Define a function to apply exp() to the innermost elements of a nested list
apply_exp_to_list <- function(nested_list) {
  lapply(nested_list, function(region) {
    lapply(region, function(sim) {
      exp(sim)
    })
  })
}

# Process all RCP lists
final_LC.sim.rcp26 <- apply_exp_to_list(final_LC.sim.rcp26)
final_LC.sim.rcp85 <- apply_exp_to_list(final_LC.sim.rcp85)

# Process all RCP lists
final_LL.sim.rcp26 <- apply_exp_to_list(final_LL.sim.rcp26)
final_LL.sim.rcp85 <- apply_exp_to_list(final_LL.sim.rcp85)


## Annualized Mortality Rates

# Function to sum every 52 rows in a dataframe
sum_by_group <- function(df, group_size = 52) {
  # Split the rows into groups of `group_size` and sum within each group
  grouped_data <- split(df, ceiling(seq_len(nrow(df)) / group_size))
  # Apply colSums to each group, then divide by group_size
  result <- do.call(rbind, lapply(grouped_data, function(group) colSums(group) / group_size))
  return(result)
}

# Apply the function to the nested list structure
annualized_simulations <- function(nested_list, group_size = 52) {
  lapply(nested_list, function(region) { # Iterate over regions
    lapply(region, function(simulation) { # Iterate over simulations
      sum_by_group(simulation, group_size) # Apply the grouping and summation
    })
  })
}

# Apply the aggregation function to your datasets
final_LC.sim.annualized.rcp26 <- annualized_simulations(final_LC.sim.rcp26, group_size = 52)
final_LC.sim.annualized.rcp85 <- annualized_simulations(final_LC.sim.rcp85, group_size = 52)

# Apply the aggregation function to your datasets
final_LL.sim.annualized.rcp26 <- annualized_simulations(final_LL.sim.rcp26, group_size = 52)
final_LL.sim.annualized.rcp85 <- annualized_simulations(final_LL.sim.rcp85, group_size = 52)


## Log Transform for Annualized Mortality Rates

# Define a function to apply log() to the innermost data.frames in the list structure
apply_log_to_list <- function(nested_list) {
  lapply(nested_list, function(region) {
    lapply(region, function(simulation) {
      log(simulation)
    })
  })
}

# Apply log() to all datasets
final_LC.sim.annualized.rcp26 <- apply_log_to_list(final_LC.sim.annualized.rcp26)
final_LC.sim.annualized.rcp85 <- apply_log_to_list(final_LC.sim.annualized.rcp85)

# Apply log() to all datasets
final_LL.sim.annualized.rcp26 <- apply_log_to_list(final_LL.sim.annualized.rcp26)
final_LL.sim.annualized.rcp85 <- apply_log_to_list(final_LL.sim.annualized.rcp85)

# Define the column names to assign
col_names <- c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")

# Function to assign column names to all data.frames in the nested list
# Function to convert simulations to data.frame and assign column names
assign_colnames <- function(nested_list, col_names) {
  lapply(nested_list, function(region) {
    lapply(region, function(simulation) {
      if (is.matrix(simulation) || is.double(simulation)) {
        # Convert to data.frame
        simulation <- as.data.frame(simulation)
        # Assign column names
        colnames(simulation) <- col_names
      }
      simulation  # Return the updated simulation
    })
  })
}

# Apply the function to your list
final_LC.sim.annualized.rcp26 <- assign_colnames(final_LC.sim.annualized.rcp26, col_names)
final_LC.sim.annualized.rcp85 <- assign_colnames(final_LC.sim.annualized.rcp85, col_names)

# Apply the function to your list
final_LL.sim.annualized.rcp26 <- assign_colnames(final_LL.sim.annualized.rcp26, col_names)
final_LL.sim.annualized.rcp85 <- assign_colnames(final_LL.sim.annualized.rcp85, col_names)


## Store Annualized Mortality Data (Log Scale)
save(final_LC.sim.annualized.rcp26,
     file = "Results/Simulated_annualized_mortality_rates/final_LC.sim.annualized.rcp26.RData")

save(final_LC.sim.annualized.rcp85,
     file = "Results/Simulated_annualized_mortality_rates/final_LC.sim.annualized.rcp85.RData")

save(final_LL.sim.annualized.rcp26, 
     file = "Results/Simulated_annualized_mortality_rates/final_LL.sim.annualized.rcp26.RData")

save(final_LL.sim.annualized.rcp85, 
     file = "Results/Simulated_annualized_mortality_rates/final_LL.sim.annualized.rcp85.RData")





