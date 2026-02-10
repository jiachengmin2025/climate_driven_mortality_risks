# Load packages
packages <- c("readxl", "ggplot2", "forecast", "dplyr", "gridExtra", "writexl", 
              "tseries", "hwwntest", "dlnm", "splines", "ISOweek", "mgcv", 
              "Metrics", "demography", "seastests", "MTS", "reshape2", "astsa",
              "MASS", "patchwork")
lapply(packages, library, character.only = TRUE)


# Load Converged LL And B-DLNM-LL Results
# Get a list of all .RData files in the folder
rdata_files <- list.files(path = "Results/Simulated_mortality_rates", pattern = "\\.RData$", full.names = TRUE)

# Load all .RData files using lapply
lapply(rdata_files, function(file) {
  load(file, envir = .GlobalEnv)
})


# DLNM-LC Results under Different Scenarios (Functions)
# Function to prepare data for plotting
prepare_plot_data <- function(sim_data, region, age_group, start_year) {
  # Extract data for the specified region and age group
  region_data <- sim_data[[region]]
  age_group_data <- lapply(region_data, function(sim) sim[[age_group]])
  
  # Combine all simulations into a matrix
  sim_matrix <- do.call(rbind, age_group_data)
  
  # Calculate mean and confidence intervals
  summary_df <- data.frame(
    time = start_year + (seq_len(ncol(sim_matrix)) - 1) / 52,  # Convert time points to years
    mean = colMeans(sim_matrix),
    lower = apply(sim_matrix, 2, quantile, probs = 0.025),
    upper = apply(sim_matrix, 2, quantile, probs = 0.975)
  )
  return(summary_df)
}

# List of datasets_final_LC and labels
datasets_final_LC <- list(
  "RCP 2.6" = final_LC.sim.rcp26,
  "RCP 8.5" = final_LC.sim.rcp85
)



## Attiki Results
# Define target region and age groups
target_region <- "Attiki"
age_groups <- c("Y20_64" = "20-64", "Y65_74" = "65-74", "Y75_84" = "75-84", "Y_GE85" = "85+")

# Set the start year (2031) and starting index (570)
start_year <- 2031
start_index <- 570

# Prepare plots for each age group
plots <- lapply(names(age_groups), function(age_group) {
  # Prepare data for each scenario
  plot_data <- lapply(names(datasets_final_LC), function(label) {
    # Prepare the data for the given scenario
    data_matrix <- do.call(rbind, lapply(datasets_final_LC[[label]][[target_region]], function(sim) sim[[age_group]]))
    # Filter the data to start from the given index
    filtered_matrix <- data_matrix[, start_index:ncol(data_matrix)]
    # Restore data by exponentiating log-transformed values
    restored_matrix <- exp(filtered_matrix)
    # Calculate mean and confidence intervals
    summary_df <- data.frame(
      mean = colMeans(restored_matrix),
      lower = apply(restored_matrix, 2, quantile, probs = 0.025),
      upper = apply(restored_matrix, 2, quantile, probs = 0.975),
      time = start_year + (seq_len(ncol(restored_matrix)) - 1) / 52  # Generate time points
    )
    summary_df$scenario <- label  # Add scenario label
    return(summary_df)
  })
  
  # Combine data from all scenarios into a single data frame
  combined_plot_data <- do.call(rbind, plot_data)
  
  # Generate the plot
  ggplot(combined_plot_data, aes(x = time, group = scenario)) +
    # Draw dashed lines for the lower and upper confidence interval boundaries instead of a ribbon
    geom_line(aes(y = lower, color = scenario), linetype = "dashed", size = 0.25) +
    geom_line(aes(y = upper, color = scenario), linetype = "dashed", size = 0.25) +
    # Draw the mean line
    geom_line(aes(y = mean, color = scenario), size = 0.25) +
    # Set labels and legends
    labs(
      title = NULL,
      subtitle = paste("Age group", age_groups[[age_group]]),
      x = "Year",
      y = "Mortality Rate",  # Reflect restored data
      color = "Scenario",
      fill = "Scenario"
    ) +
    # Manually set colors for each scenario
    scale_color_manual(values = c("RCP 2.6" = "#007bbd", 
                                  "RCP 4.5" = "#38bd48", 
                                  "RCP 6.0" = "#FF7F00", 
                                  "RCP 8.5" = "#B7282E")) +
    scale_fill_manual(values = c("RCP 2.6" = "#B3D1F5", 
                                 "RCP 4.5" = "#83d081", 
                                 "RCP 6.0" = "#FFC185", 
                                 "RCP 8.5" = "#FFE3D9")) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      plot.margin = margin(5, 5, 5, 5),
      axis.text.x = element_text(hjust = 0.5),
      axis.ticks = element_line(color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),  # Center the subtitle
      legend.position = "bottom",  # Move legend to the bottom
      legend.direction = "horizontal",  # Arrange legend elements in a row
      legend.box = "horizontal"  # Ensure elements are placed in a single row
    ) +
    scale_x_continuous(
      limits = c(start_year, 2045),  # Set x-axis limits from 2031 to 2045
      breaks = c(2031, 2035, 2040, 2045),
      expand = c(0, 0)
    ) +
    scale_y_continuous()
})

# Combine all age group plots into one layout
Attiki_plot <- wrap_plots(plots, ncol = 1) +
  plot_layout(guides = "collect") +  # Collect legends from subplots
  plot_annotation(
    #title = "Simulated Mortality Rates via DLNM-LC (Attica)",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Center the title
      legend.position = "bottom",  # Move legend to the bottom
      legend.direction = "horizontal"  # Arrange legend items horizontally
    )
  )

# Display the final plot
Attiki_plot

# Save the final plot
ggsave(filename = "Figures/RCP_Simulation/Overall/B_DLNM_LC/Attiki_all_scenarios_no_scale_line.pdf",
       plot = Attiki_plot, width = 6, height = 10, units = "in", dpi = 1000)


## Lisbon Results
# Define target region and age groups
target_region <- "Lisbon"
age_groups <- c("Y20_64" = "20-64", "Y65_74" = "65-74", "Y75_84" = "75-84", "Y_GE85" = "85+")

# Set the start year (2031) and starting index (570)
start_year <- 2031
start_index <- 570

# Prepare plots for each age group
plots <- lapply(names(age_groups), function(age_group) {
  # Prepare data for each scenario
  plot_data <- lapply(names(datasets_final_LC), function(label) {
    # Prepare the data for the given scenario
    data_matrix <- do.call(rbind, lapply(datasets_final_LC[[label]][[target_region]], function(sim) sim[[age_group]]))
    # Filter the data to start from the given index
    filtered_matrix <- data_matrix[, start_index:ncol(data_matrix)]
    # Restore data by exponentiating log-transformed values
    restored_matrix <- exp(filtered_matrix)
    # Calculate mean and confidence intervals
    summary_df <- data.frame(
      mean = colMeans(restored_matrix),
      lower = apply(restored_matrix, 2, quantile, probs = 0.025),
      upper = apply(restored_matrix, 2, quantile, probs = 0.975),
      time = start_year + (seq_len(ncol(restored_matrix)) - 1) / 52  # Generate time points
    )
    summary_df$scenario <- label  # Add scenario label
    return(summary_df)
  })
  
  # Combine data from all scenarios into a single data frame
  combined_plot_data <- do.call(rbind, plot_data)
  
  # Generate the plot
  ggplot(combined_plot_data, aes(x = time, group = scenario)) +
    # Draw dashed lines for the lower and upper confidence interval boundaries instead of a ribbon
    geom_line(aes(y = lower, color = scenario), linetype = "dashed", size = 0.25) +
    geom_line(aes(y = upper, color = scenario), linetype = "dashed", size = 0.25) +
    # Draw the mean line
    geom_line(aes(y = mean, color = scenario), size = 0.25) +
    # Set labels and legends
    labs(
      title = NULL,
      subtitle = paste("Age group", age_groups[[age_group]]),
      x = "Year",
      y = "Mortality Rate",  # Reflect restored data
      color = "Scenario",
      fill = "Scenario"
    ) +
    # Manually set colors for each scenario
    scale_color_manual(values = c("RCP 2.6" = "#007bbd", 
                                  "RCP 4.5" = "#38bd48", 
                                  "RCP 6.0" = "#FF7F00", 
                                  "RCP 8.5" = "#B7282E")) +
    scale_fill_manual(values = c("RCP 2.6" = "#B3D1F5", 
                                 "RCP 4.5" = "#83d081", 
                                 "RCP 6.0" = "#FFC185", 
                                 "RCP 8.5" = "#FFE3D9")) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      plot.margin = margin(5, 5, 5, 5),
      axis.text.x = element_text(hjust = 0.5),
      axis.ticks = element_line(color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),  # Center the subtitle
      legend.position = "bottom",  # Move legend to the bottom
      legend.direction = "horizontal",  # Arrange legend elements in a row
      legend.box = "horizontal"  # Ensure elements are placed in a single row
    ) +
    scale_x_continuous(
      limits = c(start_year, 2045),  # Set x-axis limits from 2031 to 2045
      breaks = c(2031, 2035, 2040, 2045),
      expand = c(0, 0)
    ) +
    scale_y_continuous()
})

# Combine all age group plots into one layout
Lisbon_plot <- wrap_plots(plots, ncol = 1) +
  plot_layout(guides = "collect") +  # Collect legends from subplots
  plot_annotation(
    #title = "Simulated Mortality Rates via DLNM-LC (Lisbon)",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Center the title
      legend.position = "bottom",  # Move legend to the bottom
      legend.direction = "horizontal"  # Arrange legend items horizontally
    )
  )

# Display the final plot
Lisbon_plot

# Save the final plot
ggsave(filename = "Figures/RCP_Simulation/Overall/B_DLNM_LC/Lisbon_all_scenarios_no_scale_line.pdf",
       plot = Lisbon_plot, width = 6, height = 10, units = "in", dpi = 1000)


## Roma Results

# Define target region and age groups
target_region <- "Roma"
age_groups <- c("Y20_64" = "20-64", "Y65_74" = "65-74", "Y75_84" = "75-84", "Y_GE85" = "85+")

# Set the start year (2031) and starting index (570)
start_year <- 2031
start_index <- 570

# Prepare plots for each age group
plots <- lapply(names(age_groups), function(age_group) {
  # Prepare data for each scenario
  plot_data <- lapply(names(datasets_final_LC), function(label) {
    # Prepare the data for the given scenario
    data_matrix <- do.call(rbind, lapply(datasets_final_LC[[label]][[target_region]], function(sim) sim[[age_group]]))
    # Filter the data to start from the given index
    filtered_matrix <- data_matrix[, start_index:ncol(data_matrix)]
    # Restore data by exponentiating log-transformed values
    restored_matrix <- exp(filtered_matrix)
    # Calculate mean and confidence intervals
    summary_df <- data.frame(
      mean = colMeans(restored_matrix),
      lower = apply(restored_matrix, 2, quantile, probs = 0.025),
      upper = apply(restored_matrix, 2, quantile, probs = 0.975),
      time = start_year + (seq_len(ncol(restored_matrix)) - 1) / 52  # Generate time points
    )
    summary_df$scenario <- label  # Add scenario label
    return(summary_df)
  })
  
  # Combine data from all scenarios into a single data frame
  combined_plot_data <- do.call(rbind, plot_data)
  
  # Generate the plot
  ggplot(combined_plot_data, aes(x = time, group = scenario)) +
    # Draw dashed lines for the lower and upper confidence interval boundaries instead of a ribbon
    geom_line(aes(y = lower, color = scenario), linetype = "dashed", size = 0.25) +
    geom_line(aes(y = upper, color = scenario), linetype = "dashed", size = 0.25) +
    # Draw the mean line
    geom_line(aes(y = mean, color = scenario), size = 0.25) +
    # Set labels and legends
    labs(
      title = NULL,
      subtitle = paste("Age group", age_groups[[age_group]]),
      x = "Year",
      y = "Mortality Rate",  # Reflect restored data
      color = "Scenario",
      fill = "Scenario"
    ) +
    # Manually set colors for each scenario
    scale_color_manual(values = c("RCP 2.6" = "#007bbd", 
                                  "RCP 4.5" = "#38bd48", 
                                  "RCP 6.0" = "#FF7F00", 
                                  "RCP 8.5" = "#B7282E")) +
    scale_fill_manual(values = c("RCP 2.6" = "#B3D1F5", 
                                 "RCP 4.5" = "#83d081", 
                                 "RCP 6.0" = "#FFC185", 
                                 "RCP 8.5" = "#FFE3D9")) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      plot.margin = margin(5, 5, 5, 5),
      axis.text.x = element_text(hjust = 0.5),
      axis.ticks = element_line(color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),  # Center the subtitle
      legend.position = "bottom",  # Move legend to the bottom
      legend.direction = "horizontal",  # Arrange legend elements in a row
      legend.box = "horizontal"  # Ensure elements are placed in a single row
    ) +
    scale_x_continuous(
      limits = c(start_year, 2045),  # Set x-axis limits from 2031 to 2045
      breaks = c(2031, 2035, 2040, 2045),
      expand = c(0, 0)
    ) +
    scale_y_continuous()
})

# Combine all age group plots into one layout
Roma_plot <- wrap_plots(plots, ncol = 1) +
  plot_layout(guides = "collect") +  # Collect legends from subplots
  plot_annotation(
    #title = "Simulated Mortality Rates via DLNM-LC (Rome)",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Center the title
      legend.position = "bottom",  # Move legend to the bottom
      legend.direction = "horizontal"  # Arrange legend items horizontally
    )
  )

# Display the final plot
Roma_plot

# Save the final plot
ggsave(filename = "Figures/RCP_Simulation/Overall/B_DLNM_LC/Roma_all_scenarios_no_scale_line.pdf",
       plot = Roma_plot, width = 6, height = 10, units = "in", dpi = 1000)

