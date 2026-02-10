
# Load packages
packages <- c("readxl", "ggplot2", "forecast", "dplyr", "gridExtra", "writexl", 
              "tseries", "hwwntest", "dlnm", "splines", "ISOweek", "mgcv", 
              "Metrics", "demography", "seastests", "MTS", "reshape2", "astsa",
              "MASS", "patchwork")
lapply(packages, library, character.only = TRUE)


# Load Converged LL And DLNM-LL Results (RCP 4.5, RCP 8.5)
load("Results/Simulated_mortality_rates/final_LL_rcp26.RData")
load("Results/Simulated_mortality_rates/final_LL_rcp85.RData")


# DLNM-LL Results under Different Scenarios (Functions)
# Function to prepare data for plotting
prepare_plot_data <- function(sim_data, region, age_group, start_year) {
  # Extract data for the specified region and age group
  region_data <- sim_data[[region]]
  age_group_data <- lapply(region_data, function(sim) sim[[age_group]])
  
  # Combine all simulations into a matrix
  sim_matrix <- do.call(rbind, age_group_data)
  
  # CaLLulate mean and confidence intervals
  summary_df <- data.frame(
    time = start_year + (seq_len(ncol(sim_matrix)) - 1) / 52,  # Convert time points to years
    mean = colMeans(sim_matrix),
    lower = apply(sim_matrix, 2, quantile, probs = 0.025),
    upper = apply(sim_matrix, 2, quantile, probs = 0.975)
  )
  return(summary_df)
}

# List of datasets_final_LL and labels
datasets_final_LL <- list(
  "RCP 2.6" = final_LL.sim.rcp26,
  "RCP 8.5" = final_LL.sim.rcp85
)



## Attiki Results
# parameters
target_region <- "Attiki"
age_groups <- c(Y20_64 = "20-64", Y65_74 = "65-74", Y75_84 = "75-84", Y_GE85 = "85+")
start_year  <- 2031
start_index <- 570

plots <- lapply(names(age_groups), function(age_group) {
  df <- lapply(names(datasets_final_LL), function(scenario) {
    mat <- do.call(rbind,
                   lapply(datasets_final_LL[[scenario]][[target_region]], `[[`, age_group))
    mat <- mat[, start_index:ncol(mat), drop = FALSE]
    mat <- exp(mat)
    data.frame(
      time     = start_year + (seq_len(ncol(mat)) - 1) / 52,
      mean     = colMeans(mat),
      lower    = apply(mat, 2, quantile, 0.025),
      upper    = apply(mat, 2, quantile, 0.975),
      Scenario = scenario,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  ggplot(df, aes(x = time, color = Scenario)) +
    geom_line(aes(y = lower),  linetype = "dashed", size = 0.3) +
    geom_line(aes(y = upper),  linetype = "dashed", size = 0.3) +
    geom_line(aes(y = mean),   size = 0.4) +
    scale_color_manual(
      name   = "Scenario",
      values = c(
        "RCP 2.6" = "#007bbd",
        "RCP 8.5" = "#B7282E"
      )
    ) +
    labs(
      subtitle = paste("Age group", age_groups[[age_group]]),
      x        = "Year",
      y        = "Mortality Rate"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border       = element_rect(color = "black", fill = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray80"),
      plot.subtitle      = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title         = element_text(size = 12),
      axis.text          = element_text(size = 10),
      legend.position    = "bottom",
      legend.text        = element_text(size = 13),
      legend.title       = element_text(size = 14)
    ) +
    scale_x_continuous(
      breaks = c(2031, 2036, 2041, 2046),
      limits = c(start_year, 2046),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.2, 0.2))
    )
})

final_plot <- wrap_plots(plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(final_plot)

# Save the final plot
ggsave(filename = "Figures/RCP_Simulation/Overall/B_DLNM_LL/Attiki_weekly_sim.pdf",
       plot = final_plot, width = 10, height = 6, units = "in", dpi = 1000)


## Lisbon Results

# parameters
target_region <- "Lisbon"
age_groups <- c(Y20_64 = "20-64", Y65_74 = "65-74", Y75_84 = "75-84", Y_GE85 = "85+")
start_year  <- 2031
start_index <- 570

plots <- lapply(names(age_groups), function(age_group) {
  df <- lapply(names(datasets_final_LL), function(scenario) {
    mat <- do.call(rbind,
                   lapply(datasets_final_LL[[scenario]][[target_region]], `[[`, age_group))
    mat <- mat[, start_index:ncol(mat), drop = FALSE]
    mat <- exp(mat)
    data.frame(
      time     = start_year + (seq_len(ncol(mat)) - 1) / 52,
      mean     = colMeans(mat),
      lower    = apply(mat, 2, quantile, 0.025),
      upper    = apply(mat, 2, quantile, 0.975),
      Scenario = scenario,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  ggplot(df, aes(x = time, color = Scenario)) +
    geom_line(aes(y = lower),  linetype = "dashed", size = 0.3) +
    geom_line(aes(y = upper),  linetype = "dashed", size = 0.3) +
    geom_line(aes(y = mean),   size = 0.4) +
    scale_color_manual(
      name   = "Scenario",
      values = c(
        "RCP 2.6" = "#007bbd",
        "RCP 8.5" = "#B7282E"
      )
    ) +
    labs(
      subtitle = paste("Age group", age_groups[[age_group]]),
      x        = "Year",
      y        = "Mortality Rate"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border       = element_rect(color = "black", fill = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray80"),
      plot.subtitle      = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title         = element_text(size = 12),
      axis.text          = element_text(size = 10),
      legend.position    = "bottom",
      legend.text        = element_text(size = 13),
      legend.title       = element_text(size = 14)
    ) +
    scale_x_continuous(
      breaks = c(2031, 2036, 2041, 2046),
      limits = c(start_year, 2046),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.2, 0.2))
    )
})

final_plot <- wrap_plots(plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(final_plot)

# Save the final plot
ggsave(filename = "Figures/RCP_Simulation/Overall/B_DLNM_LL/Lisbon_weekly_sim.pdf",
       plot = final_plot, width = 10, height = 6, units = "in", dpi = 1000)


## Roma Results

# parameters
target_region <- "Roma"
age_groups <- c(Y20_64 = "20-64", Y65_74 = "65-74", Y75_84 = "75-84", Y_GE85 = "85+")
start_year  <- 2031
start_index <- 570

plots <- lapply(names(age_groups), function(age_group) {
  df <- lapply(names(datasets_final_LL), function(scenario) {
    mat <- do.call(rbind,
                   lapply(datasets_final_LL[[scenario]][[target_region]], `[[`, age_group))
    mat <- mat[, start_index:ncol(mat), drop = FALSE]
    mat <- exp(mat)
    data.frame(
      time     = start_year + (seq_len(ncol(mat)) - 1) / 52,
      mean     = colMeans(mat),
      lower    = apply(mat, 2, quantile, 0.025),
      upper    = apply(mat, 2, quantile, 0.975),
      Scenario = scenario,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  ggplot(df, aes(x = time, color = Scenario)) +
    geom_line(aes(y = lower),  linetype = "dashed", size = 0.3) +
    geom_line(aes(y = upper),  linetype = "dashed", size = 0.3) +
    geom_line(aes(y = mean),   size = 0.4) +
    scale_color_manual(
      name   = "Scenario",
      values = c(
        "RCP 2.6" = "#007bbd",
        "RCP 8.5" = "#B7282E"
      )
    ) +
    labs(
      subtitle = paste("Age group", age_groups[[age_group]]),
      x        = "Year",
      y        = "Mortality Rate"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border       = element_rect(color = "black", fill = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray80"),
      plot.subtitle      = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title         = element_text(size = 12),
      axis.text          = element_text(size = 10),
      legend.position    = "bottom",
      legend.text        = element_text(size = 13),
      legend.title       = element_text(size = 14)
    ) +
    scale_x_continuous(
      breaks = c(2031, 2036, 2041, 2046),
      limits = c(start_year, 2046),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.2, 0.2))
    )
})

final_plot <- wrap_plots(plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(final_plot)

# Save the final plot
ggsave(filename = "Figures/RCP_Simulation/Overall/B_DLNM_LL/Roma_weekly_sim.pdf",
       plot = final_plot, width = 10, height = 6, units = "in", dpi = 1000)

