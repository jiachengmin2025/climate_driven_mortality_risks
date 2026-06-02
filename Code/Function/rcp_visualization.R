#' Load RCP simulation output for plotting
#'
#' @param model Model label, usually \code{"LC"} or \code{"LL"}.
#' @param scenario RCP scenario label.
#' @param data_dir Directory containing simulation \code{.RData} files.
#' @param annualized Logical; if \code{TRUE}, load annualized simulation files.
#'
#' @return Nested simulation list with display region names.
load_rcp_simulation <- function(model, scenario, data_dir, annualized = FALSE) {
  model <- toupper(model)
  file_name <- if (annualized) {
    paste0("final_", model, ".sim.annualized.", scenario, ".RData")
  } else {
    paste0("final_", model, "_", scenario, ".RData")
  }
  object_name <- if (annualized) {
    paste0("final_", model, ".sim.annualized.", scenario)
  } else {
    paste0("final_", model, ".sim.", scenario)
  }

  file_path <- file.path(data_dir, file_name)
  if (!file.exists(file_path)) {
    stop("Missing simulation output: ", file_path)
  }

  env <- new.env(parent = emptyenv())
  load(file_path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    object_name <- ls(env)[[1]]
  }

  to_output_region_names(env[[object_name]])
}

#' Summarise one RCP simulation age group for plotting
#'
#' @param sim_data Nested simulation list for one scenario.
#' @param region Display region name.
#' @param age_group Age-group column name.
#' @param scenario_label Label shown in the plot legend.
#' @param start_index First simulated time index to display.
#' @param start_year Calendar year corresponding to \code{start_index}.
#' @param annualized Logical; if \code{TRUE}, treat simulations as annualized
#'   log rates.
#'
#' @return Data.frame with time, mean, lower, upper, and scenario label columns.
summarise_rcp_age_group <- function(sim_data, region, age_group, scenario_label,
                                    start_index, start_year,
                                    annualized = FALSE) {
  region_data <- sim_data[[region]]
  if (is.null(region_data)) {
    stop("Region ", region, " is missing from simulation data.")
  }

  sim_matrix <- do.call(
    rbind,
    lapply(region_data, function(simulation) simulation[[age_group]])
  )
  sim_matrix <- sim_matrix[, start_index:ncol(sim_matrix), drop = FALSE]

  if (!annualized) {
    sim_matrix <- exp(sim_matrix)
  }

  out <- data.frame(
    time = if (annualized) {
      seq(start_year, start_year + ncol(sim_matrix) - 1)
    } else {
      start_year + (seq_len(ncol(sim_matrix)) - 1) / 52
    },
    mean = colMeans(sim_matrix),
    lower = apply(sim_matrix, 2, quantile, probs = 0.025),
    upper = apply(sim_matrix, 2, quantile, probs = 0.975),
    Scenario = scenario_label,
    stringsAsFactors = FALSE
  )

  if (annualized) {
    out$mean <- exp(out$mean)
    out$lower <- exp(out$lower)
    out$upper <- exp(out$upper)
  }

  out
}

#' Plot one age-group RCP simulation summary
#'
#' @param plot_data Data.frame returned by \code{summarise_rcp_age_group()}.
#' @param age_group Age-group column name.
#' @param config Plot configuration list.
#' @param age_labels Named vector mapping age-group column names to labels.
#' @param scenario_colors Named vector of colors for RCP scenarios.
#'
#' @return A \code{ggplot2} plot object.
plot_rcp_age_group <- function(plot_data, age_group, config,
                               age_labels, scenario_colors) {
  ggplot(plot_data, aes(x = time, color = Scenario)) +
    geom_line(aes(y = lower), linetype = "dashed", size = config$interval_size) +
    geom_line(aes(y = upper), linetype = "dashed", size = config$interval_size) +
    geom_line(aes(y = mean), size = config$mean_size) +
    scale_color_manual(values = scenario_colors, name = "Scenario") +
    labs(
      subtitle = paste("Age group", age_labels[[age_group]]),
      x = "Year",
      y = "Mortality Rate"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = config$major_y_grid,
      panel.grid.minor = element_blank(),
      plot.subtitle = element_text(
        size = config$subtitle_size,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.ticks = element_line(color = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = config$legend_text_size),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    scale_x_continuous(
      breaks = config$x_breaks,
      limits = config$x_limits,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = config$y_expand)
}

#' Plot all age groups for one region and model
#'
#' @param datasets Named list of scenario simulation outputs.
#' @param region Display region name.
#' @param config Plot configuration list including output path, size, and axis
#'   settings.
#' @param start_index First simulated time index to display.
#' @param start_year Calendar year corresponding to \code{start_index}.
#' @param age_labels Named vector mapping age-group names to labels.
#' @param scenario_colors Named vector of colors for scenarios.
#' @param annualized Logical; if \code{TRUE}, plot annualized simulations.
#' @param save_plots Logical; if \code{TRUE}, save the combined plot as PDF.
#' @param dpi Resolution passed to \code{ggplot2::ggsave()}.
#'
#' @return Combined \code{patchwork} plot object.
plot_rcp_region_model <- function(datasets, region, config,
                                  start_index, start_year,
                                  age_labels, scenario_colors,
                                  annualized = FALSE,
                                  save_plots = TRUE,
                                  dpi = 1000) {
  age_plots <- lapply(names(age_labels), function(age_group) {
    plot_data <- bind_rows(lapply(names(datasets), function(scenario_label) {
      summarise_rcp_age_group(
        sim_data = datasets[[scenario_label]],
        region = region,
        age_group = age_group,
        scenario_label = scenario_label,
        start_index = start_index,
        start_year = start_year,
        annualized = annualized
      )
    }))
    plot_rcp_age_group(plot_data, age_group, config, age_labels, scenario_colors)
  })

  combined_plot <- wrap_plots(age_plots, ncol = config$ncol) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  if (save_plots) {
    dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(
      filename = file.path(config$output_dir, config$filenames[[region]]),
      plot = combined_plot,
      width = config$width,
      height = config$height,
      units = "in",
      dpi = dpi
    )
  }

  combined_plot
}
