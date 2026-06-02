#' Load weekly RCP mortality simulation output
#'
#' @param model Model label, usually \code{"LC"} or \code{"LL"}.
#' @param scenario RCP scenario label.
#' @param simulation_dir Directory containing weekly simulation \code{.RData}
#'   files.
#'
#' @return Simulation list with display region names.
load_rcp_weekly_simulation <- function(model, scenario, simulation_dir) {
  file_path <- file.path(simulation_dir, paste0("final_", model, "_", scenario, ".RData"))
  if (!file.exists(file_path)) {
    stop("Missing weekly simulation output: ", file_path)
  }

  env <- new.env(parent = emptyenv())
  load(file_path, envir = env)

  object_name <- paste0("final_", model, ".sim.", scenario)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    object_name <- ls(env)[[1]]
  }

  to_output_region_names(env[[object_name]])
}

#' Average rows within fixed-size groups
#'
#' @param df Data.frame or matrix of weekly mortality rates.
#' @param group_size Number of rows per group, defaulting to 52 weeks.
#'
#' @return Data.frame of grouped means with standard RCP age-group names.
average_by_group <- function(df, group_size = 52) {
  groups <- split(df, ceiling(seq_len(nrow(df)) / group_size))
  out <- do.call(rbind, lapply(groups, function(group) colSums(group) / group_size))
  out <- as.data.frame(out)
  colnames(out) <- rcp_age_names
  out
}

#' Annualize weekly RCP mortality simulations
#' Converts weekly log-rate simulation paths to annualized log rates.
#'
#' @param simulation_list Nested weekly simulation list by region and simulation.
#' @param group_size Number of weekly rows per annualized group.
#'
#' @return Nested list of annualized log-rate simulation data.frames.
annualize_rcp_simulation <- function(simulation_list, group_size = 52) {
  simulation_list <- to_output_region_names(simulation_list)
  lapply(simulation_list, function(region) {
    lapply(region, function(simulation) {
      weekly_rates <- exp(as.data.frame(simulation))
      annualized_rates <- average_by_group(weekly_rates, group_size = group_size)
      log(annualized_rates)
    })
  })
}

#' Save annualized RCP simulation output
#'
#' @param annualized_list Nested annualized simulation list.
#' @param model Model label, usually \code{"LC"} or \code{"LL"}.
#' @param scenario RCP scenario label.
#' @param annualized_dir Directory where annualized \code{.RData} files are
#'   written.
#'
#' @return No direct return value. Writes an annualized simulation \code{.RData}
#'   file.
save_rcp_annualized_simulation <- function(annualized_list, model, scenario,
                                           annualized_dir) {
  object_name <- paste0("final_", model, ".sim.annualized.", scenario)
  assign(object_name, annualized_list)
  save(
    list = object_name,
    envir = environment(),
    file = file.path(
      annualized_dir,
      paste0("final_", model, ".sim.annualized.", scenario, ".RData")
    )
  )
}
