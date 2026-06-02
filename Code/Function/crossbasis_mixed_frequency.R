#' Replace missing values with neighboring means
#'
#' @param data Numeric vector or matrix containing possible \code{NA} values.
#'
#' @return Object with the same shape as \code{data}, with each missing entry
#'   replaced by the mean of its nearest available neighbors.
replace_na_with_mean <- function(data) {
  for (i in seq_along(data)) {
    if (is.na(data[i])) {
      left_index <- max(1, i - 1)
      right_index <- min(length(data), i + 1)
      data[i] <- mean(data[c(left_index, right_index)], na.rm = TRUE)
    }
  }
  data
}

#' Extract one day-of-week UTCI series from weekly groups
#'
#' @param daily_city Data.frame with \code{Week} and \code{value} columns.
#' @param day_index Integer day position within each ISO week.
#' @param rows Row indices of weekly groups to retain.
#'
#' @return Numeric vector of selected day-specific weekly UTCI values.
get_weekly_day <- function(daily_city, day_index, rows) {
  weekly_day <- dplyr::summarise(
    dplyr::group_by(daily_city, Week),
    lag = dplyr::nth(value, day_index),
    .groups = "drop"
  )
  weekly_day <- dplyr::select(weekly_day, lag)
  weekly_day <- dplyr::slice(weekly_day, rows)
  unlist(weekly_day, use.names = FALSE)
}

#' Shift a weekly vector by whole weeks
#'
#' @param x Numeric vector.
#' @param n_shift Number of leading weekly shifts to insert as \code{NA}.
#'
#' @return Numeric vector shifted by \code{n_shift} positions.
shift_weekly <- function(x, n_shift) {
  c(rep(NA_real_, n_shift), x[seq_len(length(x) - n_shift)])
}

#' Build weekly UTCI lag matrix from daily UTCI data
#' Converts daily UTCI values into the 22-column weekly lag matrix used by the
#' historical DLNM cross-basis construction.
#'
#' @param daily_mean Data.frame containing daily UTCI series and a \code{Date}
#'   column.
#' @param source_col Column name for the region-specific daily UTCI series.
#'
#' @return Matrix with 260 weekly rows and columns \code{lag0} to \code{lag21}.
build_utci_lag_matrix <- function(daily_mean, source_col) {
  daily_city <- dplyr::mutate(daily_mean, Week = ISOweek::ISOweek(Date))
  daily_city <- dplyr::transmute(daily_city, Week, value = .data[[source_col]])
  daily_city <- dplyr::mutate(
    daily_city,
    value = dplyr::if_else(dplyr::row_number() <= 4, NA_real_, value)
  )

  weekly_mean <- dplyr::summarise(
    dplyr::group_by(daily_city, Week),
    value = mean(value),
    .groups = "drop"
  )

  out <- matrix(NA_real_, nrow = 260, ncol = 22)
  colnames(out) <- paste0("lag", 0:21)
  rownames(out) <- weekly_mean$Week[2:261]

  lag0 <- get_weekly_day(daily_city, 7, 2:261)
  lag1 <- get_weekly_day(daily_city, 6, 1:260)
  lag2 <- get_weekly_day(daily_city, 5, 1:260)
  lag3 <- get_weekly_day(daily_city, 4, 1:260)
  lag4 <- get_weekly_day(daily_city, 3, 1:260)
  lag5 <- get_weekly_day(daily_city, 2, 1:260)
  lag6 <- get_weekly_day(daily_city, 1, 1:260)

  out[, 1:7] <- cbind(lag0, lag1, lag2, lag3, lag4, lag5, lag6)
  out[, 8:13] <- cbind(
    shift_weekly(lag0, 1), shift_weekly(lag1, 1),
    shift_weekly(lag2, 1), shift_weekly(lag3, 1),
    shift_weekly(lag4, 1), shift_weekly(lag5, 1)
  )
  # Preserve the original script's construction of lag13 and lag20.
  out[, 14] <- shift_weekly(lag5, 1)
  out[, 15:20] <- cbind(
    shift_weekly(lag0, 2), shift_weekly(lag1, 2),
    shift_weekly(lag2, 2), shift_weekly(lag3, 2),
    shift_weekly(lag4, 2), shift_weekly(lag5, 2)
  )
  out[, 21] <- shift_weekly(lag5, 2)
  out[, 22] <- shift_weekly(lag0, 3)

  replace_na_with_mean(out)
}

#' Save a value under a chosen object name
#'
#' @param object_name Character name to assign inside the saved \code{.RData}.
#' @param object_value Object to save.
#' @param file Output \code{.RData} path.
#'
#' @return No direct return value. Writes an \code{.RData} file.
save_named_object <- function(object_name, object_value, file) {
  env <- new.env(parent = emptyenv())
  assign(object_name, object_value, envir = env)
  save(list = object_name, file = file, envir = env)
}

#' Build and save historical UTCI cross-basis lag matrices
#' Reads daily UTCI data and saves both display-name and legacy-name
#' historical lag matrices under \code{Results/Crossbasis_matrix}.
#'
#' @param base_dir Project root containing \code{Data/UTCI_data/Daily_data}.
#' @param write_legacy_name_outputs Logical; whether to save legacy files such
#'   as \code{cb_Attiki.RData} and \code{cb_Roma.RData}.
#' @param write_display_name_outputs Logical; whether to save display-name files
#'   such as \code{cb_Athens.RData} and \code{cb_Rome.RData}.
#'
#' @return Invisibly returns a named list of historical UTCI lag matrices.
build_crossbasis_matrices <- function(base_dir,
                                      write_legacy_name_outputs = TRUE,
                                      write_display_name_outputs = TRUE) {
  daily_dir <- file.path(base_dir, "Data/UTCI_data/Daily_data")
  crossbasis_dir <- file.path(base_dir, "Results/Crossbasis_matrix")
  dir.create(crossbasis_dir, recursive = TRUE, showWarnings = FALSE)

  city_specs <- data.frame(
    display_name = c("Athens", "Lisbon", "Rome"),
    source_col = c("Attiki", "Lisbon", "Roma"),
    legacy_name = c("Attiki", "Lisbon", "Roma"),
    stringsAsFactors = FALSE
  )

  UTCI_daily_mean <- readxl::read_xlsx(file.path(daily_dir, "UTCI_daily_mean.xlsx"))
  crossbasis_list <- setNames(lapply(seq_len(nrow(city_specs)), function(i) {
    build_utci_lag_matrix(UTCI_daily_mean, city_specs$source_col[[i]])
  }), city_specs$display_name)

  UTCI_ext.Athens <- crossbasis_list$Athens
  UTCI_ext.Lisbon <- crossbasis_list$Lisbon
  UTCI_ext.Rome <- crossbasis_list$Rome
  save(crossbasis_list, UTCI_ext.Athens, UTCI_ext.Lisbon, UTCI_ext.Rome,
       file = file.path(crossbasis_dir, "Crossbasis_matrix.RData"))

  if (write_display_name_outputs) {
    for (city in names(crossbasis_list)) {
      save_named_object(
        paste0("UTCI_ext.", city),
        crossbasis_list[[city]],
        file.path(crossbasis_dir, paste0("cb_", city, ".RData"))
      )
    }
  }

  if (write_legacy_name_outputs) {
    for (i in seq_len(nrow(city_specs))) {
      save_named_object(
        paste0("UTCI_ext.", city_specs$legacy_name[[i]]),
        crossbasis_list[[city_specs$display_name[[i]]]],
        file.path(crossbasis_dir, paste0("cb_", city_specs$legacy_name[[i]], ".RData"))
      )
    }
  }

  invisible(crossbasis_list)
}

#' Ensure historical UTCI cross-basis files exist
#'
#' @param base_dir Project root.
#' @param quiet Logical; if \code{TRUE}, suppress messages from rebuilding.
#'
#' @return Invisibly returns the \code{Results/Crossbasis_matrix} directory.
ensure_crossbasis_matrix <- function(base_dir, quiet = TRUE) {
  crossbasis_dir <- file.path(base_dir, "Results/Crossbasis_matrix")
  expected_files <- file.path(
    crossbasis_dir,
    c("cb_Athens.RData", "cb_Lisbon.RData", "cb_Rome.RData",
      "cb_Attiki.RData", "cb_Roma.RData", "Crossbasis_matrix.RData")
  )

  if (all(file.exists(expected_files))) {
    return(invisible(crossbasis_dir))
  }

  if (quiet) {
    invisible(capture.output(build_crossbasis_matrices(base_dir)))
  } else {
    build_crossbasis_matrices(base_dir)
  }

  missing_files <- expected_files[!file.exists(expected_files)]
  if (length(missing_files) > 0) {
    stop("Could not generate required cross-basis files: ",
         paste(basename(missing_files), collapse = ", "))
  }

  invisible(crossbasis_dir)
}

#' Load historical UTCI lag matrix for one region
#'
#' @param base_dir Project root.
#' @param display_region Display region name, e.g. \code{"Athens"}.
#' @param legacy_region Legacy region name, e.g. \code{"Attiki"}.
#' @param quiet Logical; passed to \code{ensure_crossbasis_matrix()}.
#'
#' @return Matrix-like historical UTCI lag data for the requested region.
load_crossbasis_history <- function(base_dir, display_region, legacy_region = display_region,
                                    quiet = TRUE) {
  crossbasis_dir <- ensure_crossbasis_matrix(base_dir, quiet = quiet)
  file_candidates <- unique(file.path(
    crossbasis_dir,
    paste0("cb_", c(display_region, legacy_region), ".RData")
  ))
  cb_file <- file_candidates[file.exists(file_candidates)][1]

  if (is.na(cb_file)) {
    stop("No cross-basis history file found for ", display_region,
         ". Checked: ", paste(basename(file_candidates), collapse = ", "))
  }

  env <- new.env(parent = emptyenv())
  object_names <- load(cb_file, envir = env)
  preferred_names <- unique(paste0("UTCI_ext.", c(display_region, legacy_region)))
  preferred_names <- preferred_names[preferred_names %in% object_names]

  if (length(preferred_names) > 0) {
    return(as.matrix(env[[preferred_names[[1]]]]))
  }

  matrix_like <- object_names[vapply(object_names, function(object_name) {
    object <- env[[object_name]]
    is.matrix(object) || is.data.frame(object) || inherits(object, "crossbasis")
  }, logical(1))]

  if (length(matrix_like) == 0) {
    stop("No matrix-like UTCI lag history found in ", cb_file)
  }

  as.matrix(env[[matrix_like[[1]]]])
}
