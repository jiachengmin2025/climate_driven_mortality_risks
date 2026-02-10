## Packages
packages <- c("readxl", "dplyr", "tidyr", "zoo", "ISOweek", "ggplot2", "patchwork")
invisible(lapply(packages, library, character.only = TRUE))

## Load data
rcp26.mean = read_xlsx("Data/Simulation_data/UTCI/rcp26/UTCI_mean_rcp26.xlsx")[,-1]
rcp26.min = read_xlsx("Data/Simulation_data/UTCI/rcp26/UTCI_min_rcp26.xlsx")[,-1]
rcp26.max = read_xlsx("Data/Simulation_data/UTCI/rcp26/UTCI_max_rcp26.xlsx")[,-1]

rcp85.mean = read_xlsx("Data/Simulation_data/UTCI/rcp85/UTCI_mean_rcp85.xlsx")[,-1]
rcp85.min = read_xlsx("Data/Simulation_data/UTCI/rcp85/UTCI_min_rcp85.xlsx")[,-1]
rcp85.max = read_xlsx("Data/Simulation_data/UTCI/rcp85/UTCI_max_rcp85.xlsx")[,-1]



## Replace NA by Linear Interpolation
rcp26.min[-1] <- lapply(rcp26.min[-1], function(col) {
  na.approx(col, rule = 2)  # Linear interpolation
})

rcp26.mean[-1] <- lapply(rcp26.mean[-1], function(col) {
  na.approx(col, rule = 2)  # Linear interpolation
})

rcp26.max[-1] <- lapply(rcp26.max[-1], function(col) {
  na.approx(col, rule = 2)  # Linear interpolation
})

rcp85.min[-1] <- lapply(rcp85.min[-1], function(col) {
  na.approx(col, rule = 2)  # Linear interpolation
})

rcp85.mean[-1] <- lapply(rcp85.mean[-1], function(col) {
  na.approx(col, rule = 2)  # Linear interpolation
})

rcp85.max[-1] <- lapply(rcp85.max[-1], function(col) {
  na.approx(col, rule = 2)  # Linear interpolation
})


# RCP2.6
rcp26.mean$Week <- ISOweek(rcp26.mean$Time)
rcp26.min$Week <- ISOweek(rcp26.min$Time)
rcp26.max$Week <- ISOweek(rcp26.max$Time)


# RCP8.5
rcp85.mean$Week <- ISOweek(rcp85.mean$Time)
rcp85.min$Week <- ISOweek(rcp85.min$Time)
rcp85.max$Week <- ISOweek(rcp85.max$Time)



## Coverting the level
classify_utci <- function(x) {
  case_when(
    x < -13  ~ -1,
    x <=  32 ~  0,
    x >   32 ~  1,
    TRUE     ~ NA_real_
  )
}

# 3. Helper to compute weekly heat-/cold-wave counts
summarise_wave <- function(df_max, df_min) {
  hot <- df_max %>%
    mutate(
      Attiki_lvl = classify_utci(Attiki),
      Lisbon_lvl = classify_utci(Lisbon),
      Roma_lvl   = classify_utci(Roma),
      Attiki_hw = as.integer(Attiki_lvl == 1 & lag(Attiki_lvl)==1 & lag(Attiki_lvl,2)==1),
      Lisbon_hw = as.integer(Lisbon_lvl == 1 & lag(Lisbon_lvl)==1 & lag(Lisbon_lvl,2)==1),
      Roma_hw   = as.integer(Roma_lvl   == 1 & lag(Roma_lvl)  ==1 & lag(Roma_lvl,2)==1)
    ) %>%
    replace(is.na(.), 0) %>%
    mutate(Week = ISOweek(Date)) %>%
    group_by(Week) %>%
    summarise(
      Attiki_hot_wave3 = sum(Attiki_hw),
      Lisbon_hot_wave3 = sum(Lisbon_hw),
      Roma_hot_wave3   = sum(Roma_hw),
      .groups = "drop"
    )
  
  cold <- df_min %>%
    mutate(
      Attiki_lvl = classify_utci(Attiki),
      Lisbon_lvl = classify_utci(Lisbon),
      Roma_lvl   = classify_utci(Roma),
      Attiki_cw = as.integer(Attiki_lvl == -1 & lag(Attiki_lvl)==-1 & lag(Attiki_lvl,2)==-1),
      Lisbon_cw = as.integer(Lisbon_lvl == -1 & lag(Lisbon_lvl)==-1 & lag(Lisbon_lvl,2)==-1),
      Roma_cw   = as.integer(Roma_lvl   == -1 & lag(Roma_lvl)  ==-1 & lag(Roma_lvl,2)==-1)
    ) %>%
    replace(is.na(.), 0) %>%
    mutate(Week = ISOweek(Date)) %>%
    group_by(Week) %>%
    summarise(
      Attiki_cold_wave3 = sum(Attiki_cw),
      Lisbon_cold_wave3 = sum(Lisbon_cw),
      Roma_cold_wave3   = sum(Roma_cw),
      .groups = "drop"
    )
  
  left_join(hot, cold, by = "Week")
}



## Considering the consecutive days

# 2. Simplified UTCI classifier
classify_utci <- function(x) {
  case_when(
    x < -13  ~ -1,
    x <=  32 ~  0,
    x >   32 ~  1,
    TRUE     ~ NA_real_
  )
}

# 3. Helper to compute weekly heat-/cold-wave counts
summarise_wave <- function(df_max, df_min) {
  hot <- df_max %>%
    mutate(
      Attiki_lvl = classify_utci(Attiki),
      Lisbon_lvl = classify_utci(Lisbon),
      Roma_lvl   = classify_utci(Roma),
      Attiki_hw = as.integer(Attiki_lvl == 1 & lag(Attiki_lvl)==1 & lag(Attiki_lvl,2)==1),
      Lisbon_hw = as.integer(Lisbon_lvl == 1 & lag(Lisbon_lvl)==1 & lag(Lisbon_lvl,2)==1),
      Roma_hw   = as.integer(Roma_lvl   == 1 & lag(Roma_lvl)  ==1 & lag(Roma_lvl,2)==1)
    ) %>%
    replace(is.na(.), 0) %>%
    mutate(Week = ISOweek(Time)) %>%
    group_by(Week) %>%
    summarise(
      Attiki_hot_wave3 = sum(Attiki_hw),
      Lisbon_hot_wave3 = sum(Lisbon_hw),
      Roma_hot_wave3   = sum(Roma_hw),
      .groups = "drop"
    )
  
  cold <- df_min %>%
    mutate(
      Attiki_lvl = classify_utci(Attiki),
      Lisbon_lvl = classify_utci(Lisbon),
      Roma_lvl   = classify_utci(Roma),
      Attiki_cw = as.integer(Attiki_lvl == -1 & lag(Attiki_lvl)==-1 & lag(Attiki_lvl,2)==-1),
      Lisbon_cw = as.integer(Lisbon_lvl == -1 & lag(Lisbon_lvl)==-1 & lag(Lisbon_lvl,2)==-1),
      Roma_cw   = as.integer(Roma_lvl   == -1 & lag(Roma_lvl)  ==-1 & lag(Roma_lvl,2)==-1)
    ) %>%
    replace(is.na(.), 0) %>%
    mutate(Week = ISOweek(Time)) %>%
    group_by(Week) %>%
    summarise(
      Attiki_cold_wave3 = sum(Attiki_cw),
      Lisbon_cold_wave3 = sum(Lisbon_cw),
      Roma_cold_wave3   = sum(Roma_cw),
      .groups = "drop"
    )
  
  left_join(hot, cold, by = "Week")
}

# 4. Compute for each scenario
wave_26 <- summarise_wave(rcp26.max, rcp26.min)
wave_85 <- summarise_wave(rcp85.max, rcp85.min)

# 5. Save into a list and write out
wave_list_rcp26 <- list(
  Attiki = wave_26 %>%
    select(Attiki_hot_wave3,
           Attiki_cold_wave3),
  Lisbon = wave_26 %>%
    select(Lisbon_hot_wave3,
           Lisbon_cold_wave3),
  Rome   = wave_26 %>%
    select(Roma_hot_wave3,
           Roma_cold_wave3)
)

wave_list_rcp85 <- list(
  Attiki = wave_85 %>%
    select(Attiki_hot_wave3,
           Attiki_cold_wave3),
  Lisbon = wave_85 %>%
    select(Lisbon_hot_wave3,
           Lisbon_cold_wave3),
  Rome   = wave_85 %>%
    select(Roma_hot_wave3,
           Roma_cold_wave3)
)

save(
  wave_list_rcp26,
  file = "Data/Simulation_data/UTCI/rcp26/wave_list_rcp26.RData"
)

save(
  wave_list_rcp85,
  file = "Data/Simulation_data/UTCI/rcp85/wave_list_rcp85.RData"
)


## lag 21 function
lag21_function <- function(data_df, region_col, cb_matrix, start_day, fill_na = F) {
  # data_df:      daily data frame containing 'Week' and the region column
  # region_col:   string name of the column to process (e.g., "Attiki")
  # cb_matrix:    pre-allocated matrix (rows = number of weeks, cols = 22)
  # start_day:    start day index of the week (1, 4, or 7)
  
  # define extraction order based on start_day
  start <- switch(as.character(start_day),
                  "1" = c(1,7,6,5,4,3,2),
                  "4" = c(4,3,2,1,7,6,5),
                  "7" = c(7,6,5,4,3,2,1),
                  stop("start_day must be 1, 4, or 7")
  )
  
  # group daily values by week into a list-column
  weekly_vals <- data_df %>%
    group_by(Week) %>%
    summarize(vals = list(.data[[region_col]]), .groups = "drop")
  n_weeks <- nrow(weekly_vals)
  if (n_weeks != nrow(cb_matrix)) {
    stop("cb_matrix row count must match number of weeks")
  }
  
  # populate the matrix: current week and lags of 1-2 weeks
  for (i in seq_len(7)) {
    day_vals <- vapply(weekly_vals$vals, `[`, numeric(1), start[i])
    cb_matrix[, i]      <- day_vals
    cb_matrix[, i + 7]  <- dplyr::lag(day_vals, 1)
    cb_matrix[, i + 14] <- dplyr::lag(day_vals, 2)
  }
  # add 3-week lag for the first day only
  cb_matrix[, 22] <- dplyr::lag(vapply(weekly_vals$vals, `[`, numeric(1), start[1]), 3)
  
  # optionally carry nearest non-NA values forward/backward
  if (fill_na) {
    cb_matrix <- apply(cb_matrix, 2, function(col) {
      col <- na.locf(col,      na.rm = FALSE)
      col <- na.locf(col, fromLast = TRUE, na.rm = FALSE)
      col
    })
  }
  
  cb_matrix
}



## Crossbasis matrix (not smoothing)
n_weeks <- nrow(wave_list_rcp26$Attiki)
Attiki_prep_cb_rcp26 <- matrix(NA_real_, nrow = n_weeks, ncol = 22)

Attiki_prep_cb_rcp26 <- lag21_function(
  data_df    = rcp26.mean,
  region_col = "Attiki",
  cb_matrix  = Attiki_prep_cb_rcp26,
  start_day  = 7,
  fill_na    = TRUE
)

Lisbon_prep_cb_rcp26 <- matrix(NA_real_, nrow = n_weeks, ncol = 22)

Lisbon_prep_cb_rcp26 <- lag21_function(
  data_df    = rcp26.mean,
  region_col = "Lisbon",
  cb_matrix  = Lisbon_prep_cb_rcp26,
  start_day  = 7,
  fill_na    = TRUE
)

Roma_prep_cb_rcp26 <- matrix(NA_real_, nrow = n_weeks, ncol = 22)

Roma_prep_cb_rcp26 <- lag21_function(
  data_df    = rcp26.mean,
  region_col = "Roma",
  cb_matrix  = Roma_prep_cb_rcp26,
  start_day  = 7,
  fill_na    = TRUE
)

n_weeks <- nrow(wave_list_rcp85$Attiki)
Attiki_prep_cb_rcp85 <- matrix(NA_real_, nrow = n_weeks, ncol = 22)

Attiki_prep_cb_rcp85 <- lag21_function(
  data_df    = rcp85.mean,
  region_col = "Attiki",
  cb_matrix  = Attiki_prep_cb_rcp85,
  start_day  = 7,
  fill_na    = TRUE
)

Lisbon_prep_cb_rcp85 <- matrix(NA_real_, nrow = n_weeks, ncol = 22)

Lisbon_prep_cb_rcp85 <- lag21_function(
  data_df    = rcp85.mean,
  region_col = "Lisbon",
  cb_matrix  = Lisbon_prep_cb_rcp85,
  start_day  = 7,
  fill_na    = TRUE
)

Roma_prep_cb_rcp85 <- matrix(NA_real_, nrow = n_weeks, ncol = 22)

Roma_prep_cb_rcp85 <- lag21_function(
  data_df    = rcp85.mean,
  region_col = "Roma",
  cb_matrix  = Roma_prep_cb_rcp85,
  start_day  = 7,
  fill_na    = TRUE
)


## Save Crossbasis matrix by RData
save(Attiki_prep_cb_rcp26, file = "Data/Simulation_data/UTCI/rcp26/Attiki_cb_rcp26.RData")
save(Lisbon_prep_cb_rcp26, file = "Data/Simulation_data/UTCI/rcp26/Lisbon_cb_rcp26.RData")
save(Roma_prep_cb_rcp26, file = "Data/Simulation_data/UTCI/rcp26/Roma_cb_rcp26.RData")

save(Attiki_prep_cb_rcp85, file = "Data/Simulation_data/UTCI/rcp85/Attiki_cb_rcp85.RData")
save(Lisbon_prep_cb_rcp85, file = "Data/Simulation_data/UTCI/rcp85/Lisbon_cb_rcp85.RData")
save(Roma_prep_cb_rcp85, file = "Data/Simulation_data/UTCI/rcp85/Roma_cb_rcp85.RData")

