## Load packages
packages <- c("readxl", "dlnm", "ggplot2", "patchwork", "splines")
invisible(lapply(packages, library, character.only = TRUE))

## Load data
region_name = c("Attiki", "Lisbon", "Roma")
age = c('Y20_64', 'Y65_74', 'Y75_84', 'Y_GE85')
Y20_64 = read_xlsx("Data/Combined_data/Y20_64_combined.xlsx")
Y65_74 = read_xlsx("Data/Combined_data/Y65_74_combined.xlsx")
Y75_84 = read_xlsx("Data/Combined_data/Y75_84_combined.xlsx")
Y_GE85 = read_xlsx("Data/Combined_data/Y_GE85_combined.xlsx")
 
## LC matrix (mortality rate)

col_name = c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")

Attiki = cbind(Y20_64[,3]/Y20_64[,9],Y65_74[,3]/Y65_74[,9], 
               Y75_84[,3]/Y75_84[,9], Y_GE85[,3]/Y_GE85[,9])
colnames(Attiki) = col_name

Lisbon = cbind(Y20_64[,4]/Y20_64[,10],Y65_74[,4]/Y65_74[,10],
               Y75_84[,4]/Y75_84[,10], Y_GE85[,4]/Y_GE85[,10])
colnames(Lisbon) = col_name

Roma = cbind(Y20_64[,6]/Y20_64[,12],Y65_74[,6]/Y65_74[,12], 
             Y75_84[,6]/Y75_84[,12], Y_GE85[,6]/Y_GE85[,12])
colnames(Roma) = col_name


## Store heat/cold wave data
# Attiki
Attiki_wave = data.frame(cbind(Y20_64$Attiki_hot_wave3, Y20_64$Attiki_cold_wave3))
colnames(Attiki_wave) = c("Attiki_hot_wave3", "Attiki_cold_wave3")
# Lisbon
Lisbon_wave = data.frame(cbind(Y20_64$Lisbon_hot_wave3, Y20_64$Lisbon_cold_wave3))
colnames(Lisbon_wave) = c("Lisbon_hot_wave3", "Lisbon_cold_wave3")
# Roma
Roma_wave = data.frame(cbind(Y20_64$Roma_hot_wave3, Y20_64$Roma_cold_wave3))
colnames(Roma_wave) = c("Roma_hot_wave3", "Roma_cold_wave3")


# DLNM 
## Crossbasis matrix
load("Code/Crossbasis_matrix/cb_Attiki.RData")
load("Code/Crossbasis_matrix/cb_Lisbon.RData")
load("Code/Crossbasis_matrix/cb_Roma.RData")


# Attiki
varknots <- equalknots(UTCI_ext.Attiki, fun="ns",df = 3, degree = 3)
lagknots <- logknots(21,3)
Attiki_cb <- crossbasis(UTCI_ext.Attiki, lag=21, 
                        argvar=list(fun="ns", knots = varknots), 
                        arglag=list(knots=lagknots,df = 3))

# Lisbon
varknots <- equalknots(UTCI_ext.Lisbon,fun="ns",df = 3, degree = 3)
lagknots <- logknots(21, 3)
Lisbon_cb <- crossbasis(UTCI_ext.Lisbon, lag=21, 
                        argvar=list(fun="ns", knots = varknots), 
                        arglag=list(knots=lagknots,df = 3))

# Roma
varknots <- equalknots(UTCI_ext.Roma,fun="ns",df = 3, degree = 3)
lagknots <- logknots(21, 3)
Roma_cb <- crossbasis(UTCI_ext.Roma, lag=21, 
                      argvar=list(fun="ns", knots = varknots), 
                      arglag=list(knots=lagknots,df = 3))



## Convert to standard weekly mortality (*52)
Attiki = Attiki * 52
Lisbon = Lisbon * 52
Roma = Roma * 52
dat_list <- list(Attiki = Attiki, Lisbon = Lisbon, Roma = Roma)
wave_list <- list(Attiki = Attiki_wave, Lisbon = Lisbon_wave, Roma = Roma_wave)
cb_list <- list(Attiki = Attiki_cb, Lisbon = Lisbon_cb, Roma = Roma_cb)


# Load functions
files <- list.files("Code/Function/", pattern = "\\.R$", full.names = TRUE)
sapply(files, source)


## LC and B-DLNM-LC

# Fit B-DLNM-LC model
BDLNMLC_Attiki = B_DLNM_LC(Attiki, Attiki_cb, Attiki_wave, 1e-2, 20, "Attiki")
BDLNMLC_Lisbon = B_DLNM_LC(Lisbon, Lisbon_cb, Lisbon_wave, 1e-2, 20, "Lisbon")
BDLNMLC_Roma = B_DLNM_LC(Roma, Roma_cb, Roma_wave, 1e-2, 20, "Roma")


# kappa(t) plotting
# Function to generate kappa(t) plots for LC or DLNM-LC with shared y-scale
plot_kappa_separate <- function(region_name, data, model_type, y_limits) {
  df <- data.frame(
    time = seq(as.Date("2015-01-01"), as.Date("2019-12-30"), length.out = 260),
    kappa_t = if (model_type == "LC") as.numeric(data$baseline_LC$k_t) else as.numeric(data$final_LC$k_t),
    Model = model_type
  )
  
  ggplot(df, aes(x = time, y = kappa_t, color = Model)) +
    geom_line(size = 0.5, linetype = "solid") +
    scale_color_manual(
      name = "",
      values = if (model_type == "LC") c("LC" = "black") else c("DLNM-LC" = "darkred"),
      labels = if (model_type == "LC") c("LC") else c("DLNM-LC")
    ) +
    scale_x_date(
      limits = c(as.Date("2015-01-01"), as.Date("2020-01-01")),
      date_breaks = "1 year",
      date_labels = "%Y",
      expand = c(0, 0.01)
    ) +
    scale_y_continuous(limits = y_limits) +  # **Ensuring shared y-scale**
    labs(
      title = bquote(bold(.(region_name))),  # **Make region name bold**
      x = "Year",
      y = expression(hat(kappa)(t))
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title    = element_text(size = 16),
      axis.text     = element_text(size = 14),
      legend.text   = element_text(size = 14),
      panel.grid    = element_blank(),
      panel.border  = element_rect(color = "black", fill = NA, size = 1),
      legend.position = "none"
    )
}

# Determine global y-axis limits across all regions
all_kappa_values <- c(
  unlist(BDLNMLC_Attiki$baseline_LC$k_t), unlist(BDLNMLC_Attiki$final_LC$k_t),
  unlist(BDLNMLC_Lisbon$baseline_LC$k_t), unlist(BDLNMLC_Lisbon$final_LC$k_t),
  unlist(BDLNMLC_Roma$baseline_LC$k_t), unlist(BDLNMLC_Roma$final_LC$k_t)
)
y_limits <- range(all_kappa_values, na.rm = TRUE)  # Get min & max values for y-axis

# Generate LC model plots
plot_Athens_LC <- plot_kappa_separate("Athens", BDLNMLC_Attiki, "LC", y_limits)
plot_lisbon_LC <- plot_kappa_separate("Lisbon", BDLNMLC_Lisbon, "LC", y_limits)
plot_roma_LC <- plot_kappa_separate("Rome", BDLNMLC_Roma, "LC", y_limits)

# Generate DLNM-LC model plots
plot_Athens_DLNM <- plot_kappa_separate("Athens", BDLNMLC_Attiki, "DLNM-LC", y_limits)
plot_lisbon_DLNM <- plot_kappa_separate("Lisbon", BDLNMLC_Lisbon, "DLNM-LC", y_limits)
plot_roma_DLNM <- plot_kappa_separate("Rome", BDLNMLC_Roma, "DLNM-LC", y_limits)

# Arrange LC and DLNM-LC plots separately
LC_plots <- (plot_Athens_LC / plot_lisbon_LC / plot_roma_LC) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

DLNM_LC_plots <- (plot_Athens_DLNM / plot_lisbon_DLNM / plot_roma_DLNM) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Combine LC and DLNM-LC plots into a two-column layout with shared y-scale
combined_plot <- (LC_plots | DLNM_LC_plots) +
  plot_annotation(
    theme = theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title  = element_text(size = 16),
      axis.text   = element_text(size = 14),
      legend.text = element_text(size = 14)
    )
  ) &
  theme(
    plot.margin = unit(c(5, 20, 5, 5), "pt")
  )

# Save the final figure
ggsave("Figures/Calibration/single_model_k.pdf",
       combined_plot, width = 12, height = 6, dpi = 300)


## LL and B-DLNM-LL
BDLNMLL = B_DLNM_LL(dat_list = dat_list, cb_list = cb_list, wave_list = wave_list, tol = 1e-2, max_iter = 20)

K_final = BDLNMLL$final_LL$Kt
K_base = BDLNMLL$baseline_LL$Kt

k_Attiki_final = BDLNMLL$final_LL$kt$Attiki
k_Attiki_base = BDLNMLL$baseline_LL$kt$Attiki
k_Lisbon_final = BDLNMLL$final_LL$kt$Lisbon
k_Lisbon_base = BDLNMLL$baseline_LL$kt$Lisbon
k_Roma_final = BDLNMLL$final_LL$kt$Roma
k_Roma_base = BDLNMLL$baseline_LL$kt$Roma

# Function to generate kappa(t) plots for LL and DLNM-LL with shared y-scale
plot_kappa_BDLNMLL_separate <- function(region_name, base_kt, final_kt, model_type, y_limits) {
  df <- data.frame(
    time = seq(as.Date("2015-01-01"), as.Date("2019-12-30"), length.out = 260),
    kappa_t = if (model_type == "LL") as.numeric(base_kt) else as.numeric(final_kt),
    Model = model_type
  )
  
  ggplot(df, aes(x = time, y = kappa_t, color = Model)) +
    geom_line(size = 0.5, linetype = "solid") +
    scale_color_manual(
      name = "",
      values = if (model_type == "LL") c("LL" = "black") else c("DLNM-LL" = "darkred"),
      labels = if (model_type == "LL") c("LL") else c("DLNM-LL")
    ) +
    scale_x_date(
      limits = c(as.Date("2015-01-01"), as.Date("2020-01-01")),
      date_breaks = "1 year",
      date_labels = "%Y",
      expand = c(0, 0.01)
    ) +
    scale_y_continuous(limits = y_limits) +  # **Ensuring shared y-scale**
    labs(
      title = bquote(bold(.(region_name))),  # **Make region name bold**
      x = "Year",
      y = expression(hat(kappa)(t))
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title    = element_text(size = 16),
      axis.text     = element_text(size = 14),
      legend.text   = element_text(size = 14),
      panel.grid    = element_blank(),
      panel.border  = element_rect(color = "black", fill = NA, size = 1),
      legend.position = "none"
    )
}

# Determine global y-axis limits across all regions
all_kappa_values <- c(
  unlist(K_base), unlist(K_final),
  unlist(k_Attiki_base), unlist(k_Attiki_final),
  unlist(k_Lisbon_base), unlist(k_Lisbon_final),
  unlist(k_Roma_base), unlist(k_Roma_final)
)
y_limits <- range(all_kappa_values, na.rm = TRUE)  # Get min & max values for y-axis

# Generate LL model plots
plot_common_LL <- plot_kappa_BDLNMLL_separate("Common Trend", K_base, K_final, "LL", y_limits)
plot_attiki_LL <- plot_kappa_BDLNMLL_separate("Athens", k_Attiki_base, k_Attiki_final, "LL", y_limits)
plot_lisbon_LL <- plot_kappa_BDLNMLL_separate("Lisbon", k_Lisbon_base, k_Lisbon_final, "LL", y_limits)
plot_roma_LL <- plot_kappa_BDLNMLL_separate("Rome", k_Roma_base, k_Roma_final, "LL", y_limits)

# Generate DLNM-LL model plots
plot_common_DLNM_LL <- plot_kappa_BDLNMLL_separate("Common Trend", K_base, K_final, "DLNM-LL", y_limits)
plot_attiki_DLNM_LL <- plot_kappa_BDLNMLL_separate("Athens", k_Attiki_base, k_Attiki_final, "DLNM-LL", y_limits)
plot_lisbon_DLNM_LL <- plot_kappa_BDLNMLL_separate("Lisbon", k_Lisbon_base, k_Lisbon_final, "DLNM-LL", y_limits)
plot_roma_DLNM_LL <- plot_kappa_BDLNMLL_separate("Rome", k_Roma_base, k_Roma_final, "DLNM-LL", y_limits)

# Arrange LL and DLNM-LL plots separately
LL_plots <- (plot_common_LL / plot_attiki_LL / plot_lisbon_LL / plot_roma_LL) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

DLNM_LL_plots <- (plot_common_DLNM_LL / plot_attiki_DLNM_LL / plot_lisbon_DLNM_LL / plot_roma_DLNM_LL) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Combine LL and DLNM-LL plots into a two-column layout with shared y-scale
combined_plot <- (LL_plots | DLNM_LL_plots) +
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ) &
  theme(
    plot.margin = unit(c(5, 20, 5, 5), "pt")
  )

# Save the final figure
ggsave("Figures/Calibration/multi_model_k.pdf",
       combined_plot, width = 12, height = 8, dpi = 300)

