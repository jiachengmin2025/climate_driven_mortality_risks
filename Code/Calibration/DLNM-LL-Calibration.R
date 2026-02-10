# Load packages
packages <- c("readxl", "dlnm", "splines", "tidyr", "ggplot2", "patchwork", "scales")
invisible(lapply(packages, library, character.only = TRUE))

## Load data
region_name = c("Attiki", "Lisbon", "Roma")
age = c('Y20_64', 'Y65_74', 'Y75_84', 'Y_GE85')
setwd("Data/Combined_data/")
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

# before cross-basis
load("Code/Crossbasis_matrix/cb_Attiki.RData")
load("Code/Crossbasis_matrix/cb_Lisbon.RData")
load("Code/Crossbasis_matrix/cb_Roma.RData")

# Attiki
varknots <- equalknots(UTCI_ext.Attiki, fun="ns",df = 3, degree = 3)
lagknots <- logknots(21, fun = "ns", df = 3)
Attiki_cb <- crossbasis(UTCI_ext.Attiki, lag=21, 
                        argvar=list(fun="ns", knots = varknots), 
                        arglag=list(knots=lagknots,df = 3))

# Lisbon
varknots <- equalknots(UTCI_ext.Lisbon,fun="ns",df = 3, degree = 3)
lagknots <- logknots(21, fun = "ns", df = 3)
Lisbon_cb <- crossbasis(UTCI_ext.Lisbon, lag=21, 
                        argvar=list(fun="ns", knots = varknots), 
                        arglag=list(knots=lagknots,df = 3))

# Roma
varknots <- equalknots(UTCI_ext.Lisbon,fun="ns",df = 3, degree = 3)
lagknots <- logknots(21, fun = "ns", df = 3)
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
Time = 1:260 


## Load functions
files <- list.files("Code/Function/", pattern = "\\.R$", full.names = TRUE)
sapply(files, source)


## Fit B-DLNM-LL
## B-DLNM-LL/LC
b = B_DLNM_LL(dat_list = dat_list, wave_list = wave_list, cb_list = cb_list, tol = 1e-2, max_iter = 20)
b1 = B_DLNM_LC(Attiki, Attiki_cb, wave_data = Attiki_wave,  tol = 1e-2, max_iter = 20, region = "Attiki")
b2 = B_DLNM_LC(Lisbon, Lisbon_cb, wave_data = Lisbon_wave, tol = 1e-2, max_iter = 20, region = "Lisbon")
b3 = B_DLNM_LC(Roma, Roma_cb, wave_data = Roma_wave, tol = 1e-2, max_iter = 20, region = "Roma")

## B-DLNM-LL climate-driven mortality rates
LL_Attiki_a =  as.data.frame(matrix(rep(b$final_LL$Ax$Attiki, times = 260), ncol = 4, byrow = TRUE))
LL_Attiki_temp = exp(log(Attiki) - b$final_LL$log_fitted$Attiki) # DLNM part = log(M) - A(x) - B(x)K(t) - b(x)k(t)
colnames(LL_Attiki_temp) = colnames(Attiki)

LL_Lisbon_a =  as.data.frame(matrix(rep(b$final_LL$Ax$Lisbon, times = 260), ncol = 4, byrow = TRUE))
LL_Lisbon_temp = exp(log(Lisbon) - b$final_LL$log_fitted$Lisbon)
colnames(LL_Lisbon_temp) = colnames(Lisbon)

LL_Roma_a =  as.data.frame(matrix(rep(b$final_LL$Ax$Roma, times = 260), ncol = 4, byrow = TRUE))
LL_Roma_temp = exp(log(Roma) - b$final_LL$log_fitted$Roma)
colnames(LL_Roma_temp) = colnames(Roma)

## B-DLNM-LC climate-driven mortality components
Attiki_a = as.data.frame(matrix(rep(b1$final_LC$a_x, times = 260), ncol = 4, byrow = TRUE))
Attiki_temp = exp(log(Attiki) - Attiki_a - b1$final_LC$log_bk_fitted)
names(Attiki_temp) = colnames(Attiki)

Lisbon_a = as.data.frame(matrix(rep(b2$final_LC$a_x, times = 260), ncol = 4, byrow = TRUE))
Lisbon_temp = exp(log(Lisbon) - Lisbon_a - b2$final_LC$log_bk_fitted)
names(Lisbon_temp) = colnames(Lisbon)

Roma_a = as.data.frame(matrix(rep(b3$final_LC$a_x, times = 260), ncol = 4, byrow = TRUE))
Roma_temp = exp(log(Roma) - Roma_a - b3$final_LC$log_bk_fitted)
names(Roma_temp) = colnames(Roma_temp)




# DLNM-LC
# ratio = m_hat/m_tilde
Attiki_ratio = exp(b1$log_fitted - Attiki_a - b1$final_LC$std_mxt)
Lisbon_ratio = exp(b2$log_fitted - Lisbon_a - b2$final_LC$std_mxt)
Roma_ratio = exp(b3$log_fitted - Roma_a - b3$final_LC$std_mxt)
colnames(Attiki_ratio) = c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
colnames(Lisbon_ratio) = c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")
colnames(Roma_ratio) = c("Y20_64", "Y65_74", "Y75_84", "Y_GE85")

# DLNM-LL
#ratio = m_hat/m_tilde
Attiki_ratio1 = exp(b$log_fitted$Attiki - LL_Attiki_a - b$final_LL$std_mxti$Attiki)
Lisbon_ratio1 = exp(b$log_fitted$Lisbon - LL_Lisbon_a - b$final_LL$std_mxti$Lisbon)
Roma_ratio1 = exp(b$log_fitted$Roma - LL_Roma_a - b$final_LL$std_mxti$Roma)


theta_Attiki = 1-1/Attiki_ratio
theta_Attiki1 = 1-1/Attiki_ratio1

theta_Lisbon = 1-1/Lisbon_ratio
theta_Lisbon1 = 1-1/Lisbon_ratio1

theta_Roma = 1-1/Roma_ratio
theta_Roma1 = 1-1/Roma_ratio1

# Define age group labels
age_group_labels <- list(
  "Y20_64" = "Age 20-64",
  "Y65_74" = "Age 65-74",
  "Y75_84" = "Age 75-84",
  "Y_GE85" = "Age 85+"
)

# Find global min/max across all three regions and expand limits
all_values <- unlist(c(theta_Attiki, theta_Attiki1, theta_Lisbon, theta_Lisbon1, theta_Roma, theta_Roma1))
min_val <- min(all_values, na.rm = TRUE) - 0.05 * abs(min(all_values, na.rm = TRUE))  # Expanded lower limit
max_val <- max(all_values, na.rm = TRUE) + 0.05 * abs(max(all_values, na.rm = TRUE))  # Expanded upper limit

# Set step size for 4-5 y-axis ticks
step_size <- (max_val - min_val) / 4  

# Function to generate theta(x,t) plots for different regions and age groups
plot_theta <- function(region_name, age_group, theta_base, theta_final) {
  df <- data.frame(
    time = seq(as.Date("2015-01-01"), as.Date("2019-12-30"), length.out = 260),  # Weekly dates for 5 years
    B_DLNM_LC = as.numeric(theta_base[[age_group]]),
    B_DLNM_LL = as.numeric(theta_final[[age_group]])
  )
  
  df_long <- pivot_longer(df, cols = c("B_DLNM_LC", "B_DLNM_LL"), names_to = "Model", values_to = "Value")
  
  df_long$Model <- factor(df_long$Model, levels = c("B_DLNM_LC", "B_DLNM_LL"))
  
  ggplot(df_long, aes(x = time, y = Value, color = Model)) +
    geom_line(size = 0.3) +  
    scale_color_manual(
      name = "", 
      values = c("B_DLNM_LC" = "black", "B_DLNM_LL" = "#008F7A"),  
      labels = c("DLNM-LC", "DLNM-LL")
    ) +
    scale_x_date(
      limits = c(as.Date("2015-01-01"), as.Date("2020-01-01")),  # Ensure full range
      date_breaks = "1 year",  # One tick per year
      date_labels = "%Y",  # Show year only
      expand = c(0, 0.01)  # Prevent cutoff
    ) +
    scale_y_continuous(limits = c(min_val, max_val),  
                       breaks = seq(min_val, max_val, by = step_size),  
                       labels = number_format(accuracy = 0.01)) +  # Keep 2 decimal places
    labs(
      title = paste0(age_group_labels[[age_group]], " (", region_name, ")"),  # Ensure correct format
      x = "Year",
      y = expression(theta)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 7, face = "plain"),  # Subplot titles smaller, not bold
      panel.grid = element_blank(),  
      panel.border = element_rect(color = "black", fill = NA, size = 1),  
      legend.position = "none"
    )
}

# Generate plots for each region and age group
regions <- list(
  "Athens" = list(theta_Attiki, theta_Attiki1),
  "Lisbon" = list(theta_Lisbon, theta_Lisbon1),
  "Rome" = list(theta_Roma, theta_Roma1)
)

plots <- list()
for (region in names(regions)) {
  theta_base <- regions[[region]][[1]]
  theta_final <- regions[[region]][[2]]
  for (age_group in names(age_group_labels)) {
    plots[[paste(region, age_group)]] <- plot_theta(region, age_group, theta_base, theta_final)
  }
}

# Arrange plots into a 3x4 grid
final_plot <- (plots[["Athens Y20_64"]] | plots[["Athens Y65_74"]] | plots[["Athens Y75_84"]] | plots[["Athens Y_GE85"]]) /
  (plots[["Lisbon Y20_64"]] | plots[["Lisbon Y65_74"]] | plots[["Lisbon Y75_84"]] | plots[["Lisbon Y_GE85"]]) /
  (plots[["Rome Y20_64"]] | plots[["Rome Y65_74"]] | plots[["Rome Y75_84"]] | plots[["Rome Y_GE85"]]) +
  plot_layout(guides = "collect") & 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Ensure main title is bold
    legend.position = "bottom"
  )

# Display the final plot
print(final_plot)
ggsave("Figures/Calibration/fraction.pdf", plot = final_plot, width = 10, height = 6.18)


## Summary statistic of theta
summary(theta_Attiki) #-31~32
summary(theta_Lisbon) # -50~38
summary(theta_Roma) # -37~35

summary(theta_Attiki1) #-34~43
summary(theta_Lisbon1) # -50~40
summary(theta_Roma1) # -39~42

groups <- rep(1:(nrow(theta_Attiki)/52), each = 52)

theta_Attiki_year <- as.data.frame(
  sapply(theta_Attiki, function(col) tapply(col, groups, mean, na.rm = TRUE))
)
theta_Lisbon_year <- as.data.frame(
  sapply(theta_Lisbon, function(col) tapply(col, groups, mean, na.rm = TRUE))
)
theta_Roma_year <- as.data.frame(
  sapply(theta_Roma, function(col) tapply(col, groups, mean, na.rm = TRUE))
)

theta_Attiki_year1 <- as.data.frame(
  sapply(theta_Attiki1, function(col) tapply(col, groups, mean, na.rm = TRUE))
)
theta_Lisbon_year1 <- as.data.frame(
  sapply(theta_Lisbon1, function(col) tapply(col, groups, mean, na.rm = TRUE))
)
theta_Roma_year1 <- as.data.frame(
  sapply(theta_Roma1, function(col) tapply(col, groups, mean, na.rm = TRUE))
)

print(theta_Attiki_year)
print(theta_Lisbon_year)
print(theta_Roma_year)

print(theta_Attiki_year1)
print(theta_Lisbon_year1)
print(theta_Roma_year1)



# Fit DLNM on each region each age group (on the whole dataset) for DLNM-LL
## Overall Effect
## Attiki
DLNM20.Attiki = glm(LL_Attiki_temp$Y20_64 ~ Attiki_cb + Attiki_wave$Attiki_hot_wave3 + Attiki_wave$Attiki_cold_wave3, family = quasipoisson())
DLNM65.Attiki = glm(LL_Attiki_temp$Y65_74   ~ Attiki_cb + Attiki_wave$Attiki_hot_wave3 + Attiki_wave$Attiki_cold_wave3, family = quasipoisson())
DLNM75.Attiki = glm(LL_Attiki_temp$Y75_84   ~ Attiki_cb + Attiki_wave$Attiki_hot_wave3 + Attiki_wave$Attiki_cold_wave3, family = quasipoisson())
DLNM85.Attiki = glm(LL_Attiki_temp$Y_GE85  ~ Attiki_cb + Attiki_wave$Attiki_hot_wave3 + Attiki_wave$Attiki_cold_wave3, family = quasipoisson())

pdf("Figures/DLNM/Overall/Athens.pdf", width = 10, height = 6.18)
par(
  mfrow   = c(2, 2),
  mar     = c(5, 5, 2, 1),
  oma     = c(1, 1, 3, 1),
  cex.lab = 2,
  cex.axis= 2, 
  cex.main= 2
)

pred20.Attiki <- crosspred(Attiki_cb, DLNM20.Attiki, cen = 22, at = -5:40)
plot(
  pred20.Attiki, "overall",
  xlab   = "UTCI",
  ylab   = "RR",
  xlim   = c(-5, 40),
  ylim   = c(0.9, 1.8),
  ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
  main   = "Age 20-64"
)

pred65.Attiki <- crosspred(Attiki_cb, DLNM65.Attiki, cen = 22, at = -5:40)
plot(
  pred65.Attiki, "overall",
  xlab   = "UTCI",
  ylab   = "RR",
  xlim   = c(-5, 40),
  ylim   = c(0.9, 1.8),
  ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
  main   = "Age 65-74"
)

pred75.Attiki <- crosspred(Attiki_cb, DLNM75.Attiki, cen = 20, at = -5:40)
plot(
  pred75.Attiki, "overall",
  xlab   = "UTCI",
  ylab   = "RR",
  xlim   = c(-5, 40),
  ylim   = c(0.9, 1.8),
  ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
  main   = "Age 75-84"
)

pred85.Attiki <- crosspred(Attiki_cb, DLNM85.Attiki, cen = 20, at = -5:40)
plot(
  pred85.Attiki, "overall",
  xlab   = "UTCI",
  ylab   = "RR",
  xlim   = c(-5, 40),
  ylim   = c(0.9, 1.8),
  ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
  main   = "Age 85+"
)
dev.off()


## Lisbon
DLNM20.Lisbon = glm(LL_Lisbon_temp$Y20_64 ~ Lisbon_cb + Lisbon_wave$Lisbon_hot_wave3 + Lisbon_wave$Lisbon_cold_wave3, family = quasipoisson())
DLNM65.Lisbon = glm(LL_Lisbon_temp$Y65_74  ~ Lisbon_cb+ Lisbon_wave$Lisbon_hot_wave3 + Lisbon_wave$Lisbon_cold_wave3, family = quasipoisson())
DLNM75.Lisbon = glm(LL_Lisbon_temp$Y75_84  ~ Lisbon_cb+ Lisbon_wave$Lisbon_hot_wave3 + Lisbon_wave$Lisbon_cold_wave3, family = quasipoisson())
DLNM85.Lisbon = glm(LL_Lisbon_temp$Y_GE85  ~ Lisbon_cb+ Lisbon_wave$Lisbon_hot_wave3 + Lisbon_wave$Lisbon_cold_wave3, family = quasipoisson())



pdf("Figures/DLNM/Overall/Lisbon.pdf", width = 10, height = 6.18)
par(
  mfrow   = c(2, 2),
  mar     = c(5, 5, 2, 1),
  oma     = c(1, 1, 3, 1),
  cex.lab = 2,
  cex.axis= 2, 
  cex.main= 2
)
pred20.Lisbon = crosspred(Lisbon_cb, DLNM20.Lisbon, cen = 22, at = -5:40)
plot(pred20.Lisbon, "overall", xlab = "UTCI", ylab = "RR", xlim = c(-5, 40), ylim = c(0.9, 1.8),
     ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
     main = "Age 20-64")

pred65.Lisbon = crosspred(Lisbon_cb, DLNM65.Lisbon, cen = 22, at = -5:40)
plot(pred65.Lisbon, "overall", xlab = "UTCI", ylab = "RR", xlim = c(-5, 40), ylim = c(0.9, 1.8),
     ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
     main = "Age 65-74")

pred75.Lisbon = crosspred(Lisbon_cb, DLNM75.Lisbon, cen = 20, at = -5:40)
plot(pred75.Lisbon, "overall", xlab = "UTCI", ylab = "RR", xlim = c(-5, 40), ylim = c(0.9, 1.8),
     ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
     main = "Age 75-84")

pred85.Lisbon = crosspred(Lisbon_cb, DLNM85.Lisbon, cen = 20, at = -5:40)
plot(pred85.Lisbon, "overall", xlab = "UTCI", ylab = "RR", xlim = c(-5, 40), ylim = c(0.9, 3.5),
     ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
     main = "Age 85+")
dev.off()


## Roma
DLNM20.Roma = glm(LL_Roma_temp$Y20_64 ~ Roma_cb + Roma_wave$Roma_hot_wave3 + Roma_wave$Roma_cold_wave3, family = quasipoisson())
DLNM65.Roma = glm(LL_Roma_temp$Y65_74  ~ Roma_cb + Roma_wave$Roma_hot_wave3 + Roma_wave$Roma_cold_wave3, family = quasipoisson())
DLNM75.Roma = glm(LL_Roma_temp$Y75_84  ~ Roma_cb + Roma_wave$Roma_hot_wave3 + Roma_wave$Roma_cold_wave3, family = quasipoisson())
DLNM85.Roma = glm(LL_Roma_temp$Y_GE85  ~ Roma_cb + Roma_wave$Roma_hot_wave3 + Roma_wave$Roma_cold_wave3, family = quasipoisson())

pdf("Figures/DLNM/Overall/Rome.pdf", width = 10, height = 6.18)
par(
  mfrow   = c(2, 2),
  mar     = c(5, 5, 2, 1),
  oma     = c(1, 1, 3, 1),
  cex.lab = 2,
  cex.axis= 2, 
  cex.main= 2
)
pred20.Roma = crosspred(Roma_cb, DLNM20.Roma, cen = 18, at = -5:40)
plot(pred20.Roma, "overall", xlab = "UTCI", ylab = "RR", xlim = c(-5, 40), ylim = c(0.9, 1.8),
     ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
     main = "Age 20-64")

pred65.Roma = crosspred(Roma_cb, DLNM65.Roma, cen = 20, at = -5:40)
plot(pred65.Roma, "overall", xlab = "UTCI", ylab = "RR", xlim = c(-5, 40), ylim = c(0.9, 1.8),
     ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
     main = "Age 65-74")

pred75.Roma = crosspred(Roma_cb, DLNM75.Roma, cen = 20, at = -5:40)
plot(pred75.Roma, "overall", xlab = "UTCI", ylab = "RR", xlim = c(-5, 40), ylim = c(0.9, 1.8),
     ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
     main = "Age 75-84")

pred85.Roma = crosspred(Roma_cb, DLNM85.Roma, cen = 20, at = -5:40)
plot(pred85.Roma, "overall", xlab = "UTCI", ylab = "RR", xlim = c(-5, 40), ylim = c(0.9, 1.8),
     ci.arg = list(lty = 2, lwd = 2, col = "#FFBE7A"),
     main = "Age 85+")
dev.off()


## Lag Effect
## Attiki
pdf("Figures/DLNM/Lag/Athens_lag.pdf", width = 16, height = 8)

par(
  mfcol   = c(2, 4),
  mar     = c(5, 5, 2, 1),
  oma     = c(1, 1, 3, 1),
  ps      = 16,
  cex.lab = 1.2,
  cex.axis= 1.2,
  cex.main= 1.4,
  font.main = 2
)

plot(pred20.Attiki, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 20-64", ylim = c(0.98, 1.05))

plot(pred20.Attiki, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 20-64", ylim = c(0.98, 1.05))

plot(pred65.Attiki, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 65-74", ylim = c(0.98, 1.05))

plot(pred65.Attiki, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 65-74", ylim = c(0.98, 1.05))

plot(pred75.Attiki, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 75-84", ylim = c(0.98, 1.05))

plot(pred75.Attiki, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 75-84", ylim = c(0.98, 1.05))

plot(pred85.Attiki, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 85+", ylim = c(0.98, 1.05))

plot(pred85.Attiki, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 85+", ylim = c(0.98, 1.05))
dev.off()


## Lisbon
pdf("Figures/DLNM/Lag/Lisbon_lag.pdf", width = 16, height = 8)

par(
  mfcol   = c(2, 4),
  mar     = c(5, 5, 2, 1),
  oma     = c(1, 1, 3, 1),
  ps      = 16,
  cex.lab = 1.2,
  cex.axis= 1.2,
  cex.main= 1.4,
  font.main = 2
)

plot(pred20.Lisbon, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 20-64", ylim = c(0.95, 1.2))

plot(pred20.Lisbon, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 20-64", ylim = c(0.95, 1.2))

plot(pred65.Lisbon, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 65-74", ylim = c(0.95, 1.2))

plot(pred65.Lisbon, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 65-74", ylim = c(0.95, 1.2))

plot(pred75.Lisbon, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 75-84", ylim = c(0.95, 1.2))

plot(pred75.Lisbon, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 75-84", ylim = c(0.95, 1.2))

plot(pred85.Lisbon, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 85+", ylim = c(0.95, 1.2))

plot(pred85.Lisbon, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 85+", ylim = c(0.95, 1.2))
dev.off()


## Roma
pdf("Figures/DLNM/Lag/Roma_lag.pdf", width = 16, height = 8)

par(
  mfcol   = c(2, 4),
  mar     = c(5, 5, 2, 1),
  oma     = c(1, 1, 3, 1),
  ps      = 16,
  cex.lab = 1.2,
  cex.axis= 1.2,
  cex.main= 1.4,
  font.main = 2
)

plot(pred20.Roma, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 20-64", ylim = c(0.95, 1.1))

plot(pred20.Roma, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 20-64", ylim = c(0.95, 1.1))

plot(pred65.Roma, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 65-74", ylim = c(0.95, 1.1))

plot(pred65.Roma, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 65-74", ylim = c(0.95, 1.1))

plot(pred75.Roma, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 75-84", ylim = c(0.95, 1.1))

plot(pred75.Roma, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 75-84", ylim = c(0.95, 1.1))

plot(pred85.Roma, "slices", var = -5, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = -5, Age 85+", ylim = c(0.95, 1.1))

plot(pred85.Roma, "slices", var = 35, col = 1, ylab = "RR", 
     ci.arg = list(density = 50, lwd = 2, col = "#AFEEEE"),
     main = "UTCI = 35, Age 85+", ylim = c(0.95, 1.1))
dev.off()


## Coefficients - LL
round(unname(c(coef(DLNM20.Attiki)[11:12], coef(DLNM20.Lisbon)[11:12], coef(DLNM20.Roma)[11:12])), 4)
round(unname(c(coef(DLNM65.Attiki)[11:12], coef(DLNM65.Lisbon)[11:12], coef(DLNM65.Roma)[11:12])), 4)
round(unname(c(coef(DLNM75.Attiki)[11:12], coef(DLNM75.Lisbon)[11:12], coef(DLNM75.Roma)[11:12])), 4)
round(unname(c(coef(DLNM85.Attiki)[11:12], coef(DLNM85.Lisbon)[11:12], coef(DLNM85.Roma)[11:12])), 4)

summary(DLNM65.Roma)


# Save single DLNM coefficeints
# Regions and age‐groups
regions   <- c("Attiki","Lisbon","Roma")
age_groups <- c("Y20_64","Y65_74","Y75_84","Y_GE85")

# Initialize the nested list
DLNM_coef_list <- vector("list", length(regions))
names(DLNM_coef_list) <- regions

# Loop over regions
for(region in regions) {
  # Retrieve the temp data, cb matrix name, and wave data for this region
  temp_data   <- get(paste0("LL_", region, "_temp"))
  cb_vector   <- get(paste0(region, "_cb"))
  wave_data   <- get(paste0(region, "_wave"))
  
  # Fit one model per age‐group
  models <- lapply(age_groups, function(age) {
    glm(
      temp_data[[age]] ~ 
        cb_vector +
        wave_data[[paste0(region, "_hot_wave3")]] +
        wave_data[[paste0(region, "_cold_wave3")]],
      family = quasipoisson()
    )
  })
  names(models) <- age_groups
  
  # Extract coefficients and vcov
  coefs <- lapply(models, coef)
  vcovs  <- lapply(models, vcov)
  
  # Store in the master list
  DLNM_coef_list[[region]] <- list(
    coef = coefs,
    vcov = vcovs
  )
}
