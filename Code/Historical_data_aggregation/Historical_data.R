## Load packages
packages <- c("readxl", "ggplot2", "ggpubr", "reshape2")
invisible(lapply(packages, library, character.only = TRUE))

## LC matrix (mortality rate)
Y20_64 = read_xlsx("Data/Combined_data/Y20_64_combined.xlsx")
Y65_74 = read_xlsx("Data/Combined_data/Y65_74_combined.xlsx")
Y75_84 = read_xlsx("Data/Combined_data/Y75_84_combined.xlsx")
Y_GE85 = read_xlsx("Data/Combined_data/Y_GE85_combined.xlsx")
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


## Convert to standard weekly mortality (*52)
Attiki = Attiki * 52
Lisbon = Lisbon * 52
Roma = Roma * 52
dat_list = list(Attiki = Attiki, Lisbon = Lisbon, Roma = Roma)


## Heat/cold wave data
# Attiki
Attiki_wave = data.frame(cbind(Y20_64$Attiki_hot_wave3, Y20_64$Attiki_cold_wave3))
colnames(Attiki_wave) = c("Attiki_hot_wave3", "Attiki_cold_wave3")
# Lisbon
Lisbon_wave = data.frame(cbind(Y20_64$Lisbon_hot_wave3, Y20_64$Lisbon_cold_wave3))
colnames(Lisbon_wave) = c("Lisbon_hot_wave3", "Lisbon_cold_wave3")
# Roma
Roma_wave = data.frame(cbind(Y20_64$Roma_hot_wave3, Y20_64$Roma_cold_wave3))
colnames(Roma_wave) = c("Roma_hot_wave3", "Roma_cold_wave3")

wave_list = list(Attiki = Attiki_wave, Lisbon = Lisbon_wave, Roma = Roma_wave)


## Save data
save(dat_list, wave_list,
     file = "Code/Historical_data_aggregation/Historical_data.RData")



## Mortality Rate Plot

age_group_labels <- c(
  "Y20_64" = "20-64",
  "Y65_74" = "65-74",
  "Y75_84" = "75-84",
  "Y_GE85" = "85+"
)

line_colors <- c("Y20_64" = "#9ecae1",  
                 "Y65_74" = "#6baed6", 
                 "Y75_84" = "#2171b5", 
                 "Y_GE85" = "#084594")
smooth_color <- "red"  


## Attiki
Attiki$Time <- seq(from = as.Date("2015-01-01"), by = "week", length.out = nrow(Attiki))

Attiki_long <- reshape2::melt(Attiki, id.vars="Time",
                              variable.name="AgeGroup", value.name="value")

# Create the plot with log scale
Attiki_plot <- ggplot(Attiki_long) +
  geom_line(aes(x = Time, y = value, color = AgeGroup), size = 1.2, alpha = 0.8) +
  geom_smooth(aes(x = Time, y = value, group = AgeGroup), 
              method = "loess", span = 0.05, size = 0.5, se = TRUE, 
              color = smooth_color) +
  labs(
    x = "Year",
    y = "Mortality",
    color = "Age Group"
  ) +
  scale_color_manual(
    name = "Age Group",
    values = line_colors, 
    labels = age_group_labels
  ) +
  scale_x_date(
    limits = c(as.Date("2015-01-01"), as.Date("2020-01-01")),  
    date_breaks = "1 year",  
    date_labels = "%y",    
    expand = c(0, 0.01)   
  ) +
  scale_y_log10() +  
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_line(color = "gray80", size = 0.5),      
    panel.grid.minor = element_line(color = "gray90", size = 0.25),  
    legend.position = "right",                                        
    plot.title = element_text(hjust = 0.5),                         
    text = element_text(size = 15)                                 
  )


## Lisbon
Lisbon$Time <- seq(from = as.Date("2015-01-01"), by = "week", length.out = nrow(Lisbon))

Lisbon_long <- reshape2::melt(Lisbon, id.vars="Time",
                              variable.name="AgeGroup", value.name="value")

Lisbon_plot <- ggplot(Lisbon_long) +
  geom_line(aes(x = Time, y = value, color = AgeGroup), size = 1.2, alpha = 0.8) +
  geom_smooth(aes(x = Time, y = value, group = AgeGroup), 
              method = "loess", span = 0.03, size = 0.5, se = TRUE, 
              color = smooth_color) +
  labs(
    x = "Year",
    y = "Mortality",
    color = "Age Group"
  ) +
  scale_color_manual(
    name = "Age Group",
    values = line_colors,  
    labels = age_group_labels
  ) +
  scale_x_date(
    limits = c(as.Date("2015-01-01"), as.Date("2020-01-01")),  
    date_breaks = "1 year", 
    date_labels = "%y",      
    expand = c(0, 0.01)      
  ) +
  scale_y_log10() +  
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_line(color = "gray80", size = 0.5),      
    panel.grid.minor = element_line(color = "gray90", size = 0.25),     
    legend.position = "right",                                      
    plot.title = element_text(hjust = 0.5),                            
    text = element_text(size = 15)                           
  )


## Roma
Roma$Time <- seq(from = as.Date("2015-01-01"), by = "week", length.out = nrow(Roma))

Roma_long <- reshape2::melt(Roma, id.vars="Time",
                              variable.name="AgeGroup", value.name="value")

Roma_plot <- ggplot(Roma_long) +
  geom_line(aes(x = Time, y = value, color = AgeGroup), size = 1.2, alpha = 0.8) +
  geom_smooth(aes(x = Time, y = value, group = AgeGroup), 
              method = "loess", span = 0.03, size = 0.5, se = TRUE, 
              color = smooth_color) +
  labs(
    x = "Year",
    y = "Mortality",
    color = "Age Group"
  ) +
  scale_color_manual(
    name = "Age Group",
    values = line_colors, 
    labels = age_group_labels
  ) +
  scale_x_date(
    breaks      = as.Date(c("2015-01-01","2016-01-01",
                            "2017-01-01","2018-01-01")),
    date_labels = "%y",                      
    limits      = c(as.Date("2015-01-01"),
                    as.Date("2018-12-31")),  
    expand      = c(0, 0.01)
  ) +
  scale_y_log10() + 
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_line(color = "gray80", size = 0.5),      
    panel.grid.minor = element_line(color = "gray90", size = 0.25),   
    legend.position = "right",                                        
    plot.title = element_text(hjust = 0.5),                           
    text = element_text(size = 15)                               
  )



# Plot Historical Mortality Rates
# Ensure Date column is correctly formatted
Attiki$Time <- seq(from = as.Date("2015-01-01"), by = "week", length.out = nrow(Attiki))
Lisbon$Time <- seq(from = as.Date("2015-01-01"), by = "week", length.out = nrow(Lisbon))
Roma$Time <- seq(from = as.Date("2015-01-01"), by = "week", length.out = nrow(Roma))

# Convert to long format
Attiki_long <- melt(Attiki, id.vars = "Time", variable.name = "AgeGroup", value.name = "value")
Lisbon_long <- melt(Lisbon, id.vars = "Time", variable.name = "AgeGroup", value.name = "value")
Roma_long <- melt(Roma, id.vars = "Time", variable.name = "AgeGroup", value.name = "value")

y_min <- min(c(Attiki_long$value, Lisbon_long$value, Roma_long$value), na.rm = TRUE)
y_max <- max(c(Attiki_long$value, Lisbon_long$value, Roma_long$value), na.rm = TRUE)
year_lines <- seq(as.Date("2016-01-01"), as.Date("2019-01-01"), by = "year")

plot_mortality <- function(data, region) {
  ggplot(data, aes(x = Time)) +
    geom_vline(xintercept = year_lines,
               linetype   = "dashed",
               color      = "gray50",
               size       = 0.3) +
    geom_line(aes(y = value, color = AgeGroup), size = 1.2, alpha = 0.8) + 
    geom_smooth(aes(y = value, group = AgeGroup),
                method = "loess", span = 0.05, size = 0.5, se = TRUE,
                color = smooth_color) +
    labs(title = region,
         x     = "Year", 
         y     = "Mortality Rate",
         color = "Age Group") +
    scale_color_manual(name   = "Age Group",
                       values = line_colors,
                       labels = age_group_labels) +
    scale_x_date(
      breaks      = seq(as.Date("2015-01-01"), as.Date("2020-01-01"), by = "1 year"),
      date_labels = "%y",                           
      limits      = c(as.Date("2015-01-01"), 
                      as.Date("2020-01-01")),
      expand      = c(0, 0.01)
    ) +
    scale_y_log10(limits = c(y_min, y_max)) +
    theme_minimal() +
    theme(
      legend.position  = "none",
      panel.border     = element_rect(color = "black", fill = NA, size = 1),
      panel.grid       = element_blank(),
      plot.title       = element_text(hjust = 0.5),
      text             = element_text(size = 12)
    )
}

legend_df <- data.frame(
  AgeGroup = factor(names(line_colors), levels = names(line_colors)),
  y        = 1
)

legend_plot <- ggplot(legend_df, aes(x = AgeGroup, y = y, color = AgeGroup)) +
  geom_point() +
  scale_color_manual(
    name   = "Age Group",
    values = line_colors,
    labels = age_group_labels
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 10)
  )

shared_leg <- get_legend(legend_plot)

plot_Athens <- plot_mortality(Attiki_long, "Athens")
plot_Lisbon <- plot_mortality(Lisbon_long, "Lisbon")
plot_Rome   <- plot_mortality(Roma_long,   "Rome")

pa <- plot_Athens + theme(plot.margin = margin(2,5,2,2, "mm"))
pl <- plot_Lisbon + theme(plot.margin = margin(2,5,2,2, "mm"))
pr <- plot_Rome   + theme(plot.margin = margin(2,5,2,2, "mm"))

Mortality_plot <- ggarrange(
  pa, pl, pr,
  ncol          = 3,
  common.legend = TRUE,
  legend        = "bottom"
)

ggsave(
  "Figures/Mortality/Historical_combined.pdf",
  Mortality_plot,
  width  = 7.5,
  height = 3,
  units  = "in",
  dpi    = 300
)



## Load UTCI data

UTCI.mean = read_xlsx("Data/UTCI_data/Daily_data/UTCI_daily_mean.xlsx")
UTCI.max = read_xlsx("Data/UTCI_data/Daily_data/UTCI_daily_max.xlsx")
UTCI.min = read_xlsx("Data/UTCI_data/Daily_data/UTCI_daily_min.xlsx")

UTCI_Attiki = data.frame(Date = UTCI.mean$Date, Mean = UTCI.mean$Attiki, 
                         Min = UTCI.min$Attiki, Max = UTCI.max$Attiki)
UTCI_Lisbon = data.frame(Date = UTCI.mean$Date, Mean = UTCI.mean$Lisbon, 
                         Min = UTCI.min$Lisbon, Max = UTCI.max$Lisbon)
UTCI_Roma = data.frame(Date = UTCI.mean$Date, Mean = UTCI.mean$Roma, 
                       Min = UTCI.min$Roma, Max = UTCI.max$Roma)


## Plot UTCI data

# Ensure Date column is correctly formatted
UTCI_Attiki$Date <- as.Date(UTCI_Attiki$Date)
UTCI_Lisbon$Date <- as.Date(UTCI_Lisbon$Date)
UTCI_Roma$Date <- as.Date(UTCI_Roma$Date)

# Define the plotting function
plot_utci <- function(data, region) {
  ggplot(data, aes(x = Date)) +
    geom_line(aes(y = Mean,    color = "Mean"),    size = 0.5) +
    geom_line(aes(y = Min,     color = "Minimum"), size = 0.5) +
    geom_line(aes(y = Max,     color = "Maximum"), size = 0.5) +
    labs(title = region,
         x     = "Year",
         y     = "UTCI",
         color = "UTCI") +
    scale_color_manual(
      values = c("Minimum" = "blue", "Mean" = "black", "Maximum" = "red"), 
      breaks = c("Minimum", "Mean", "Maximum")
    ) +
    scale_y_continuous(
      breaks = c(-40, -20, 0, 20, 40, 60),
      limits = c(-40, 60),
      expand = expansion(mult = 0, add = 0)
    ) +
    scale_x_date(
      limits      = c(as.Date("2015-01-01"), as.Date("2020-01-01")),
      date_breaks = "1 year",
      date_labels = "%Y",
      expand      = c(0, 0.02)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border  = element_rect(color = "black", fill = NA, size = 1),
      panel.grid    = element_blank(),
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "bottom",
      plot.margin     = margin(t = 1, r = 20, b = 1, l = 5, unit = "pt")
    )
}

plot_Athens <- plot_utci(UTCI_Attiki, "Athens")
plot_lisbon <- plot_utci(UTCI_Lisbon, "Lisbon")
plot_rome   <- plot_utci(UTCI_Roma,   "Rome")

UTCI_plot <- ggarrange(
  plot_utci(UTCI_Attiki, "Athens"),
  plot_utci(UTCI_Lisbon, "Lisbon"),
  plot_utci(UTCI_Roma,   "Rome"),
  ncol          = 1,
  nrow          = 3,
  common.legend = TRUE,
  legend        = "bottom"
)


# Save the plot
ggsave(filename = "Figures/UTCI/Historical_UTCI.pdf",
       plot = UTCI_plot, width = 9, height = 6, units = "in", dpi = 300)

