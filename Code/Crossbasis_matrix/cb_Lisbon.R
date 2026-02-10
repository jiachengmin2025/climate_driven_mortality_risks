## Used Packages
packages <- c("readxl", "ISOweek", "dplyr")
invisible(lapply(packages, library, character.only = TRUE))

## Daily mean UTCI data
UTCI_daily_mean = read_xlsx("Data/UTCI_data/Daily_data/UTCI_daily_mean.xlsx")
Week <- ISOweek(UTCI_daily_mean$Date)
UTCI_daily_mean = UTCI_daily_mean %>% mutate(Week) %>% select(Week, Lisbon) # 2 cols
UTCI_weekly_mean = UTCI_daily_mean %>% group_by(Week) %>%  summarise(Lisbon = mean(Lisbon))


UTCI_daily_mean = UTCI_daily_mean %>%
  mutate(Lisbon = if_else(row_number() <= 4, NA_real_, Lisbon))

UTCI_ext_3 = matrix(nrow = 260, ncol = 22)
colnames(UTCI_ext_3) = paste("lag", 0:21, sep = "")
rownames(UTCI_ext_3) = UTCI_weekly_mean$Week[2:261]


## 21 Lags (Start at the last day of the week) 
lag0 = UTCI_daily_mean %>% group_by(Week) %>% 
  summarise(lag = Lisbon[7]) %>% select(lag) %>% slice(2:261)
lag0 = unlist(lag0[,1])
UTCI_ext_3[,1] = lag0

lag7 = c(NA, lag0[1:259])
UTCI_ext_3[,8] = lag7

lag14 = c(NA, NA, lag0[1:258])
UTCI_ext_3[,15] = lag14

lag21 = c(NA, NA, NA, lag0[1:257])
UTCI_ext_3[,22] = lag21

lag1 = UTCI_daily_mean %>% group_by(Week) %>% 
  summarise(lag = Lisbon[6]) %>% select(lag) %>% slice(1:260)
lag1 = unlist(lag1[,1])
UTCI_ext_3[,2] = lag1

lag8 = c(NA, lag1[1:259])
UTCI_ext_3[,9] = lag8

lag15 = c(NA, NA, lag1[1:258])
UTCI_ext_3[,16] = lag15

lag2 = UTCI_daily_mean %>% group_by(Week) %>% 
  summarise(lag = Lisbon[5]) %>% select(lag) %>% slice(1:260)
lag2 = unlist(lag2[,1])
UTCI_ext_3[,3] = lag2

lag9 = c(NA, lag2[1:259])
UTCI_ext_3[,10] = lag9

lag16 = c(NA, NA, lag2[1:258])
UTCI_ext_3[,17] = lag16

lag3 = UTCI_daily_mean %>% group_by(Week) %>% 
  summarise(lag = Lisbon[4]) %>% select(lag) %>% slice(1:260)
lag3 = unlist(lag3[,1])
UTCI_ext_3[,4] = lag3

lag10 = c(NA, lag3[1:259])
UTCI_ext_3[,11] = lag10

lag17 = c(NA, NA, lag3[1:258])
UTCI_ext_3[,18] = lag17

lag4 = UTCI_daily_mean %>% group_by(Week) %>% 
  summarise(lag = Lisbon[3]) %>% select(lag) %>% slice(1:260)
lag4 = unlist(lag4[,1])
UTCI_ext_3[,5] = lag4

lag11 = c(NA, lag4[1:259])
UTCI_ext_3[,12] = lag11

lag18 = c(NA, NA, lag4[1:258])
UTCI_ext_3[,19] = lag18


lag5 = UTCI_daily_mean %>% group_by(Week) %>% 
  summarise(lag = Lisbon[2]) %>% select(lag) %>% slice(1:260)
lag5 = unlist(lag5[,1])
UTCI_ext_3[,6] = lag5

lag12 = c(NA, lag5[1:259])
UTCI_ext_3[,13] = lag12

lag19 = c(NA, NA, lag5[1:258])
UTCI_ext_3[,20] = lag19

lag6 = UTCI_daily_mean %>% group_by(Week) %>% 
  summarise(lag = Lisbon[1]) %>% select(lag) %>% slice(1:260)
lag6 = unlist(lag6[,1])
UTCI_ext_3[,7] = lag6

lag13 = c(NA, lag5[1:259])
UTCI_ext_3[,14] = lag13


lag20 = c(NA, NA, lag5[1:258])
UTCI_ext_3[,21] = lag20


## Replace with mean
replace_na_with_mean <- function(data) {
  for (i in 1:length(data)) {
    if (is.na(data[i])) {
      left_index <- max(1, i - 1)
      right_index <- min(length(data), i + 1)
      avg <- mean(data[c(left_index, right_index)], na.rm = TRUE)
      data[i] <- avg
    }
  }
  return(data)
}
UTCI_ext.Lisbon = replace_na_with_mean(UTCI_ext_3)


## Save matrix
save(UTCI_ext.Lisbon, file = "Code/Crossbasis_matrix/cb_Lisbon.RData")
