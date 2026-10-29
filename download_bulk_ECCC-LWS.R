# Install and load the required packages
install.packages("devtools") 
devtools::install_github("ropensci/weathercan")
library(weathercan)
library(dplyr)

# Optional: look up the data for stations of interest
infoStations <- stations_search("sachs harbour", interval = "hour")

# Download Sachs Harbour data
sachs_a <- weather_dl(station_ids = 1794, start = "1980-09-01", end = "2020-10-01", time_disp = "UTC")  # 1980-2015
sachs_a2 <- weather_dl(station_ids = 53326, start = "1980-09-01", end = "2020-10-01", time_disp = "UTC") # from 2015 onwards
sachs_a3 <- weather_dl(station_ids = 53328, start = "1980-09-01", end = "2020-10-01", time_disp = "UTC") # from 2015 onwards for weather observations
sachs_a2$weather <- sachs_a3$weather # Add the weather observations
sachs_all = rbind(sachs_a, sachs_a2) # merge both datasets
write.csv(sachs_all, file='./Sachs_A_8020.csv')

#Download Thomsen River data
thomsen <- weather_dl(station_ids = 10598, start = "1980-09-01", end = "2020-10-01", time_disp = "UTC")
write.csv(thomsen, file='./Thomsen_River_8020.csv')

# Download Nangmagvik Lake data
nangmagvik <- weather_dl(station_ids = 27118, start = "1980-09-01", end = "2020-10-01", time_disp = "UTC")
write.csv(nangmagvik, file='./Nangmagvik_Lake_8020.csv')


#df <- weather_dl(station_ids = infoStations$station_id, start = date_start, end = date_end, interval = "hour", time_disp  = "UTC")
#dataStations <- df[(df$temp >= tempmin) & (df$temp <= tempmax) & (df$wind_spd < wind_spd) & (df$rel_hum <= rel_hum),] 

