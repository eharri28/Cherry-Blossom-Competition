library(tidyverse)
library(lubridate)
library(randomForest)
library(jsonlite)
library(httr)

system("git --version")

#Load and Clean Blossom Data
file_paths <- list(
  "Japan" = "/Users/emmaharris/Desktop/Blossom data/data/japan.csv",
  "Kyoto" = "/Users/emmaharris/Desktop/Blossom data/data/kyoto.csv",
  "Liestal" = "/Users/emmaharris/Desktop/Blossom data/data/liestal.csv",
  "MeteoSwiss" = "/Users/emmaharris/Desktop/Blossom data/data/meteoswiss.csv",
  "NYC" = "/Users/emmaharris/Desktop/Blossom data/data/nyc.csv",
  "SouthKorea" = "/Users/emmaharris/Desktop/Blossom data/data/south_korea.csv",
  "Vancouver" = "/Users/emmaharris/Desktop/Blossom data/data/vancouver.csv",
  "WashingtonDC" = "/Users/emmaharris/Desktop/Blossom data/data/washingtondc.csv"
)

datasets <- lapply(file_paths, read_csv)
combined_df <- bind_rows(datasets, .id = "location")

colnames(combined_df) <- tolower(gsub(" ", "_", colnames(combined_df)))

combined_df <- combined_df %>%
  mutate(
    bloom_date = as.Date(bloom_date, format = "%Y-%m-%d"),
    bloom_doy = as.numeric(bloom_doy)
  ) %>%
  drop_na(bloom_date, bloom_doy, year) %>%
  filter(year >= 2009)

print(paste("Filtered combined_df to years 2009-2024: Rows =", nrow(combined_df)))

#Load and Clean Observed Solar Data
observed_url <- "https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json"
observed_data <- fromJSON(observed_url)

observed_df <- as.data.frame(observed_data) %>%
  select(`time-tag`, ssn, f10.7) %>%
  mutate(
    year = as.numeric(substr(`time-tag`, 1, 4)),
    sunspot_number = as.numeric(ssn),
    radio_flux = as.numeric(f10.7)
  ) %>%
  select(year, sunspot_number, radio_flux) %>%
  group_by(year) %>%
  summarize(
    sunspot_number = last(sunspot_number),
    radio_flux = last(radio_flux)
  )

print(paste("Rows in observed_df:", nrow(observed_df)))

predicted_url <- "https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json"
predicted_data <- fromJSON(predicted_url)

predicted_df <- as.data.frame(predicted_data) %>%
  select(`time-tag`, predicted_ssn, predicted_f10.7) %>%
  mutate(
    year = as.numeric(substr(`time-tag`, 1, 4)),
    predicted_sunspot_number = as.numeric(predicted_ssn),
    predicted_radio_flux = as.numeric(predicted_f10.7)
  ) %>%
  select(year, predicted_sunspot_number, predicted_radio_flux)

print(paste("Rows in predicted_df:", nrow(predicted_df)))

# Fetch GOES X-ray Flux Data (last 7 days)
xray_url <- "https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json"
xray_data <- fromJSON(xray_url)

xray_df <- xray_data %>%
  select(time_tag, flux) %>%
  mutate(
    timestamp = as.POSIXct(time_tag, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
    flux = as.numeric(flux)
  ) %>%
  drop_na()

xray_daily <- xray_df %>%
  mutate(date = as.Date(timestamp)) %>%
  group_by(date) %>%
  summarize(daily_xray_flux = mean(flux))

print(paste("Rows in xray_daily:", nrow(xray_daily)))

#Load global temperature anomaly data from NASA GISTEMP
temp_url <- "https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv"
temp_data <- read.csv(temp_url, skip = 1, header = TRUE)

# Check column names
print(colnames(temp_data))

#Clean and filter temperature data
temp_data <- temp_data %>%
  rename(year = Year) %>%
  select(year, J.D) %>%  # Use correct column name for global anomalies
  filter(year >= 2009) %>%
  rename(temp_anomaly = J.D) %>%
  mutate(temp_anomaly = as.numeric(gsub("[^0-9.-]", "", temp_anomaly)))

print(paste("Rows in temp_data:", nrow(temp_data)))

#Merge All Data
final_df <- combined_df %>%
  left_join(observed_df, by = "year") %>%
  left_join(predicted_df, by = "year") %>%
  left_join(xray_daily, by = c("bloom_date" = "date")) %>%
  left_join(temp_data, by = "year") %>%
  filter(!is.na(bloom_doy))

print(paste("Final dataset rows after merging:", nrow(final_df)))

#Define weight factors for each influence
temp_anomaly_weight <- 0.5  
sunspot_weight <- 0.3  
base_temp_weight <- 0.2

#Adjust bloom DOY using weighted factors
final_df <- final_df %>%
  mutate(
    adjusted_bloom_doy = bloom_doy - (temp_anomaly_weight * temp_anomaly) +
      (sunspot_weight * sunspot_number / 100)
  )

print(head(final_df))

#Predict the bloom DOY for 2025 based on climate and solar activity
#Function to predict bloom DOY for each location
predict_bloom_loess <- function(df, location) {
  location_df <- df %>% filter(location == !!location)  # Filter data for the location
  
  #Ensure enough data points for LOESS
  if (nrow(location_df) < 5) {
    #If too few data points, use median instead
    pred_doy <- median(location_df$adjusted_bloom_doy, na.rm = TRUE)
  } else {
    #Fit LOESS model
    loess_model <- tryCatch(
      loess(adjusted_bloom_doy ~ year, data = location_df, span = 0.5),
      error = function(e) return(NULL)  # Return NULL if LOESS fails
    )
    
    if (!is.null(loess_model)) {
      pred_doy <- predict(loess_model, newdata = data.frame(year = 2025), se = FALSE)
      
      #If LOESS prediction fails, use linear interpolation
      if (is.na(pred_doy)) {
        linear_model <- lm(adjusted_bloom_doy ~ year, data = location_df)
        pred_doy <- predict(linear_model, newdata = data.frame(year = 2025))
      }
    } else {
      #Fall back to median if both methods fail
      pred_doy <- median(location_df$adjusted_bloom_doy, na.rm = TRUE)
    }
  }
  
  #Convert DOY to Date
  pred_date <- as.Date(paste0("2025-01-01")) + round(pred_doy)
  
  return(tibble(predicted_bloom_doy = pred_doy, predicted_date = pred_date))
}

#Apply prediction to each location
predicted_blooms_2025 <- final_df %>%
  group_by(location) %>%  # Group by location instead of rowwise()
  summarize(predicted_bloom_doy = predict_bloom_loess(cur_data(), location)) %>%
  unnest_wider(predicted_bloom_doy)  # Avoid duplicate column names

#Print result
print(predicted_blooms_2025)
