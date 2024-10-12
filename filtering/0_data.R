# Genome 2 M&V analysis
# written by Aoyu



#### SETUP ####

# use pacman
require(pacman)

# load libraries
pacman::p_load(tidyverse, lubridate, scales, slider, cowplot, patchwork, RColorBrewer, # general
               broom, ggpmisc, ggpubr, # linear models
               RMV2.0, # TOWT model
               sprtt, effsize) # sequential testing

# turn off scientific notation
options(scipen = 999)

# set default theme for ggplot
theme_set(theme_minimal())

# define base ggplot theme
theme_update(plot.title = element_text(size = 14, colour = "grey20", face = "bold", hjust = 0.5),
             plot.subtitle = element_text(size = 10, colour = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
             plot.caption = element_text(size = 8, colour = "grey20", face = "italic", hjust = 0.5),
             plot.background = element_rect(fill = "white", colour = NA),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             strip.text = element_text(size = 10, color = "grey20", face = "bold"),
             strip.background = element_blank())

# T test
perform_t_test <- function(df_pre, df_post, building) {
  t.test(df_pre[[building]], df_post[[building]], alternative = "two.sided")
}

perform_mean_diff <- function(df_pre, df_post, building) {
  abs(mean(df_pre[[building]], na.rm = T) - mean(df_post[[building]], na.rm = T)) / mean(df_pre[[building]], na.rm = T) * 100
}

perform_mean_est <- function(df_pre, df_post, building) {
  0.5 * (sum(df_pre[[building]], na.rm = T) + sum(df_post[[building]], na.rm = T))
}

process_csv <- function(file) {

  part_a <- read_csv(file) %>% 
    slice(1:5)
  
  part_b <- read_csv(file, skip = 5)
  
  gea <- part_a %>% pull(gea) %>% .[1]
  scenario <- part_a %>% pull(Scenario) %>% .[1]
  year <- as.numeric(part_a %>% pull(t) %>% .[1])
  
  # Add the extracted value as a new column to Part B (or manipulate Part B as needed)
  part_b <- part_b %>% 
    mutate(gea = gea, 
           year = year, 
           scenario = scenario, 
           datetime = as.POSIXct(timestamp_local, format = "%m/%d/%y %H:%M", tz = "UTC")) %>% 
    select(scenario, 
           gea, 
           year, 
           datetime, 
           aer = aer_load_co2e)
  
  # Combine Part B with the extracted value added back in
  return(part_b)
}


#### READ DATA ####
data_path <- "../../../Genome/buds-lab-building-data-genome-project-2/data/meters/cleaned/"
meta_path <- "../../../Genome/buds-lab-building-data-genome-project-2/data/metadata/"
weather_path <- "../../../Genome/buds-lab-building-data-genome-project-2/data/weather/"
output_path <- "../readfiles/"

df_elec <- read_csv(paste0(data_path, "electricity_cleaned.csv"))
df_meta <- read_csv(paste0(meta_path, "metadata.csv"))
df_weather <- read_csv(paste0(weather_path, "weather.csv"))




#### FILTER ####
# Focus
remove_type <- c()
remove_site <- c("Eagle", "Bobcat", "Swan", "Hog", "Gator", "Robin", "Lamb", "Moose", "Wolf", "Shrew", "Mouse", "Crow")

# NAs
na_counts <- sapply(df_elec[-1], function(x) sum(is.na(x)))
cols_to_keep <- names(na_counts[na_counts <= 1440])
cols_to_keep <- c("timestamp", cols_to_keep)

df_elec <- df_elec[, cols_to_keep]

# 0 means
means <- sapply(df_elec[-1], mean, na.rm = TRUE)
cols_to_keep <- names(means[means > 1])
cols_to_keep <- c("timestamp", cols_to_keep)

df_elec <- df_elec[, cols_to_keep]

# yearly difference
df_elec$timestamp = as.POSIXct(df_elec$timestamp, format="%Y-%m-%d %H:%M:%OS")

df_pre <- df_elec %>%
  filter(year(timestamp) == 2016)

df_post <- df_elec %>%
  filter(year(timestamp) == 2017)

buildings <- names(df_elec)[names(df_elec) != "timestamp"]

p_values <- sapply(buildings, function(building) {
  t_test_result <- perform_t_test(df_pre, df_post, building)
  t_test_result$p.value
})

mean_diff <- sapply(buildings, function(building) {
  perform_mean_diff(df_pre, df_post, building)
})

mean_est <- sapply(buildings, function(building) {
  perform_mean_est(df_pre, df_post, building)
})

# electricity intensity
eui_max <- 500

results <- data.frame(Building = buildings, P_Value = p_values, Dev = mean_diff, mean = mean_est) %>% 
  filter(Dev < 25) %>% 
  left_join(df_meta, by = c("Building" = "building_id")) %>% 
  mutate(elec_eui = mean / sqm) %>% 
  filter(elec_eui <= eui_max)

cols_to_keep <- c("timestamp", results %>% .$Building)


df_elec <- df_elec[, cols_to_keep]


# process corresponding metadata
df_meta <- df_meta %>% 
  filter(building_id %in% results$Building) %>% 
  select(building_id, sqm, sqft, eui) %>% 
  separate(building_id, into = c("site", "type", "name"), sep = "_") %>% 
  filter(!type %in% remove_type) %>% 
  filter(!site %in% remove_site)

# process corresponding weather data
df_weather <- df_weather %>% 
  select(timestamp, 
         site = site_id, 
         t_out = airTemperature) %>% 
  mutate(timestamp = as.POSIXct(timestamp, format="%Y-%m-%d %H:%M:%OS")) %>% 
  filter(site %in% df_meta$site)





#### Cambium ####
data_path <- "../../Cambium/"
df_annual <- list.files(path = data_path, pattern = "annual", full.names = TRUE) %>% 
  map_dfr(~ read_csv(.x, skip = 5)) %>% 
  select(scenario, 
         gea, 
         year = t, 
         aer = aer_load_co2e)

df_month_hour <- list.files(path = data_path, pattern = "month-hour", full.names = TRUE) %>% 
  map_dfr(~ read_csv(.x, skip = 5)) %>% 
  select(scenario, 
         gea, 
         year = t, 
         month = month, 
         hour = hour, 
         aer = aer_load_co2e)

df_season <- df_month_hour %>% 
  mutate(season = ifelse(month >= 3 & month <= 5, "spring", 
                         ifelse(month >= 6 & month <= 8, "summer", 
                                ifelse(month >= 9 & month <= 11, "fall", 
                                       "winter")))) %>% 
  group_by(scenario, gea, year, season) %>% 
  summarise(aer = mean(aer)) %>% 
  ungroup()
  

df_season_hour <- df_month_hour %>% 
  mutate(season = ifelse(month >= 3 & month <= 5, "spring", 
                          ifelse(month >= 6 & month <= 8, "summer", 
                                 ifelse(month >= 9 & month <= 11, "fall", 
                                        "winter")))) %>% 
  group_by(scenario, gea, year, season, hour) %>% 
  summarise(aer = mean(aer)) %>% 
  ungroup()

df_tod <- list.files(path = data_path, pattern = "tod", full.names = TRUE) %>% 
  map_dfr(~ read_csv(.x, skip = 5)) %>% 
  select(scenario, 
         gea, 
         year = t, 
         hour = hour, 
         aer = aer_load_co2e)

# hourly Cambium dataset
subfolder_list <- list.dirs(path = paste0(data_path, "hourly"), recursive = FALSE)

s <- 1
for (subfolder in subfolder_list) {
  name <- paste0("df_s", s, "_hourly")
  file_list <- list.files(path = subfolder, pattern = "\\.csv$", full.names = TRUE)
  
  subfolder_data <- file_list %>%
    map_dfr(process_csv)
  
  assign(name, subfolder_data)
  
  s <- s + 1
  
}






#### OUTPUT ####
df_energy <- df_elec %>% 
  pivot_longer(c(-timestamp), names_to = "buildings", values_to = "eload") %>% 
  separate(buildings, into = c("site", "type", "name"), sep = "_") %>% 
  filter(!type %in% remove_type) %>% 
  filter(!site %in% remove_site) 

write_rds(df_energy, paste0(output_path, "df_energy.rds"), compress = "gz")
write_rds(df_meta, paste0(output_path, "df_meta.rds"), compress = "gz")
write_rds(df_weather, paste0(output_path, "df_weather.rds"), compress = "gz")

write_rds(df_annual, paste0(output_path, "df_annual.rds"), compress = "gz")
write_rds(df_season, paste0(output_path, "df_season.rds"), compress = "gz")
write_rds(df_month_hour, paste0(output_path, "df_month_hour.rds"), compress = "gz")
write_rds(df_season_hour, paste0(output_path, "df_season_hour.rds"), compress = "gz")
write_rds(df_tod, paste0(output_path, "df_tod.rds"), compress = "gz")
write_rds(df_s1_hourly, paste0(output_path, "df_s1_hourly.rds"), compress = "gz")
write_rds(df_s2_hourly, paste0(output_path, "df_s2_hourly.rds"), compress = "gz")
write_rds(df_s3_hourly, paste0(output_path, "df_s3_hourly.rds"), compress = "gz")
write_rds(df_s4_hourly, paste0(output_path, "df_s4_hourly.rds"), compress = "gz")
write_rds(df_s5_hourly, paste0(output_path, "df_s5_hourly.rds"), compress = "gz")
write_rds(df_s6_hourly, paste0(output_path, "df_s6_hourly.rds"), compress = "gz")
write_rds(df_s7_hourly, paste0(output_path, "df_s7_hourly.rds"), compress = "gz")
write_rds(df_s8_hourly, paste0(output_path, "df_s8_hourly.rds"), compress = "gz")

