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



#### FUNCTION ####
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
           er = lrmer_co2e)
  
  # Combine Part B with the extracted value added back in
  return(part_b)
}


#### READ DATA ####
data_path <- "../../Cambium/"
output_path <- "../readfiles/moer/"
df_annual <- list.files(path = data_path, pattern = "annual", full.names = TRUE) %>% 
  map_dfr(~ read_csv(.x, skip = 5)) %>% 
  select(scenario, 
         gea, 
         year = t, 
         er = lrmer_co2e)

df_month_hour <- list.files(path = data_path, pattern = "month-hour", full.names = TRUE) %>% 
  map_dfr(~ read_csv(.x, skip = 5)) %>% 
  select(scenario, 
         gea, 
         year = t, 
         month = month, 
         hour = hour, 
         er = lrmer_co2e)

df_season <- df_month_hour %>% 
  mutate(season = ifelse(month >= 3 & month <= 5, "spring", 
                         ifelse(month >= 6 & month <= 8, "summer", 
                                ifelse(month >= 9 & month <= 11, "fall", 
                                       "winter")))) %>% 
  group_by(scenario, gea, year, season) %>% 
  summarise(er = mean(er)) %>% 
  ungroup()


df_season_hour <- df_month_hour %>% 
  mutate(season = ifelse(month >= 3 & month <= 5, "spring", 
                         ifelse(month >= 6 & month <= 8, "summer", 
                                ifelse(month >= 9 & month <= 11, "fall", 
                                       "winter")))) %>% 
  group_by(scenario, gea, year, season, hour) %>% 
  summarise(er = mean(er)) %>% 
  ungroup()

df_tod <- list.files(path = data_path, pattern = "tod", full.names = TRUE) %>% 
  map_dfr(~ read_csv(.x, skip = 5)) %>% 
  select(scenario, 
         gea, 
         year = t, 
         hour = hour, 
         er = lrmer_co2e)

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

