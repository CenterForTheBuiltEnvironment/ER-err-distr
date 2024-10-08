# AER error distribution analysis
# written by Aoyu



#### SETUP ####

# use pacman
require(pacman)

# load libraries
pacman::p_load(tidyverse, lubridate, scales, slider, cowplot, patchwork, RColorBrewer, # general
               broom, ggpmisc, ggpubr, # linear models
               sprtt, effsize) # sequential testing

# turn off scientific notation
options(scipen = 999)

# set default theme for ggplot
theme_set(theme_minimal())

# define base ggplot theme
theme_update(plot.title = element_text(size = 14, color = "grey20", face = "bold", hjust = 0.5),
             plot.subtitle = element_text(size = 10, color = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
             plot.caption = element_text(size = 8, color = "grey20", face = "italic", hjust = 0.5),
             plot.background = element_rect(fill = "white", color = NA),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             strip.text = element_text(size = 10, color = "grey20", face = "bold"),
             strip.background = element_blank())





#### READ ####
# genome dataset
readfile_path <- "../readfiles/"
figs_path <- "../figs/"
output_path <- paste0(readfile_path, "results/")

df_energy <- read_rds(paste0(readfile_path, "df_energy.rds")) %>% 
  filter(year(timestamp) == 2016)

df_compr <- read_rds(paste0(readfile_path, "df_energy.rds")) %>% 
  filter(year(timestamp) == 2017)

df_meta <- read_rds(paste0(readfile_path, "df_meta.rds"))

df_weather <- read_rds(paste0(readfile_path, "df_weather.rds")) %>% 
  filter(year(timestamp) == 2016)

all_sites <- df_energy %>%
  select(site) %>%
  distinct() %>%
  arrange(site)

all_types <- df_energy %>%
  select(type) %>%
  mutate(type = as.factor(type)) %>%
  distinct()

# carbon emissions dataset
gea_map <- read_csv(paste0(readfile_path, "gea_map.csv")) %>% 
  filter(site %in% all_sites$site)

df_annual <- read_rds(paste0(readfile_path, "df_annual.rds")) 

df_month_hour <- read_rds(paste0(readfile_path, "df_month_hour.rds")) 

df_tod <- read_rds(paste0(readfile_path, "df_tod.rds")) 

file_list <- paste0(readfile_path, "df_s", 1:8, "_hourly.rds")
df_s_hourly <- map(file_list, readRDS)

all_names <- df_energy %>%
  select(name, site) %>% 
  distinct()

all_scenario <- df_annual %>% 
  select(scenario) %>% 
  distinct()

all_year <- df_annual %>% 
  select(year) %>% 
  distinct()




#### PROCESS ####
for (s in all_sites$site){
  
  ifelse(!dir.exists(file.path(str_glue(paste0(output_path, "{s}")))), dir.create(file.path(str_glue(paste0(output_path, "{s}")))), FALSE)
  
  df_hourly <- df_energy %>% 
    filter(site == s) %>% 
    select(timestamp, name, eload) %>% 
    mutate(eload = replace_na(eload, 0) / 1000) %>% 
    pivot_wider(id_cols = timestamp,names_from = name, values_from = eload) %>% 
    select(-timestamp) 
  
  g <- gea_map %>% 
    filter(site == s) %>%
    .$gea_name
  
  # Annual average rate
  aer_annual <- df_annual %>% 
    filter(gea == g) %>% 
    pivot_wider(names_from = c(scenario, year), values_from = aer) %>% 
    slice(rep(1, 8784)) %>% 
    select(-gea)
  
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_annual) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{s}/", "annual_result.rds")), compress = "gz")
  
  # Month-hour average rate
  month_days <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  year_hours <- do.call(rbind, lapply(1:12, function(m) {
    data.frame(
      month = m,
      day = 1:month_days[m],
      hour = rep(0:23, each = month_days[m])
    )
  }))
  
  aer_month_hour <- year_hours %>%
    left_join(df_month_hour %>%  
                filter(gea == g) %>% 
                pivot_wider(names_from = c(scenario, year), values_from = aer), 
              by = c("month", "hour")) %>% 
    arrange(month, day, hour) %>% 
    select(-c(month, day, hour, gea))
    
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_month_hour) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{s}/", "month_hour_result.rds")), compress = "gz")
  
  # Time-of-day average rate
  aer_tod <- year_hours %>%
    left_join(df_tod %>%  
                filter(gea == g) %>% 
                pivot_wider(names_from = c(scenario, year), values_from = aer), 
              by = "hour") %>% 
    arrange(month, day, hour) %>% 
    select(-c(month, day, hour, gea))
  
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_month_hour) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{s}/", "tod_result.rds")), compress = "gz")
  
  # Hourly
  for (z in 1:length(all_scenario$scenario)){
    scenario <- unique(df_s_hourly[[z]] %>% .$scenario)
    aer_hour <- df_s_hourly[[z]] %>% 
      mutate(toy = format(datetime, "%m-%d %H:%M:%S")) %>% 
      select(-datetime) %>% 
      filter(gea == g) %>% 
      pivot_wider(names_from = c(scenario, year), values_from = aer) %>% 
      drop_na(toy) %>% 
      select(-c(gea, toy)) 

      
    
    df_hourly <- df_energy %>% 
      filter(site == s) %>% 
      select(timestamp, name, eload) %>% 
      mutate(eload = replace_na(eload, 0) / 1000) %>% 
      filter(date(timestamp) != as.Date("2016-02-29")) %>% 
      pivot_wider(id_cols = timestamp,names_from = name, values_from = eload) %>% 
      select(-timestamp) 
    
    result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_hour) / 1000)
    
    write_rds(result, str_glue(paste0(output_path, "{s}/", "{scenario}_hourly_result.rds")), compress = "gz")
  }
}
