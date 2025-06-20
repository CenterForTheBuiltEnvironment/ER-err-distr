# ER error distribution analysis
# written by Aoyu



#### SETUP ####

# use pacman
require(pacman)

# load libraries
pacman::p_load(tidyverse, lubridate, scales, slider, cowplot, patchwork, RColorBrewer, lvplot, # general
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

ls_colors <- c("Annual avg." = "#E68B81", 
               "Season avg." = "#EAAA60", 
               "Time-of-day avg." = "#B7B2D0",
               "Season-hour avg." = "#7DA6C6",
               "Month-hour avg." = "#84C3B7"
               
)

# emissions <- "aer"
emissions <- "moer"


#### FUNCTION ####
run_interpo <- function(df_all, site_day){
  
  na_counts <- df_all %>%
    group_by(name, date = as.Date(timestamp)) %>%
    summarize(na_hours = sum(is.na(eload))) %>% 
    ungroup()
  
  # Filter out days with more than half of the hours having NAs
  valid_days <- na_counts %>%
    group_by(name) %>% 
    filter(na_hours <= 12) %>%
    ungroup() %>% 
    select(name, date)
  
  df_filtered <- df_all %>%
    mutate(date = as.Date(timestamp)) %>% 
    right_join(valid_days, by = c("name", "date")) %>% 
    drop_na(date)
  
  df_filtered <- df_filtered %>%
    mutate(across(c(eload), ~ zoo::na.approx(., na.rm = FALSE)))
  
  return(df_filtered)
}





#### READ ####
# genome dataset
readfile_path <- "../readfiles/"
output_path <- paste0(readfile_path, "results/", str_glue("{emissions}/"))

gis <- read_csv(paste0(readfile_path, "gis.csv"))
site_map <- read_csv(paste0(readfile_path, "site_map.csv"))
solar_map <- read_csv(paste0(readfile_path, "solar_map.csv"))

df_energy <- read_rds(paste0(readfile_path, "df_energy.rds")) %>% 
  left_join(site_map, by = "site") %>% 
  filter(year(timestamp) == 2016, 
         date(timestamp) != "2016-02-29") %>% 
  select(site = location, 
         timestamp, 
         type, 
         name, 
         eload)

df_compr <- read_rds(paste0(readfile_path, "df_energy.rds")) %>% 
  left_join(site_map, by = "site") %>% 
  filter(year(timestamp) == 2017) %>% 
  select(site = location, 
         timestamp, 
         type, 
         name, 
         eload)

df_meta <- read_rds(paste0(readfile_path, "df_meta.rds"))

df_weather <- read_rds(paste0(readfile_path, "/weather/", "df_weather.rds")) %>% 
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
df_annual <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_annual.rds")) 

df_month_hour <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_month_hour.rds")) 

df_season_hour <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_season_hour.rds"))

df_season <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_season.rds"))

df_tod <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_tod.rds")) 

file_list <- paste0(readfile_path, str_glue("/{emissions}/"), "df_s", 1:8, "_hourly.rds")
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

all_region <- unique(gis %>% drop_na() %>% .$gea)



#### SOLARPV ####
df_hourly <- df_energy %>% 
  run_interpo() %>% 
  select(timestamp, name, eload) %>% 
  pivot_wider(names_from = name, values_from = eload) %>% 
  select(-timestamp)


for (g in all_region){
  
  all_files <- list.files(paste0(readfile_path, "/solar/"), pattern = "\\.csv$", full.names = TRUE)
  
  file_to_read <- all_files[grep(g, all_files)]
  
    
  # convert W/m^2 to kW/m^2
  df_solar <- read_csv(file_to_read[1]) %>% 
    mutate(PV = PV / 1000)
  
  solar_rate <- sum(df_solar$PV)
  
  solar_req <- data.frame()
  
  for (perc in seq(0, 100, by = 25)){
    
    solar_req <- df_hourly %>% 
      summarise_all(~ mean(., na.rm = TRUE) * 8760 * perc / 100) %>% 
      pivot_longer(everything(), names_to = "names", values_to = "values")
    
    PV_area <- solar_req$values / solar_rate
    
    solar_gen <- as.matrix(df_solar$PV) %*% t(as.matrix(PV_area))
    
    df_hourly_net <- as.data.frame(as.matrix(df_hourly) - solar_gen)
      
    }
  
}





#### PROCESS ####
df_hourly <- df_energy %>% 
  run_interpo() %>% 
  select(timestamp, name, eload) %>% 
  pivot_wider(names_from = name, values_from = eload) %>% 
  select(-timestamp)

for (g in all_region){
  
  ifelse(!dir.exists(file.path(str_glue(paste0(output_path, "operational/{g}")))), dir.create(file.path(str_glue(paste0(output_path, "operational/{g}")))), FALSE)
  ifelse(!dir.exists(file.path(str_glue(paste0(output_path, "avoided/{g}")))), dir.create(file.path(str_glue(paste0(output_path, "avoided/{g}")))), FALSE)
  
  all_files <- list.files(paste0(readfile_path, "/solar/"), pattern = "\\.csv$", full.names = TRUE)
  file_to_read <- all_files[grep(g, all_files)]
  
  # convert W/m^2 to kW/m^2
  df_solar <- read_csv(file_to_read[1]) %>% 
    mutate(PV = PV / 1000)
  
  solar_rate <- sum(df_solar$PV)
  
  solar_req <- data.frame()
  
  for (perc in seq(0, 100, by = 25)){
    
    ifelse(!dir.exists(file.path(str_glue(paste0(output_path, "operational/{g}/{perc}/")))), dir.create(file.path(str_glue(paste0(output_path, "operational/{g}/{perc}/")))), FALSE)
    ifelse(!dir.exists(file.path(str_glue(paste0(output_path, "avoided/{g}/{perc}/")))), dir.create(file.path(str_glue(paste0(output_path, "avoided/{g}/{perc}/")))), FALSE)
    
    solar_req <- df_hourly %>% 
      summarise_all(~ mean(., na.rm = TRUE) * 8760 * perc / 100) %>% 
      pivot_longer(everything(), names_to = "names", values_to = "values")
    
    PV_area <- solar_req$values / solar_rate
    
    solar_gen <- as.matrix(df_solar$PV) %*% t(as.matrix(PV_area))
    
    # convert electricity usage from kWh to MWh
    df <- (as.matrix(df_hourly) - solar_gen) / 1000
    df_hourly_net <- as.data.frame(pmax(df, 0))
    df_hourly_avd <- as.data.frame(-ifelse(df < 0, df, 0))
    
    # Annual average rate
    er_annual <- df_annual %>% 
      filter(gea == g) %>% 
      pivot_wider(names_from = c(scenario, year), values_from = er,  names_sep = "-") %>% 
      slice(rep(1, 8760)) %>% 
      select(-gea)
    
    # convert carbon equivalent emissions from kg to tons
    oe_result <- as.data.frame(t(as.matrix(replace(df_hourly_net, is.na(df_hourly_net), 0))) %*% 
                                 as.matrix(er_annual) / 1000)
    
    ae_result <- as.data.frame(t(as.matrix(replace(df_hourly_avd, is.na(df_hourly_avd), 0))) %*% 
                                 as.matrix(er_annual) / 1000)
    
    write_rds(oe_result, str_glue(paste0(output_path, "operational/{g}/{perc}/", "annual.rds")), compress = "gz")
    write_rds(ae_result, str_glue(paste0(output_path, "avoided/{g}/{perc}/", "annual.rds")), compress = "gz")
    
    # Month-hour average rate
    month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    year_hours <- do.call(rbind, lapply(1:12, function(m) {
      data.frame(
        month = m,
        day = 1:month_days[m],
        hour = rep(0:23, each = month_days[m])
      )
    }))
    
    er_month_hour <- year_hours %>%
      left_join(df_month_hour %>%  
                  filter(gea == g) %>% 
                  pivot_wider(names_from = c(scenario, year), values_from = er,  names_sep = "-"), 
                by = c("month", "hour")) %>% 
      arrange(month, day, hour) %>% 
      select(-c(month, day, hour, gea))
    
    oe_result <- as.data.frame(t(as.matrix(replace(df_hourly_net, is.na(df_hourly_net), 0))) %*% 
                                 as.matrix(er_month_hour) / 1000)
    ae_result <- as.data.frame(t(as.matrix(replace(df_hourly_avd, is.na(df_hourly_avd), 0))) %*% 
                                 as.matrix(er_month_hour) / 1000)
    
    write_rds(oe_result, str_glue(paste0(output_path, "operational/{g}/{perc}/", "month_hour.rds")), compress = "gz")
    write_rds(ae_result, str_glue(paste0(output_path, "avoided/{g}/{perc}/", "month_hour.rds")), compress = "gz")
    
    # Season-hour average rate
    year_season_hours <- year_hours %>% 
      mutate(season = ifelse(month >= 3 & month <= 5, "spring", 
                             ifelse(month >= 6 & month <= 8, "summer", 
                                    ifelse(month >= 9 & month <= 11, "fall", 
                                           "winter"))))
    
    er_season_hour <- year_season_hours %>%
      left_join(df_season_hour %>%  
                  filter(gea == g) %>% 
                  pivot_wider(names_from = c(scenario, year), values_from = er,  names_sep = "-"), 
                by = c("season", "hour")) %>% 
      arrange(month, day, hour) %>% 
      select(-c(month, day, hour, gea, season))
    
    oe_result <- as.data.frame(t(as.matrix(replace(df_hourly_net, is.na(df_hourly_net), 0))) %*% 
                                 as.matrix(er_season_hour) / 1000)
    ae_result <- as.data.frame(t(as.matrix(replace(df_hourly_avd, is.na(df_hourly_avd), 0))) %*% 
                                 as.matrix(er_season_hour) / 1000)
    
    write_rds(oe_result, str_glue(paste0(output_path, "operational/{g}/{perc}/", "season_hour.rds")), compress = "gz")
    write_rds(ae_result, str_glue(paste0(output_path, "avoided/{g}/{perc}/", "season_hour.rds")), compress = "gz")
    
    # Season average rate
    er_season <- year_season_hours %>%
      left_join(df_season %>%  
                  filter(gea == g) %>% 
                  pivot_wider(names_from = c(scenario, year), values_from = er,  names_sep = "-"), 
                by = c("season")) %>% 
      arrange(month, day, hour) %>% 
      select(-c(month, day, hour, gea, season))
    
    oe_result <- as.data.frame(t(as.matrix(replace(df_hourly_net, is.na(df_hourly_net), 0))) %*% 
                                 as.matrix(er_season) / 1000)
    ae_result <- as.data.frame(t(as.matrix(replace(df_hourly_avd, is.na(df_hourly_avd), 0))) %*% 
                                 as.matrix(er_season) / 1000)
    
    write_rds(oe_result, str_glue(paste0(output_path, "operational/{g}/{perc}/", "season.rds")), compress = "gz")
    write_rds(ae_result, str_glue(paste0(output_path, "avoided/{g}/{perc}/", "season.rds")), compress = "gz")
    
    # Time-of-day average rate
    er_tod <- year_hours %>%
      left_join(df_tod %>%  
                  filter(gea == g) %>% 
                  pivot_wider(names_from = c(scenario, year), values_from = er,  names_sep = "-"), 
                by = "hour") %>% 
      arrange(month, day, hour) %>% 
      select(-c(month, day, hour, gea))
    
    oe_result <- as.data.frame(t(as.matrix(replace(df_hourly_net, is.na(df_hourly_net), 0))) %*% 
                                 as.matrix(er_tod) / 1000)
    ae_result <- as.data.frame(t(as.matrix(replace(df_hourly_avd, is.na(df_hourly_avd), 0))) %*% 
                                 as.matrix(er_tod) / 1000)
    
    write_rds(oe_result, str_glue(paste0(output_path, "operational/{g}/{perc}/", "tod.rds")), compress = "gz")
    write_rds(ae_result, str_glue(paste0(output_path, "avoided/{g}/{perc}/", "tod.rds")), compress = "gz")
    
    # Hourly
    n_building <- length(df_hourly_net)
    oe_result <- data.frame(matrix(nrow = n_building, ncol = 0))
    ae_result <- data.frame(matrix(nrow = n_building, ncol = 0))
    
    for (z in 1:length(all_scenario$scenario)){
      
      scenario <- unique(df_s_hourly[[z]] %>% .$scenario)
      er_hour <- df_s_hourly[[z]] %>% 
        mutate(toy = format(datetime, "%m-%d %H:%M:%S")) %>% 
        select(-datetime) %>% 
        filter(gea == g) %>% 
        pivot_wider(names_from = c(scenario, year), values_from = er,  names_sep = "-") %>% 
        drop_na(toy) %>% 
        select(-c(gea, toy)) 
      
      oe_s_result <- as.data.frame(t(as.matrix(replace(df_hourly_net, is.na(df_hourly_net), 0))) %*% 
                                     as.matrix(er_hour) / 1000)
      oe_result <- bind_cols(oe_result, oe_s_result)
      
      ae_s_result <- as.data.frame(t(as.matrix(replace(df_hourly_avd, is.na(df_hourly_avd), 0))) %*% 
                                     as.matrix(er_hour) / 1000)
      ae_result <- bind_cols(ae_result, ae_s_result)
    }
    
    write_rds(oe_result, str_glue(paste0(output_path, "operational/{g}/{perc}/", "hourly.rds")), compress = "gz")
    write_rds(ae_result, str_glue(paste0(output_path, "avoided/{g}/{perc}/", "hourly.rds")), compress = "gz")
  
  }
  
}

rm(oe_result, ae_result, oe_s_result, ae_s_result, er_hour, er_tod, er_month_hour, er_annual)





#### ANALYSIS ####
readfile_path <- str_glue("../readfiles/results/{emissions}/")
figs_path <- str_glue("../figs/{emissions}/")

# operational
for (g in all_region){
  
  subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
  
  for (perc in seq(0, 100, by = 25)){
    
    plot_list <- list()
    z_index <- 1
    
    for (z in all_scenario$scenario){
      
      er_hourly <- read_rds(str_glue(paste0(readfile_path, "operational/{g}/{perc}/hourly.rds"))) %>% 
        select(contains(z))
      
      colname <- colnames(er_hourly)
      
      er_annual <- read_rds(str_glue(paste0(readfile_path, "operational/{g}/{perc}/annual.rds"))) %>% 
        select(all_of(colname))
      
      er_tod <- read_rds(str_glue(paste0(readfile_path, "operational/{g}/{perc}/tod.rds"))) %>% 
        select(all_of(colname))
      
      er_month_hour <- read_rds(str_glue(paste0(readfile_path, "operational/{g}/{perc}/month_hour.rds"))) %>% 
        select(all_of(colname))
      
      er_season <- read_rds(str_glue(paste0(readfile_path, "operational/{g}/{perc}/season.rds"))) %>% 
        select(all_of(colname))
      
      er_season_hour <- read_rds(str_glue(paste0(readfile_path, "operational/{g}/{perc}/season_hour.rds"))) %>% 
        select(all_of(colname))
      
      # Calculate annual error distribution
      annual_err <- abs((er_annual - er_hourly) / er_hourly * 100)
      
      mean <- annual_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "annual", 
               gea = g, 
               solar = perc)
      
      # annual_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -5, y = 0.2, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-10, 10)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("annual_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Calculate tod error distribution
      tod_err <- abs((er_tod - er_hourly) / er_hourly * 100)
      
      mean <- tod_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "tod", 
               gea = g, 
               solar = perc)
      
      # tod_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -2.5, y = 0.5, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-5, 5)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("tod_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Calculate month-hour error distribution
      month_hour_err <- abs((er_month_hour - er_hourly) / er_hourly * 100)
      
      mean <- month_hour_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "month-hour", 
               gea = g, 
               solar = perc)
      
      # month_hour_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -1, y = 0.75, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-2.5, 2.5)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("month-hour_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Calculate season error distribution
      season_err <- abs((er_season - er_hourly) / er_hourly * 100)
      
      mean <- season_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "season", 
               gea = g, 
               solar = perc)
      
      # season_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -5, y = 0.2, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-2.5, 2.5)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("season_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Calculate season-hour error distribution
      season_hour_err <- abs((er_season_hour - er_hourly) / er_hourly * 100)
      
      mean <- season_hour_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "season-hour", 
               gea = g, 
               solar = perc)
      
      # season_hour_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -1, y = 0.75, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-2.5, 2.5)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("season-hour_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Combine plots for all
      df_error <- bind_rows(
        annual_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Annual avg."), 
        season_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Season avg."), 
        tod_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Time-of-day avg."),
        season_hour_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Season-hour avg."),
        month_hour_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Month-hour avg."))
      
      p <- df_error %>% 
        mutate(scenario = as.factor(scenario), 
               type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
        filter(year %in% c(2025, 2030, 2050)) %>% 
        ggplot(aes(x = year, y = error, fill = type)) +
        geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
        geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
        scale_y_continuous(expand = c(0, 0), 
                           breaks = seq(0, 50, by = 10), 
                           labels = number_format(suffix = " %")) +
        coord_cartesian(ylim = c(0, 60)) +
        scale_fill_manual(values = ls_colors) +
        scale_color_manual(values = ls_colors)
      
      
      if (z_index == 1 |  z_index == 5){
        
        p <- p + 
          labs(x = NULL, 
               y = "Error distribution", 
               color = NULL, 
               fill = NULL, 
               subtitle = str_glue("{z}")) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      } else {
        
        p <- p +
          labs(x = NULL, 
               y = NULL, 
               color = NULL, 
               fill = NULL, 
               subtitle = str_glue("{z}")) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                axis.text.y = element_blank(), 
                plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      }
      plot_list[[z_index]] <- p
      
      z_index <- z_index + 1
    }
    
    ggarrange(plotlist = plot_list, 
              nrow = 2, ncol = 4, 
              align = "hv",
              common.legend = TRUE,
              legend = "bottom") +
      plot_annotation(title = str_glue("{g} operational carbon emissions\nwith {perc} % solar PV generation"))
    
    ggsave(filename = str_glue("operational/{g}_{perc}_summary.png"), path = figs_path, units = "in", height = 10, width = 16, dpi = 300)
    
  }
  
}

# avoided
for (g in all_region){
  
  subfigs_path <- paste0(figs_path, str_glue("avoided/{g}/"))
  
  for (perc in seq(25, 100, by = 25)){
    
    plot_list <- list()
    z_index <- 1
    
    for (z in all_scenario$scenario){
      
      er_hourly <- read_rds(str_glue(paste0(readfile_path, "avoided/{g}/{perc}/hourly.rds"))) %>% 
        select(contains(z))
      
      colname <- colnames(er_hourly)
      
      er_annual <- read_rds(str_glue(paste0(readfile_path, "avoided/{g}/{perc}/annual.rds"))) %>% 
        select(all_of(colname))
      
      er_tod <- read_rds(str_glue(paste0(readfile_path, "avoided/{g}/{perc}/tod.rds"))) %>% 
        select(all_of(colname))
      
      er_month_hour <- read_rds(str_glue(paste0(readfile_path, "avoided/{g}/{perc}/month_hour.rds"))) %>% 
        select(all_of(colname))
      
      er_season <- read_rds(str_glue(paste0(readfile_path, "avoided/{g}/{perc}/season.rds"))) %>% 
        select(all_of(colname))
      
      er_season_hour <- read_rds(str_glue(paste0(readfile_path, "avoided/{g}/{perc}/season_hour.rds"))) %>% 
        select(all_of(colname))
      
      # Calculate annual error distribution
      annual_err <- abs((er_annual - er_hourly) / er_hourly * 100)
      
      mean <- annual_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "annual", 
               gea = g, 
               solar = perc)
      
      # annual_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -5, y = 0.2, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-10, 10)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("annual_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Calculate tod error distribution
      tod_err <- abs((er_tod - er_hourly) / er_hourly * 100)
      
      mean <- tod_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "tod", 
               gea = g, 
               solar = perc)
      
      # tod_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -2.5, y = 0.5, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-5, 5)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("tod_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Calculate month-hour error distribution
      month_hour_err <- abs((er_month_hour - er_hourly) / er_hourly * 100)
      
      mean <- month_hour_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "month-hour", 
               gea = g, 
               solar = perc)
      
      # month_hour_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -1, y = 0.75, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-2.5, 2.5)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("month-hour_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Calculate season error distribution
      season_err <- abs((er_season - er_hourly) / er_hourly * 100)
      
      mean <- season_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "season", 
               gea = g, 
               solar = perc)
      
      # season_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -5, y = 0.2, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-2.5, 2.5)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("season_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Calculate season-hour error distribution
      season_hour_err <- abs((er_season_hour - er_hourly) / er_hourly * 100)
      
      mean <- season_hour_err %>% 
        pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
        group_by(scenario) %>% 
        summarise(mean = mean(values)) %>% 
        ungroup() %>% 
        mutate(type = "season-hour", 
               gea = g, 
               solar = perc)
      
      # season_hour_err %>% 
      #   pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      #   mutate(scenario = as.factor(scenario)) %>% 
      #   ggplot(aes(x = values, group = scenario)) +
      #   geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      #   geom_density(alpha=.2, fill="#FF6666") +
      #   geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      #   geom_text(data = mean, 
      #             aes(x = -1, y = 0.75, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      #   facet_wrap(~scenario, nrow = 6) +
      #   coord_cartesian(xlim = c(-2.5, 2.5)) +
      #   theme(panel.grid.major.y = element_line(color = "grey80"),
      #         legend.direction = "horizontal",
      #         legend.position = "bottom",
      #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      # 
      # ggsave(filename = str_glue("season-hour_{z}_{perc}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
      
      # Combine plots for all
      df_error <- bind_rows(
        annual_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Annual avg."), 
        season_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Season avg."), 
        tod_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Time-of-day avg."),
        season_hour_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Season-hour avg."),
        month_hour_err %>% 
          pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
          separate(year, into = c("scenario", "year"), sep = "-") %>% 
          mutate(type = "Month-hour avg."))
      
      p <- df_error %>% 
        mutate(scenario = as.factor(scenario), 
               type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
        filter(year %in% c(2025, 2030, 2050)) %>% 
        ggplot(aes(x = year, y = error, fill = type)) +
        geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
        geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
        scale_y_continuous(expand = c(0, 0), 
                           breaks = seq(0, 100, by = 25), 
                           labels = number_format(suffix = " %")) +
        coord_cartesian(ylim = c(0, 120)) +
        scale_fill_manual(values = ls_colors) +
        scale_color_manual(values = ls_colors)
      
      
      if (z_index == 1 |  z_index == 5){
        
        p <- p + 
          labs(x = NULL, 
               y = "Error distribution", 
               color = NULL, 
               fill = NULL, 
               subtitle = str_glue("{z}")) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      } else {
        
        p <- p +
          labs(x = NULL, 
               y = NULL, 
               color = NULL, 
               fill = NULL, 
               subtitle = str_glue("{z}")) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                axis.text.y = element_blank(), 
                plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
      }
      plot_list[[z_index]] <- p
      
      z_index <- z_index + 1
    }
    
    ggarrange(plotlist = plot_list, 
              nrow = 2, ncol = 4, 
              align = "hv",
              common.legend = TRUE,
              legend = "bottom") +
      plot_annotation(title = str_glue("{g} avoided carbon emissions\nwith {perc} % solar PV generation"))
    
    ggsave(filename = str_glue("avoided/{g}_{perc}_summary.png"), path = figs_path, units = "in", height = 10, width = 16, dpi = 300)
    
  }
  
}
