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
output_path <- paste0(readfile_path, "results/")

gea_map <- read_csv(paste0(readfile_path, "gea_map.csv"))

df_energy <- read_rds(paste0(readfile_path, "df_energy.rds")) %>% 
  filter(year(timestamp) == 2016) %>% 
  left_join(gea_map, by = "site") %>% 
  select(site = location, 
         timestamp, 
         type, 
         name, 
         eload)

df_compr <- read_rds(paste0(readfile_path, "df_energy.rds")) %>% 
  filter(year(timestamp) == 2017) %>% 
  left_join(gea_map, by = "site") %>% 
  select(site = location, 
         timestamp, 
         type, 
         name, 
         eload)

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
    filter(location == s) %>%
    .$gea_name
  
  # Annual average rate
  aer_annual <- df_annual %>% 
    filter(gea == g) %>% 
    pivot_wider(names_from = c(scenario, year), values_from = aer) %>% 
    slice(rep(1, 8784)) %>% 
    select(-gea)
  
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_annual) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{s}/", "annual.rds")), compress = "gz")
  
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
  
  write_rds(result, str_glue(paste0(output_path, "{s}/", "month_hour.rds")), compress = "gz")
  
  # Time-of-day average rate
  aer_tod <- year_hours %>%
    left_join(df_tod %>%  
                filter(gea == g) %>% 
                pivot_wider(names_from = c(scenario, year), values_from = aer), 
              by = "hour") %>% 
    arrange(month, day, hour) %>% 
    select(-c(month, day, hour, gea))
  
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_tod) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{s}/", "tod.rds")), compress = "gz")
  
  # Hourly
  n_building <- length(df_hourly)
  result <- data.frame(matrix(nrow = n_building, ncol = 0))
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
    
    s_result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_hour) / 1000)
    
    result <- bind_cols(result, s_result)
  }
  
  write_rds(result, str_glue(paste0(output_path, "{s}/", "hourly.rds")), compress = "gz")
}

rm(result, s_result, aer_hour, aer_tod, aer_month_hour, aer_annual)




#### ANALYSIS ####
readfile_path <- "../readfiles/results/"
figs_path <- "../figs/"

for (s in all_sites$site){
  
  g <- gea_map %>% 
    filter(location == s) %>%
    .$gea_name
  
  for (z in 1:length(all_scenario$scenario)){
    
    aer_hour <- df_s_hourly[[z]] %>% 
      mutate(toy = format(datetime, "%m-%d %H:%M:%S")) %>% 
      filter(gea == g) 
    
  }
}

all_median <- data.frame()

for (s in all_sites$site){
  figs_path <- str_glue("../figs/{s}/")
  
  g <- gea_map %>% 
    filter(location == s) %>%
    .$gea_name
  
  for (z in all_scenario$scenario){
    
    aer_hourly <- read_rds(str_glue(paste0(readfile_path, "{s}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(aer_hourly)
    
    aer_annual <- read_rds(str_glue(paste0(readfile_path, "{s}/annual.rds"))) %>% 
      select(colname)
    
    aer_tod <- read_rds(str_glue(paste0(readfile_path, "{s}/tod.rds"))) %>% 
      select(colname)
    
    aer_month_hour <- read_rds(str_glue(paste0(readfile_path, "{s}/month_hour.rds"))) %>% 
      select(colname)
    
    # Calculate annual error distribution
    annual_err <- (aer_annual - aer_hourly) / aer_hourly * 100
    
    median <- annual_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(median = median(values)) %>% 
      ungroup() %>% 
      mutate(type = "annual", 
             site = s, 
             gea = g)
    
    all_median <- bind_rows(all_median, median)
    
    annual_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      mutate(scenario = as.factor(scenario)) %>% 
      ggplot(aes(x = values, group = scenario)) +
      geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666") +
      geom_vline(data = median, aes(xintercept = median),linewidth = 1) +
      geom_text(data = median, 
                aes(x = -5, y = 0.2, label = paste0("Median: ", round(median, digits = 2), " tons CO2e")), check_overlap = T) +
      facet_wrap(~scenario, nrow = 6) +
      coord_cartesian(xlim = c(-10, 10)) +
      theme(panel.grid.major.y = element_line(color = "grey80"),
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
    
    
    ggsave(filename = str_glue("annual_{z}_dist.png"), path = str_glue(figs_path), units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate tod error distribution
    tod_err <- (aer_tod - aer_hourly) / aer_hourly * 100

    median <- tod_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(median = median(values)) %>% 
      ungroup() %>% 
      mutate(type = "tod", 
             site = s, 
             gea = g)
    
    all_median <- bind_rows(all_median, median)
    
    tod_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      mutate(scenario = as.factor(scenario)) %>% 
      ggplot(aes(x = values, group = scenario)) +
      geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666") +
      geom_vline(data = median, aes(xintercept = median),linewidth = 1) +
      geom_text(data = median, 
                aes(x = -2.5, y = 0.5, label = paste0("Median: ", round(median, digits = 2), " tons CO2e")), check_overlap = T) +
      facet_wrap(~scenario, nrow = 6) +
      coord_cartesian(xlim = c(-5, 5)) +
      theme(panel.grid.major.y = element_line(color = "grey80"),
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
    
    
    ggsave(filename = str_glue("tod_{z}_dist.png"), path = str_glue(figs_path), units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate month-hour error distribution
    month_hour_err <- (aer_hourly - aer_month_hour) / aer_hourly * 100
    
    median <- month_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(median = median(values)) %>% 
      ungroup() %>% 
      mutate(type = "month-hour", 
             site = s, 
             gea = g)
    
    all_median <- bind_rows(all_median, median)
    
    month_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      mutate(scenario = as.factor(scenario)) %>% 
      ggplot(aes(x = values, group = scenario)) +
      geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666") +
      geom_vline(data = median, aes(xintercept = median),linewidth = 1) +
      geom_text(data = median, 
                aes(x = -1, y = 0.75, label = paste0("Median: ", round(median, digits = 2), " tons CO2e")), check_overlap = T) +
      facet_wrap(~scenario, nrow = 6) +
      coord_cartesian(xlim = c(-2.5, 2.5)) +
      theme(panel.grid.major.y = element_line(color = "grey80"),
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
    
    
    ggsave(filename = str_glue("month-hour_{z}_dist.png"), path = str_glue(figs_path), units = "in", height = 8, width = 8, dpi = 300)
    
  }

}

all_median <- bind_rows(all_median)

for (s in all_sites$site){
  
  all_median %>% 
    filter(site == s) %>% 
    mutate(type = as_factor(type), 
           type = recode_factor(type, 
                                "annual" = "Annual avg.", 
                                "month-hour" = "Month-hour avg.", 
                                "tod" = "Time-of-day avg.")) %>% 
    separate(scenario, into = c("scenario", "year"), sep = "_") %>% 
    ggplot() +
    geom_bar(aes(x = year, y = median, fill = type), stat = "identity", position = "dodge") +
    facet_wrap(~ scenario, nrow = 4, scales = "free_y") +
    labs(x = "Projected year", 
         y = "Error distribution median") +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " %")) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = str_glue("{s}"), 
         subtitle = str_glue("{g}"), 
         fill = NULL) +
    theme(panel.grid.major.y = element_line(color = "grey80"),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggsave(filename = str_glue("{s}_overall_summary.png"), path = str_glue("../figs/"), units = "in", height = 10, width = 8, dpi = 300)
}
