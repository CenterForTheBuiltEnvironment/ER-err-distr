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
             plot.subtitle = element_text(size = 12, color = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
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

alpha_values <- c(1, 0.7, 0.4)

emissions <- "aer"
# emissions <- "moer"


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
figs_path <- "../figs/graphical/"

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




#### CAISO ####
gea_example <- "CAISO"
sce_example <- "Decarb100by2035"
# index <- 

df_annual %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Annual avg."), 
         year = factor(year)) %>% 
  ggplot() +
  geom_col(aes(x = year, y = er, fill = type, alpha = year)) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_fill_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 250, by = 50)) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = NULL,
       y = "Emissions rate",
       fill = NULL,
       title = "Annual averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "annual.svg", path = figs_path, units = "in", height = 8, width = 8, dpi = 300)

df_season %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Season avg."), 
         year = as.factor(year), 
         season = factor(season, levels = c("spring", "summer", "fall", "winter"))) %>% 
  ggplot() +
  geom_col(aes(x = season, y = er, fill = type, alpha = year), position = "dodge") +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_fill_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 250, by = 50)) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = NULL,
       y = "Emissions rate",
       fill = NULL,
       title = "Season averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "seaons.svg", path = figs_path, units = "in", height = 8, width = 8, dpi = 300)

df_tod %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Time-of-day avg."), 
         year = as.factor(year)) %>% 
  ggplot() +
  geom_line(aes(x = hour, y = er, color = type, alpha = year), linewidth = 3) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 250, by = 50)) +
  coord_cartesian(ylim = c(0, 250)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = c(6, 12, 18), 
                     labels = c("6 AM", "12 PM", "6 PM")) +
  labs(x = NULL,
       y = "Emissions rate",
       color = NULL,
       title = "Time-of-day averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "tod.svg", path = figs_path, units = "in", height = 8, width = 8, dpi = 300)

df_season_hour %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Season-hour avg."), 
         year = as.factor(year), 
         season = factor(season, levels = c("spring", "summer", "fall", "winter"))) %>% 
  ggplot() +
  geom_line(aes(x = hour, y = er, color = type, alpha = year), linewidth = 2) +
  facet_wrap(~season) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 350, by = 100)) +
  coord_cartesian(ylim = c(0, 350)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = c(6, 12, 18), 
                     labels = c("6 AM", "12 PM", "6 PM")) +
  labs(x = NULL,
       y = "Emissions rate",
       color = NULL,
       title = "Season-hour averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "season_hour.svg", path = figs_path, units = "in", height = 8, width = 12, dpi = 300)

df_month_hour %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Month-hour avg."), 
         year = as.factor(year), 
         month = factor(month.abb[month],levels = month.abb)) %>% 
  ggplot() +
  geom_line(aes(x = hour, y = er, color = type, alpha = year), linewidth = 1.5) +
  facet_wrap(~month, nrow = 2) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 350, by = 100)) +
  coord_cartesian(ylim = c(0, 350)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = c(6, 18), 
                     labels = c("6 AM", "6 PM")) +
  labs(x = NULL,
       y = "Emissions rate",
       color = NULL,
       title = "Month-hour averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "month_hour.svg", path = figs_path, units = "in", height = 8, width = 18, dpi = 300)

df_s_hourly[[1]] %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Hourly avg."), 
         year = as.factor(year)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = er, color = type, group = year), size = 0.1, alpha = 0.2) +
  geom_smooth(aes(x = datetime, y = er, color = type, group = year), size = 0.8) +
  facet_wrap(~year, nrow = 3, scales = "free_x") +
  scale_color_manual(values = ls_colors) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 350, by = 100)) +
  coord_cartesian(ylim = c(0, 350)) +
  labs(x = NULL,
       y = "Emissions rate",
       fill = NULL,
       title = "Time-of-day averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "hour.svg", path = figs_path, units = "in", height = 8, width = 8, dpi = 300)

  
# all scenarios at hourly resolution


# Summary
sce_example <- c("MidCase", "HighNGPrice", "Decarb95by2050", "Decarb100by2035")
plot_list <- list()
n <- 1
for (i in index){
  
  sce <- sce_example[[n]]
  
  # range
  p <- df_s_hourly[[i]] %>% 
    filter(gea == gea_example, 
           scenario %in% sce_example, 
           year %in% c(2025, 2030, 2050)) %>% 
    mutate(type = as.factor("Hourly avg."), 
           year = as.factor(year)) %>% 
    ggplot() +
    geom_lv(aes(x = year, y = er), alpha = 0.4, k = 4, outlier.size = 0.4) +
    geom_boxplot(aes(x = year, y = er), outlier.alpha = 0, coef = 0, fill = "#00000000") +
    geom_point(data = df_annual %>% 
                 filter(gea == gea_example, 
                        scenario == sce, 
                        year %in% c(2025, 2030, 2050)) %>% 
                 mutate(type = as.factor("Annual avg."), 
                        year = as.factor(year)), 
               aes(x = year, y = er, color = "Annual avg."), 
               size = 3, 
               alpha = 0.5, 
               position = position_nudge(x = -0.1)) +
    geom_point(data = df_season %>% 
                 filter(gea == gea_example, 
                        scenario == sce, 
                        year %in% c(2025, 2030, 2050)) %>% 
                 mutate(type = as.factor("Season avg."), 
                        year = as.factor(year)), 
               aes(x = year, y = er, color = "Season avg."), 
               size = 3, 
               alpha = 0.5, 
               position = position_nudge(x = 0.1)) +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = seq(0, 500, by = 100)) +
    scale_color_manual(values = ls_colors) +
    coord_cartesian(ylim = c(-20, 500)) 
  
  if (n %in% c(1, 3)){
    p <- p + 
      labs(x = NULL,
           y = "Emissions rate (gCO2e/kWh)",
           color = NULL,
           subtitle = str_glue("{sce}")) +
      theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
            axis.text = element_text(size = 12), 
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  } else {
    p <- p + 
      labs(x = NULL,
           y = NULL,
           color = NULL,
           subtitle = str_glue("{sce}")) +
      theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
            axis.text = element_text(size = 12), 
            axis.text.y = element_blank(), 
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  }
  
  plot_list[[n]] <- p
  
  # hour-exceeding
  df_s_hourly[[i]] %>% 
    filter(gea == gea_example, 
           scenario == sce_example, 
           year %in% c(2025, 2030, 2050)) %>% 
    select(datetime, year, hourly = er) %>% 
    left_join(df_annual %>% 
                filter(gea == gea_example, 
                       scenario == sce,
                       year %in% c(2025, 2030, 2050)) %>% 
                select(year, annual = er), 
              by = c("year")) %>% 
    mutate(month = month(datetime), 
           season = ifelse(month >= 3 & month <= 5, "spring", 
                           ifelse(month >= 6 & month <= 8, "summer", 
                                  ifelse(month >= 9 & month <= 11, "fall", 
                                         "winter")))) %>% 
    left_join(df_season %>% 
                filter(gea == gea_example, 
                       scenario == sce,
                       year %in% c(2025, 2030, 2050)) %>% 
                select(year, season, seasonal = er), 
              by = c("season", "year")) %>% 
    select(-c(month, season)) %>% 
    pivot_longer(c(annual, seasonal), names_to = "alt", values_to = "er") %>% 
    mutate(exceed = abs(hourly - er)) %>% 
    group_by(year, alt) %>% 
    summarise(exceed_dev = sd(exceed), 
              exceed_med = median(exceed), 
              exceed_frac = exceed_dev / exceed_med) %>%
    ungroup() %>% 
    print()
  
  n <- n + 1
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          legend = "bottom") +
  plot_annotation(title = "Illustration of differences in emissions factors")


# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in seq(0, 100, by = 50)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/season_hour.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate annual error distribution
    annual_err <- abs((er_annual - er_hourly) / er_hourly * 100)
    
    # Calculate tod error distribution
    tod_err <- abs((er_tod - er_hourly) / er_hourly * 100)
    
    # Calculate month-hour error distribution
    month_hour_err <- abs((er_month_hour - er_hourly) / er_hourly * 100)
    
    # Calculate season error distribution
    season_err <- abs((er_season - er_hourly) / er_hourly * 100)
    
    # Calculate season-hour error distribution
    season_hour_err <- abs((er_season_hour - er_hourly) / er_hourly * 100)
    
    # Combine plots for all
    df_error <- bind_rows(
      annual_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Annual avg."), 
      season_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Season avg."), 
      tod_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Time-of-day avg."),
      season_hour_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Season-hour avg."),
      month_hour_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
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
    
    if (z_index %% 4 == 1){
      
      if (z_index / 4 <= 1){
        
        p <- p + 
          labs(x = NULL, 
               y = str_glue("{perc}% PV generation"), 
               color = NULL, 
               fill = NULL, 
               subtitle = str_glue("{z}")) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      } else {
        
        p <- p + 
          labs(x = NULL, 
               y = str_glue("{perc}% PV generation"), 
               color = NULL, 
               fill = NULL) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      }
    
      
    } else {
      
      if (z_index / 4 <= 1){
        
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
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      } else {
      
        p <- p +
          labs(x = NULL, 
               y = NULL, 
               color = NULL, 
               fill = NULL) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                axis.text.y = element_blank(), 
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
      }
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 4, nrow = 3,
          align = "hv",
          common.legend = TRUE,
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions", 
                  subtitle = str_glue("{gea_example}"))

# ggsave(filename = str_glue("operational/{gea_example}_{perc}_summary.png"), path = figs_path, units = "in", height = 10, width = 16, dpi = 300)


# avoided

# subfigs_path <- paste0(figs_path, str_glue("avoided/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate annual error distribution
    annual_err <- abs((er_annual - er_hourly) / er_hourly * 100)
    
    # Calculate tod error distribution
    tod_err <- abs((er_tod - er_hourly) / er_hourly * 100)
    
    # Calculate month-hour error distribution
    month_hour_err <- abs((er_month_hour - er_hourly) / er_hourly * 100)
    
    # Calculate season error distribution
    season_err <- abs((er_season - er_hourly) / er_hourly * 100)
    
    # Calculate season-hour error distribution
    season_hour_err <- abs((er_season_hour - er_hourly) / er_hourly * 100)
    
    # Combine plots for all
    df_error <- bind_rows(
      annual_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Annual avg."), 
      season_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Season avg."), 
      tod_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Time-of-day avg."),
      season_hour_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Season-hour avg."),
      month_hour_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Month-hour avg."))
    
    breaks <- if (perc == 25) seq(0, 4000, by = 1000) else seq(0, 300, by = 100)
    ylim <- if (perc == 25) c(0, 4200) else c(0, 320)
      
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      filter(year %in% c(2025, 2030, 2050)) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = breaks, 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 4 == 1){
      
      if (z_index / 4 <= 1){
        
        p <- p + 
          labs(x = NULL, 
               y = str_glue("{perc}% PV generation"), 
               color = NULL, 
               fill = NULL, 
               subtitle = str_glue("{z}")) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      } else {
        
        p <- p + 
          labs(x = NULL, 
               y = str_glue("{perc}% PV generation"), 
               color = NULL, 
               fill = NULL) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      }
      
      
    } else {
      
      if (z_index / 4 <= 1){
        
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
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      } else {
        
        p <- p +
          labs(x = NULL, 
               y = NULL, 
               color = NULL, 
               fill = NULL) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                axis.text.y = element_blank(), 
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
      }
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}
  
ggarrange(plotlist = plot_list, 
          ncol = 4, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          legend = "bottom") +
  plot_annotation(title = "Error distribution of avoided carbon emissions", 
                  subtitle = str_glue("{gea_example}"))

# ggsave(filename = str_glue("avoided/{gea_example}_{perc}_summary.png"), path = figs_path, units = "in", height = 10, width = 16, dpi = 300)





#### ERCOT ####
gea_example <- "ERCOT"
sce_example <- "Decarb95by2050"
# index <- 

df_annual %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Annual avg."), 
         year = factor(year)) %>% 
  ggplot() +
  geom_col(aes(x = year, y = er, fill = type, alpha = year)) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_fill_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 250, by = 50)) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = NULL,
       y = "Emissions rate",
       fill = NULL,
       title = "Annual averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "annual.svg", path = figs_path, units = "in", height = 8, width = 8, dpi = 300)

df_season %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Season avg."), 
         year = as.factor(year), 
         season = factor(season, levels = c("spring", "summer", "fall", "winter"))) %>% 
  ggplot() +
  geom_col(aes(x = season, y = er, fill = type, alpha = year), position = "dodge") +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_fill_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 250, by = 50)) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = NULL,
       y = "Emissions rate",
       fill = NULL,
       title = "Season averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "seaons.svg", path = figs_path, units = "in", height = 8, width = 8, dpi = 300)

df_tod %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Time-of-day avg."), 
         year = as.factor(year)) %>% 
  ggplot() +
  geom_line(aes(x = hour, y = er, color = type, alpha = year), linewidth = 3) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 250, by = 50)) +
  coord_cartesian(ylim = c(0, 250)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = c(6, 12, 18), 
                     labels = c("6 AM", "12 PM", "6 PM")) +
  labs(x = NULL,
       y = "Emissions rate",
       color = NULL,
       title = "Time-of-day averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "tod.svg", path = figs_path, units = "in", height = 8, width = 8, dpi = 300)

df_season_hour %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Season-hour avg."), 
         year = as.factor(year), 
         season = factor(season, levels = c("spring", "summer", "fall", "winter"))) %>% 
  ggplot() +
  geom_line(aes(x = hour, y = er, color = type, alpha = year), linewidth = 2) +
  facet_wrap(~season) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 350, by = 100)) +
  coord_cartesian(ylim = c(0, 350)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = c(6, 12, 18), 
                     labels = c("6 AM", "12 PM", "6 PM")) +
  labs(x = NULL,
       y = "Emissions rate",
       color = NULL,
       title = "Season-hour averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "season_hour.svg", path = figs_path, units = "in", height = 8, width = 12, dpi = 300)

df_month_hour %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Month-hour avg."), 
         year = as.factor(year), 
         month = factor(month.abb[month],levels = month.abb)) %>% 
  ggplot() +
  geom_line(aes(x = hour, y = er, color = type, alpha = year), linewidth = 1.5) +
  facet_wrap(~month, nrow = 2) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 350, by = 100)) +
  coord_cartesian(ylim = c(0, 350)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = c(6, 18), 
                     labels = c("6 AM", "6 PM")) +
  labs(x = NULL,
       y = "Emissions rate",
       color = NULL,
       title = "Month-hour averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "month_hour.svg", path = figs_path, units = "in", height = 8, width = 18, dpi = 300)

df_s_hourly[[2]] %>% 
  filter(gea == gea_example, 
         scenario == sce_example, 
         year %in% c(2025, 2030, 2050)) %>% 
  mutate(type = as.factor("Hourly avg."), 
         year = as.factor(year)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = er, color = type, group = year), size = 0.1, alpha = 0.2) +
  geom_smooth(aes(x = datetime, y = er, color = type, group = year), size = 0.8) +
  facet_wrap(~year, nrow = 3, scales = "free_x") +
  scale_color_manual(values = ls_colors) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 350, by = 100)) +
  coord_cartesian(ylim = c(0, 350)) +
  labs(x = NULL,
       y = "Emissions rate",
       fill = NULL,
       title = "Time-of-day averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "hour.svg", path = figs_path, units = "in", height = 8, width = 8, dpi = 300)

# summary
sce_example <- c("MidCase", "HighNGPrice", "Decarb95by2050", "Decarb100by2035")
n <- 1
for (i in index){
  
  sce <- sce_example[[n]]
  
  # range
  df_s_hourly[[i]] %>% 
    filter(gea == gea_example, 
           scenario %in% sce_example, 
           year %in% c(2025, 2030, 2050)) %>% 
    mutate(type = as.factor("Hourly avg."), 
           year = as.factor(year)) %>% 
    ggplot() +
    geom_lv(aes(x = year, y = er), alpha = 0.4, k = 4, outlier.size = 0.4) +
    geom_boxplot(aes(x = year, y = er), outlier.alpha = 0, coef = 0, fill = "#00000000") +
    geom_point(data = df_annual %>% 
                 filter(gea == gea_example, 
                        scenario == sce, 
                        year %in% c(2025, 2030, 2050)) %>% 
                 mutate(type = as.factor("Annual avg."), 
                        year = as.factor(year)), 
               aes(x = year, y = er, color = "Annual avg."), 
               size = 4, 
               alpha = 0.5, 
               position = position_nudge(x = -0.1)) +
    geom_point(data = df_season %>% 
                 filter(gea == gea_example, 
                        scenario == sce, 
                        year %in% c(2025, 2030, 2050)) %>% 
                 mutate(type = as.factor("Season avg."), 
                        year = as.factor(year)), 
               aes(x = year, y = er, color = "Season avg."), 
               size = 4, 
               alpha = 0.5, 
               position = position_nudge(x = 0.1)) +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = seq(0, 500, by = 100), 
                       labels = number_format(suffix = " gCO2e/kWh")) +
    scale_color_manual(values = ls_colors) +
    coord_cartesian(ylim = c(-20, 500)) +
    labs(x = NULL,
         y = "Emissions rate",
         color = NULL,
         title = "Time-of-day averaged emissions rate", 
         subtitle = str_glue("{sce}")) +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          axis.text = element_text(size = 12), 
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  # hour-exceeding
  df_s_hourly[[i]] %>% 
    filter(gea == gea_example, 
           scenario == sce_example, 
           year %in% c(2025, 2030, 2050)) %>% 
    select(datetime, year, hourly = er) %>% 
    left_join(df_annual %>% 
                filter(gea == gea_example, 
                       scenario == sce,
                       year %in% c(2025, 2030, 2050)) %>% 
                select(year, annual = er), 
              by = c("year")) %>% 
    mutate(month = month(datetime), 
           season = ifelse(month >= 3 & month <= 5, "spring", 
                           ifelse(month >= 6 & month <= 8, "summer", 
                                  ifelse(month >= 9 & month <= 11, "fall", 
                                         "winter")))) %>% 
    left_join(df_season %>% 
                filter(gea == gea_example, 
                       scenario == sce,
                       year %in% c(2025, 2030, 2050)) %>% 
                select(year, season, seasonal = er), 
              by = c("season", "year")) %>% 
    select(-c(month, season)) %>% 
    pivot_longer(c(annual, seasonal), names_to = "alt", values_to = "er") %>% 
    mutate(exceed = abs(hourly - er)) %>% 
    group_by(year, alt) %>% 
    summarise(exceed_dev = sd(exceed), 
              exceed_med = median(exceed), 
              exceed_frac = exceed_dev / exceed_med) %>% 
    ungroup() %>% 
    print()
  
  n <- n + 1
  
}


# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in seq(0, 100, by = 50)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/operational/{gea_example}/{perc}/season_hour.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate annual error distribution
    annual_err <- abs((er_annual - er_hourly) / er_hourly * 100)
    
    # Calculate tod error distribution
    tod_err <- abs((er_tod - er_hourly) / er_hourly * 100)
    
    # Calculate month-hour error distribution
    month_hour_err <- abs((er_month_hour - er_hourly) / er_hourly * 100)
    
    # Calculate season error distribution
    season_err <- abs((er_season - er_hourly) / er_hourly * 100)
    
    # Calculate season-hour error distribution
    season_hour_err <- abs((er_season_hour - er_hourly) / er_hourly * 100)
    
    # Combine plots for all
    df_error <- bind_rows(
      annual_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Annual avg."), 
      season_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Season avg."), 
      tod_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Time-of-day avg."),
      season_hour_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Season-hour avg."),
      month_hour_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
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
    
    if (z_index %% 4 == 1){
      
      if (z_index / 4 <= 1){
        
        p <- p + 
          labs(x = NULL, 
               y = str_glue("{perc}% PV generation"), 
               color = NULL, 
               fill = NULL, 
               subtitle = str_glue("{z}")) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      } else {
        
        p <- p + 
          labs(x = NULL, 
               y = str_glue("{perc}% PV generation"), 
               color = NULL, 
               fill = NULL) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      }
      
      
    } else {
      
      if (z_index / 4 <= 1){
        
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
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      } else {
        
        p <- p +
          labs(x = NULL, 
               y = NULL, 
               color = NULL, 
               fill = NULL) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                axis.text.y = element_blank(), 
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
      }
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 4, nrow = 3,
          align = "hv",
          common.legend = TRUE,
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions", 
                  subtitle = str_glue("{gea_example}"))

# ggsave(filename = str_glue("operational/{gea_example}_{perc}_summary.png"), path = figs_path, units = "in", height = 10, width = 16, dpi = 300)


# avoided

# subfigs_path <- paste0(figs_path, str_glue("avoided/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate annual error distribution
    annual_err <- abs((er_annual - er_hourly) / er_hourly * 100)
    
    # Calculate tod error distribution
    tod_err <- abs((er_tod - er_hourly) / er_hourly * 100)
    
    # Calculate month-hour error distribution
    month_hour_err <- abs((er_month_hour - er_hourly) / er_hourly * 100)
    
    # Calculate season error distribution
    season_err <- abs((er_season - er_hourly) / er_hourly * 100)
    
    # Calculate season-hour error distribution
    season_hour_err <- abs((er_season_hour - er_hourly) / er_hourly * 100)
    
    # Combine plots for all
    df_error <- bind_rows(
      annual_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Annual avg."), 
      season_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Season avg."), 
      tod_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Time-of-day avg."),
      season_hour_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Season-hour avg."),
      month_hour_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "_") %>% 
        mutate(type = "Month-hour avg."))
    
    breaks <- if (perc == 25) seq(0, 2000, by = 500) else seq(0, 450, by = 100)
    ylim <- if (perc == 25) c(0, 2200) else c(0, 470)
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      filter(year %in% c(2025, 2030, 2050)) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = breaks, 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 4 == 1){
      
      if (z_index / 4 <= 1){
        
        p <- p + 
          labs(x = NULL, 
               y = str_glue("{perc}% PV generation"), 
               color = NULL, 
               fill = NULL, 
               subtitle = str_glue("{z}")) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      } else {
        
        p <- p + 
          labs(x = NULL, 
               y = str_glue("{perc}% PV generation"), 
               color = NULL, 
               fill = NULL) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      }
      
      
    } else {
      
      if (z_index / 4 <= 1){
        
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
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
        
      } else {
        
        p <- p +
          labs(x = NULL, 
               y = NULL, 
               color = NULL, 
               fill = NULL) +
          theme(panel.grid.major.y = element_line(color = "grey80"),
                legend.direction = "horizontal",
                legend.position = "bottom",
                axis.text.y = element_blank(), 
                plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"))
      }
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 4, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          legend = "bottom") +
  plot_annotation(title = "Error distribution of avoided carbon emissions", 
                  subtitle = str_glue("{gea_example}"))
