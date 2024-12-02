# AER error distribution analysis
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





#### READ ####
# genome dataset
readfile_path <- "../readfiles/"
output_path <- paste0(readfile_path, "results/")

gea_map <- read_csv(paste0(readfile_path, "gea_map.csv"))
site_map <- read_csv(paste0(readfile_path, "site_map.csv"))

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
df_annual <- read_rds(paste0(readfile_path, "/aer/", "df_annual.rds")) 

df_month_hour <- read_rds(paste0(readfile_path, "/aer/", "df_month_hour.rds")) 

df_season_hour <- read_rds(paste0(readfile_path, "/aer/", "df_season_hour.rds"))

df_season <- read_rds(paste0(readfile_path, "/aer/", "df_season.rds"))

df_tod <- read_rds(paste0(readfile_path, "/aer/", "df_tod.rds")) 
  
file_list <- paste0(readfile_path, "/aer/", "df_s", 1:8, "_hourly.rds")
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




#### SITE ####
set3 <- colorRampPalette(brewer.pal('Set3',n=12))
type_colors <- setNames(set3(16), all_types$type)

# site and type number summary
p1 <- df_energy %>%
  group_by(type) %>%
  distinct(name) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(type = fct_reorder(type, n, .desc = F)) %>%
  ggplot(aes(x = 1, y = n, fill = as.factor(type))) +
  geom_col() +
  geom_text(aes(label = ifelse(n > 10, as.character(n), "")), color = "black", position = position_stack(vjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = breaks_pretty(n = 4)) +
  scale_x_discrete(expand = c(0, 0.1)) +
  scale_fill_manual(values = type_colors) +
  labs(x = NULL,
       y = NULL,
       subtitle = "Across all sites",
       fill = NULL) +
  theme(axis.text.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        legend.direction = "horizontal",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- df_energy %>%
  group_by(site, type) %>%
  distinct(name) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n),
         total = sum(n),
         ymax = cumsum(proportion),
         ymin = c(0, head(ymax, n = -1))) %>%
  mutate(label_pos = (ymax + ymin) / 2) %>%
  ungroup() %>%
  group_by(type) %>%
  mutate(order = sum(n)) %>%
  ungroup() %>%
  mutate(type = fct_reorder(type, order, .desc = F)) %>%
  ggplot(aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = type)) +
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  facet_wrap(~ site, nrow = 2) +
  labs(x = NULL,
       y = NULL,
       subtitle = "For each site",
       fill = NULL) +
  scale_fill_manual(values = type_colors) +
  geom_text(aes(x = 3.5, y = label_pos, label = ifelse(n > 1, as.character(n), ""))) +
  geom_text(aes(x = 2, y = 0, label = paste0("Total\n", total)), color = "black") +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p1, p2,
          ncol = 2, nrow = 1,
          widths = c(0.25, 1),
          common.legend = TRUE,
          legend = "bottom") +
  plot_annotation(title = "Case study building type summary")





#### PROCESS ####
for (area in gea_map$gea_name){
  
  g <- gea_map %>% filter(gea_name == area) %>% .$gea
  
  ifelse(!dir.exists(file.path(str_glue(paste0(output_path, "{g}")))), dir.create(file.path(str_glue(paste0(output_path, "{g}")))), FALSE)
  
  # convert electricity usage from kWh to MWh
  df_hourly <- df_energy %>% 
    select(timestamp, name, eload) %>% 
    mutate(eload = replace_na(eload, 0) / 1000) %>% 
    pivot_wider(id_cols = timestamp,names_from = name, values_from = eload) %>% 
    select(-timestamp) 
  
  # Annual average rate
  aer_annual <- df_annual %>% 
    filter(gea == area) %>% 
    pivot_wider(names_from = c(scenario, year), values_from = aer) %>% 
    slice(rep(1, 8760)) %>% 
    select(-gea)
  
  # convert carbon equivalent emissions from kg to tons
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_annual) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{g}/", "annual.rds")), compress = "gz")
  
  # Month-hour average rate
  month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  year_hours <- do.call(rbind, lapply(1:12, function(m) {
    data.frame(
      month = m,
      day = 1:month_days[m],
      hour = rep(0:23, each = month_days[m])
    )
  }))
  
  aer_month_hour <- year_hours %>%
    left_join(df_month_hour %>%  
                filter(gea == area) %>% 
                pivot_wider(names_from = c(scenario, year), values_from = aer), 
              by = c("month", "hour")) %>% 
    arrange(month, day, hour) %>% 
    select(-c(month, day, hour, gea))
    
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_month_hour) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{g}/", "month_hour.rds")), compress = "gz")
  
  # Season-hour average rate
  year_season_hours <- year_hours %>% 
    mutate(season = ifelse(month >= 3 & month <= 5, "spring", 
                           ifelse(month >= 6 & month <= 8, "summer", 
                                  ifelse(month >= 9 & month <= 11, "fall", 
                                         "winter"))))
  
  aer_season_hour <- year_season_hours %>%
    left_join(df_season_hour %>%  
                filter(gea == area) %>% 
                pivot_wider(names_from = c(scenario, year), values_from = aer), 
              by = c("season", "hour")) %>% 
    arrange(month, day, hour) %>% 
    select(-c(month, day, hour, gea, season))
  
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_season_hour) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{g}/", "season_hour.rds")), compress = "gz")
  
  # Season average rate
  aer_season <- year_season_hours %>%
    left_join(df_season %>%  
                filter(gea == area) %>% 
                pivot_wider(names_from = c(scenario, year), values_from = aer), 
              by = c("season")) %>% 
    arrange(month, day, hour) %>% 
    select(-c(month, day, hour, gea, season))
  
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_season) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{g}/", "season.rds")), compress = "gz")
  
  # Time-of-day average rate
  aer_tod <- year_hours %>%
    left_join(df_tod %>%  
                filter(gea == area) %>% 
                pivot_wider(names_from = c(scenario, year), values_from = aer), 
              by = "hour") %>% 
    arrange(month, day, hour) %>% 
    select(-c(month, day, hour, gea))
  
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_tod) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{g}/", "tod.rds")), compress = "gz")
  
  # Hourly
  n_building <- length(df_hourly)
  result <- data.frame(matrix(nrow = n_building, ncol = 0))
  for (z in 1:length(all_scenario$scenario)){
    
    scenario <- unique(df_s_hourly[[z]] %>% .$scenario)
    aer_hour <- df_s_hourly[[z]] %>% 
      mutate(toy = format(datetime, "%m-%d %H:%M:%S")) %>% 
      select(-datetime) %>% 
      filter(gea == area) %>% 
      pivot_wider(names_from = c(scenario, year), values_from = aer) %>% 
      drop_na(toy) %>% 
      select(-c(gea, toy)) 
    
    df_hourly <- df_energy %>% 
      select(timestamp, name, eload) %>% 
      mutate(eload = replace_na(eload, 0) / 1000) %>% 
      pivot_wider(id_cols = timestamp,names_from = name, values_from = eload) %>% 
      select(-timestamp) 
    
    s_result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_hour) / 1000)
    
    result <- bind_cols(result, s_result)
  }
  
  write_rds(result, str_glue(paste0(output_path, "{g}/", "hourly.rds")), compress = "gz")
  
  
}

rm(result, s_result, aer_hour, aer_tod, aer_month_hour, aer_annual)





#### ANALYSIS ####
readfile_path <- "../readfiles/results/"
figs_path <- "../figs/"

for (area in gea_map$gea_name){
  
  g <- gea_map %>% filter(gea_name == area) %>% .$gea
  
  subfigs_path <- paste0(figs_path, str_glue("{g}/"))
  
  plot_list <- list()
  z_index <- 1
  
  for (z in all_scenario$scenario){
    
    aer_hourly <- read_rds(str_glue(paste0(readfile_path, "{g}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(aer_hourly)
    
    aer_annual <- read_rds(str_glue(paste0(readfile_path, "{g}/annual.rds"))) %>% 
      select(all_of(colname))
    
    aer_tod <- read_rds(str_glue(paste0(readfile_path, "{g}/tod.rds"))) %>% 
      select(all_of(colname))
    
    aer_month_hour <- read_rds(str_glue(paste0(readfile_path, "{g}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    aer_season <- read_rds(str_glue(paste0(readfile_path, "{g}/season.rds"))) %>% 
      select(all_of(colname))
    
    aer_season_hour <- read_rds(str_glue(paste0(readfile_path, "{g}/season_hour.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate annual error distribution
    annual_err <- (aer_annual - aer_hourly) / aer_hourly * 100
    
    mean <- annual_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values)) %>% 
      ungroup() %>% 
      mutate(type = "annual", 
             gea = area)
    
    annual_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      mutate(scenario = as.factor(scenario)) %>% 
      ggplot(aes(x = values, group = scenario)) +
      geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666") +
      geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      geom_text(data = mean, 
                aes(x = -5, y = 0.2, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      facet_wrap(~scenario, nrow = 6) +
      coord_cartesian(xlim = c(-10, 10)) +
      theme(panel.grid.major.y = element_line(color = "grey80"),
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
    
    ggsave(filename = str_glue("annual_{z}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate tod error distribution
    tod_err <- (aer_tod - aer_hourly) / aer_hourly * 100
    
    mean <- tod_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values)) %>% 
      ungroup() %>% 
      mutate(type = "tod", 
             gea = area)
    
    tod_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      mutate(scenario = as.factor(scenario)) %>% 
      ggplot(aes(x = values, group = scenario)) +
      geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666") +
      geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      geom_text(data = mean, 
                aes(x = -2.5, y = 0.5, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      facet_wrap(~scenario, nrow = 6) +
      coord_cartesian(xlim = c(-5, 5)) +
      theme(panel.grid.major.y = element_line(color = "grey80"),
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
    
    ggsave(filename = str_glue("tod_{z}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate month-hour error distribution
    month_hour_err <- (aer_month_hour - aer_hourly) / aer_hourly * 100
    
    mean <- month_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values)) %>% 
      ungroup() %>% 
      mutate(type = "month-hour", 
             gea = area)
    
    month_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      mutate(scenario = as.factor(scenario)) %>% 
      ggplot(aes(x = values, group = scenario)) +
      geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666") +
      geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      geom_text(data = mean, 
                aes(x = -1, y = 0.75, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      facet_wrap(~scenario, nrow = 6) +
      coord_cartesian(xlim = c(-2.5, 2.5)) +
      theme(panel.grid.major.y = element_line(color = "grey80"),
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
    
    ggsave(filename = str_glue("month-hour_{z}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate season error distribution
    season_err <- (aer_season - aer_hourly) / aer_hourly * 100
    
    mean <- season_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values)) %>% 
      ungroup() %>% 
      mutate(type = "season", 
             gea = area)
    
    season_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      mutate(scenario = as.factor(scenario)) %>% 
      ggplot(aes(x = values, group = scenario)) +
      geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666") +
      geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      geom_text(data = mean, 
                aes(x = -5, y = 0.2, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      facet_wrap(~scenario, nrow = 6) +
      coord_cartesian(xlim = c(-2.5, 2.5)) +
      theme(panel.grid.major.y = element_line(color = "grey80"),
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
    
    ggsave(filename = str_glue("season_{z}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate season-hour error distribution
    season_hour_err <- (aer_season_hour - aer_hourly) / aer_hourly * 100
    
    mean <- season_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values)) %>% 
      ungroup() %>% 
      mutate(type = "season-hour", 
             gea = area)
    
    season_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      mutate(scenario = as.factor(scenario)) %>% 
      ggplot(aes(x = values, group = scenario)) +
      geom_histogram(bins = 100, aes(y=after_stat(density)), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666") +
      geom_vline(data = mean, aes(xintercept = mean),linewidth = 1) +
      geom_text(data = mean, 
                aes(x = -1, y = 0.75, label = paste0("mean: ", round(mean, digits = 2), " tons CO2e")), check_overlap = T) +
      facet_wrap(~scenario, nrow = 6) +
      coord_cartesian(xlim = c(-2.5, 2.5)) +
      theme(panel.grid.major.y = element_line(color = "grey80"),
            legend.direction = "horizontal",
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
    
    ggsave(filename = str_glue("season-hour_{z}_dist.png"), path = subfigs_path, units = "in", height = 8, width = 8, dpi = 300)
    
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
                         breaks = breaks_pretty(n = 4), 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = c(-40, 40)) +
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
    plot_annotation(title = str_glue("{area}"))
  
  ggsave(filename = str_glue("{area}_summary.png"), path = figs_path, units = "in", height = 8, width = 18, dpi = 300)

}
