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
site_map <- read_csv(paste0(readfile_path, "site_map.csv"))

df_energy <- read_rds(paste0(readfile_path, "df_energy.rds")) %>% 
  left_join(site_map, by = "site") %>% 
  filter(year(timestamp) == 2016) %>% 
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

df_season_hour <- read_rds(paste0(readfile_path, "df_season_hour.rds"))

df_season <- read_rds(paste0(readfile_path, "df_season.rds"))

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
  
  df_hourly <- df_energy %>% 
    select(timestamp, name, eload) %>% 
    mutate(eload = replace_na(eload, 0) / 1000) %>% 
    pivot_wider(id_cols = timestamp,names_from = name, values_from = eload) %>% 
    select(-timestamp) 
  
  # Annual average rate
  aer_annual <- df_annual %>% 
    filter(gea == area) %>% 
    pivot_wider(names_from = c(scenario, year), values_from = aer) %>% 
    slice(rep(1, 8784)) %>% 
    select(-gea)
  
  result <- as.data.frame(t(as.matrix(df_hourly)) %*% as.matrix(aer_annual) / 1000)
  
  write_rds(result, str_glue(paste0(output_path, "{g}/", "annual.rds")), compress = "gz")
  
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
      filter(date(timestamp) != as.Date("2016-02-29")) %>% 
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

all_mean <- data.frame()
all_ci <- data.frame()

for (area in gea_map$gea_name){
  g <- gea_map %>% filter(gea_name == area) %>% .$gea
  
  figs_path <- str_glue("../figs/{g}/")
  
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
    
    all_mean <- bind_rows(all_mean, mean)
    
    ci <- annual_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values, na.rm = TRUE),
                se = sd(values, na.rm = TRUE) / sqrt(n()),
                low = mean - 1.96 * se,
                high = mean + 1.96 * se) %>% 
      ungroup() %>% 
      mutate(type = "annual", 
             gea = area)
    
    all_ci <- bind_rows(all_ci, ci)
    
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
    
    
    ggsave(filename = str_glue("annual_{z}_dist.png"), path = str_glue(figs_path), units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate tod error distribution
    tod_err <- (aer_tod - aer_hourly) / aer_hourly * 100

    mean <- tod_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values)) %>% 
      ungroup() %>% 
      mutate(type = "tod", 
             gea = area)
    
    all_mean <- bind_rows(all_mean, mean)
    
    ci <- tod_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values, na.rm = TRUE),
                se = sd(values, na.rm = TRUE) / sqrt(n()),
                low = mean - 1.96 * se,
                high = mean + 1.96 * se) %>% 
      ungroup() %>% 
      mutate(type = "tod", 
             gea = area)
    
    all_ci <- bind_rows(all_ci, ci)
    
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
    
    
    ggsave(filename = str_glue("tod_{z}_dist.png"), path = str_glue(figs_path), units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate month-hour error distribution
    month_hour_err <- (aer_hourly - aer_month_hour) / aer_hourly * 100
    
    mean <- month_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values)) %>% 
      ungroup() %>% 
      mutate(type = "month-hour", 
             gea = area)
    
    all_mean <- bind_rows(all_mean, mean)
    
    ci <- month_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values, na.rm = TRUE),
                se = sd(values, na.rm = TRUE) / sqrt(n()),
                low = mean - 1.96 * se,
                high = mean + 1.96 * se) %>% 
      ungroup() %>% 
      mutate(type = "month-hour", 
             gea = area)
    
    all_ci <- bind_rows(all_ci, ci)
    
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
    
    
    ggsave(filename = str_glue("month-hour_{z}_dist.png"), path = str_glue(figs_path), units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate season error distribution
    season_err <- (aer_hourly - aer_season) / aer_hourly * 100
    
    mean <- season_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values)) %>% 
      ungroup() %>% 
      mutate(type = "season", 
             gea = area)
    
    all_mean <- bind_rows(all_mean, mean)
    
    ci <- season_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values, na.rm = TRUE),
                se = sd(values, na.rm = TRUE) / sqrt(n()),
                low = mean - 1.96 * se,
                high = mean + 1.96 * se) %>% 
      ungroup() %>% 
      mutate(type = "season", 
             gea = area)
    
    all_ci <- bind_rows(all_ci, ci)
    
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
    
    ggsave(filename = str_glue("season_{z}_dist.png"), path = str_glue(figs_path), units = "in", height = 8, width = 8, dpi = 300)
    
    # Calculate season-hour error distribution
    season_hour_err <- (aer_hourly - aer_season_hour) / aer_hourly * 100
    
    mean <- season_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values)) %>% 
      ungroup() %>% 
      mutate(type = "season-hour", 
             gea = area)
    
    all_mean <- bind_rows(all_mean, mean)
    
    ci <- season_hour_err %>% 
      pivot_longer(cols = everything(), names_to = "scenario", values_to = "values") %>% 
      group_by(scenario) %>% 
      summarise(mean = mean(values, na.rm = TRUE),
                se = sd(values, na.rm = TRUE) / sqrt(n()),
                low = mean - 1.96 * se,
                high = mean + 1.96 * se) %>% 
      ungroup() %>% 
      mutate(type = "season-hour", 
             gea = area)
    
    all_ci <- bind_rows(all_ci, ci)
    
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
    
    
    ggsave(filename = str_glue("season-hour_{z}_dist.png"), path = str_glue(figs_path), units = "in", height = 8, width = 8, dpi = 300)
    
  }

}

# Overall plot for distribution mean
all_mean <- bind_rows(all_mean)

for (area in gea_map$gea_name){
  
  g <- gea_map %>% filter(gea_name == area) %>% .$gea
  
  all_mean %>% 
    filter(gea == area) %>% 
    mutate(type = as_factor(type), 
           type = recode_factor(type, 
                                "annual" = "Annual avg.", 
                                "season" = "Season avg.", 
                                "tod" = "Time-of-day avg.", 
                                "season-hour" = "Season-hour avg.",
                                "month-hour" = "Month-hour avg."
                                )) %>% 
    separate(scenario, into = c("scenario", "year"), sep = "_") %>% 
    filter(year %in% c(2025, 2030, 2050)) %>% 
    ggplot() +
    geom_bar(aes(x = year, y = mean, fill = type), stat = "identity", position = "dodge") +
    geom_vline(xintercept = seq(1.5, 5, by = 1), lty = "dashed", color = "grey80") +
    facet_wrap(~ scenario, nrow = 2) +
    labs(x = "Projected year", 
         y = "Error distribution mean") +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " %")) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = str_glue("{area}"), 
         fill = NULL) +
    theme(panel.grid.major.y = element_line(color = "grey80"),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggsave(filename = str_glue("{g}_mean_summary.png"), path = str_glue("../figs/"), units = "in", height = 6, width = 12, dpi = 300)
}

# Overall plot for distribution confidence interval
all_ci <- bind_rows(all_ci)

for (area in gea_map$gea_name){
  
  g <- gea_map %>% filter(gea_name == area) %>% .$gea
  
  all_ci %>% 
    filter(gea == area) %>% 
    mutate(type = as_factor(type), 
           type = recode_factor(type, 
                                "annual" = "Annual avg.", 
                                "season" = "Season avg.", 
                                "tod" = "Time-of-day avg.", 
                                "season-hour" = "Season-hour avg.",
                                "month-hour" = "Month-hour avg.")) %>% 
    separate(scenario, into = c("scenario", "year"), sep = "_") %>% 
    filter(year == 2025 | year == 2050) %>% 
    ggplot() +
    geom_point(aes(x = year, y = mean, color = type), 
               size = 2, 
               position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(x = year, ymin = low, ymax = high, color = type), 
                  width = 0.75, 
                  alpha = 0.8, 
                  position = position_dodge(width = 0.75)) +
    geom_vline(xintercept = 1.5, lty = "dashed", color = "grey80") +
    facet_wrap(~ scenario, nrow = 2) +
    labs(x = "Projected year", 
         y = "Error distribution mean and 95% CI") +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " %")) +
    scale_color_brewer(palette = "Set2") +
    labs(title = str_glue("{area}"), 
         color = NULL) +
    theme(panel.grid.major.y = element_line(color = "grey80"),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggsave(filename = str_glue("{g}_ci_summary.png"), path = str_glue("../figs/"), units = "in", height = 6, width = 12, dpi = 300)
}




#### CLUSTER ####
# k <- 6
# df_energyT <- df_energy %>% 
#   select(name, timestamp, eload) %>% 
#   arrange(name, timestamp) %>% 
#   group_by(name) %>% 
#   mutate(across(eload, ~ zoo::na.approx(., na.rm = FALSE))) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = timestamp, values_from = eload) %>% 
#   select(-name) %>% 
#   drop_na()
# 
# kmeans_result <- kmeans(df_energyT, centers = k, nstart = 25)
# 
# # Plot the clusters
# df_cluster <- as.data.frame(kmeans_result$centers) %>%
#   mutate(cluster = as.factor(row_number())) %>% 
#   pivot_longer(-cluster, names_to = "datetime", values_to = "eload") 
# 
# df_cluster %>% 
#   mutate(datetime = as.POSIXct(datetime)) %>% 
#   ggplot(aes(x = datetime, y = eload, group = cluster, color = cluster)) +
#   geom_smooth() +
#   labs(title = "K-means Clustering of Genome Dataset",
#        x = "Datetime",
#        y = "Eload") +
#   scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
#   theme(panel.grid.major.y = element_line(color = "grey80"),
#         legend.direction = "horizontal",
#         legend.position = "bottom",
#         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
