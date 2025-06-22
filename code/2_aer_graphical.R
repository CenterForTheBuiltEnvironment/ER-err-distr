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
               "Month-hour avg." = "#84C3B7", 
               "PV" = "#feb24c", 
               "Wind" = "#9ecae1", 
               "Others" = "#045a8d", 
               "Fossil" = "grey20"
               
)

alpha_values <- c(1, 0.7, 0.4)

emissions <- "aer"


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
df_annual <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_annual.rds")) %>% 
  filter(year %in% c(2025, 2035, 2050))

df_month_hour <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_month_hour.rds")) %>% 
  filter(year %in% c(2025, 2035, 2050))

df_season_hour <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_season_hour.rds")) %>% 
  filter(year %in% c(2025, 2035, 2050))

df_season <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_season.rds")) %>% 
  filter(year %in% c(2025, 2035, 2050))

df_tod <- read_rds(paste0(readfile_path, str_glue("/{emissions}/"), "df_tod.rds")) %>% 
  filter(year %in% c(2025, 2035, 2050))

file_list <- paste0(readfile_path, str_glue("/{emissions}/"), "df_s", 1:8, "_hourly.rds")

df_hourly <- map(file_list, readRDS) %>% 
  bind_rows() %>% 
  filter(year %in% c(2025, 2035, 2050))

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





#### GENERATION ####
df_gen <- read_rds(paste0(readfile_path, "df_gen.rds")) %>% 
  filter(gea %in% c("CAISO", "ERCOT", "PJM_East"), 
         scenario %in% c("MidCase", "LowRECost", "HighNGPrice", "LowRECost_HighNGPrice"), 
         year %in% c(2025, 2035, 2050)) %>% 
  mutate(solar = rowSums(across(c("distpv_MWh", "upv_MWh"))), 
         solar_rt = round(rowSums(across(c("distpv_MWh", "upv_MWh"))) / generation * 100, digits = 0), 
         wind = rowSums(across(c("wind-ons_MWh", "wind-ofs_MWh"))), 
         wind_rt = round(rowSums(across(c("wind-ons_MWh", "wind-ofs_MWh"))) / generation * 100, digits = 0),
         foss = rowSums(across(c("coal_MWh", "gas-cc_MWh", "gas-ct_MWh"))), 
         foss_rt = round(rowSums(across(c("coal_MWh", "gas-cc_MWh", "gas-ct_MWh"))) / generation * 100, digits = 0),
         other = generation - solar - wind - foss, 
         solar_pos = wind + foss + other + solar / 2, 
         wind_pos = foss + other + wind / 2, 
         foss_pos = foss / 2) %>% 
  mutate(scenario = factor(scenario, levels = c("MidCase", "LowRECost", "HighNGPrice", "LowRECost_HighNGPrice")), 
         gea = as.factor(gea), 
         year = as.factor(year))

p1 <- df_gen %>% 
  filter(gea == "CAISO") %>% 
  select(c(gea, year, scenario, other, solar, solar_rt, wind, wind_rt, foss, foss_rt, solar_pos, wind_pos, foss_pos)) %>% 
  pivot_longer(c(solar, wind, foss, other), names_to = "gen", values_to = "value") %>% 
  mutate(gen = factor(gen, levels = c("solar", "wind", "other", "foss")), 
         gen = recode_factor(gen, 
                             "solar" = "PV", 
                             "wind" = "Wind", 
                             "other" = "Others", 
                             "foss" = "Fossil")) %>% 
  ggplot() +
  geom_col(aes(x = year, y = value / 10000000, fill = gen), position = "stack") +
  geom_text(aes(x = year, y = solar_pos / 10000000, label = paste0(solar_rt, "%")), check_overlap = T, size = 3) +
  geom_text(aes(x = year, y = wind_pos / 10000000, label = paste0(wind_rt, "%")), check_overlap = T, size = 3) +
  geom_text(aes(x = year, y = foss_pos / 10000000, label = paste0(foss_rt, "%")), check_overlap = T, color = "white", size = 3) +
  facet_wrap(~scenario, nrow = 1) +
  labs(x = NULL,
       y = NULL,
       fill = NULL,
       subtitle = "CAISO") +
  scale_fill_manual(values = ls_colors) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- df_gen %>% 
  filter(gea == "ERCOT") %>% 
  select(c(gea, year, scenario, other, solar, solar_rt, wind, wind_rt, foss, foss_rt, solar_pos, wind_pos, foss_pos)) %>% 
  pivot_longer(c(solar, wind, foss, other), names_to = "gen", values_to = "value") %>% 
  mutate(gen = factor(gen, levels = c("solar", "wind", "other", "foss")), 
         gen = recode_factor(gen, 
                             "solar" = "PV", 
                             "wind" = "Wind", 
                             "other" = "Others", 
                             "foss" = "Fossil")) %>% 
  ggplot() +
  geom_col(aes(x = year, y = value / 10000000, fill = gen), position = "stack") +
  geom_text(aes(x = year, y = solar_pos / 10000000, label = paste0(solar_rt, "%")), check_overlap = T, size = 3) +
  geom_text(aes(x = year, y = wind_pos / 10000000, label = paste0(wind_rt, "%")), check_overlap = T, size = 3) +
  geom_text(aes(x = year, y = foss_pos / 10000000, label = paste0(foss_rt, "%")), check_overlap = T, color = "white", size = 3, position = position_nudge(y = 0.5)) +
  facet_wrap(~scenario, nrow = 1) +
  labs(x = NULL,
       y = NULL,
       fill = NULL,
       subtitle = "ERCOT") +
  scale_fill_manual(values = ls_colors) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        axis.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p3 <- df_gen %>% 
  filter(gea == "PJM_East") %>% 
  select(c(gea, year, scenario, other, solar, solar_rt, wind, wind_rt, foss, foss_rt, solar_pos, wind_pos, foss_pos)) %>% 
  pivot_longer(c(solar, wind, foss, other), names_to = "gen", values_to = "value") %>% 
  mutate(gen = factor(gen, levels = c("solar", "wind", "other", "foss")), 
         gen = recode_factor(gen, 
                             "solar" = "PV", 
                             "wind" = "Wind", 
                             "other" = "Others", 
                             "foss" = "Fossil")) %>% 
  ggplot() +
  geom_col(aes(x = year, y = value / 10000000, fill = gen), position = "stack") +
  geom_text(aes(x = year, y = solar_pos / 10000000, label = paste0(solar_rt, "%")), check_overlap = T, size = 3, position = position_nudge(y = 0.25)) +
  geom_text(aes(x = year, y = wind_pos / 10000000, label = paste0(wind_rt, "%")), check_overlap = T, size = 3, position = position_nudge(y = -0.25)) +
  geom_text(aes(x = year, y = foss_pos / 10000000, label = paste0(foss_rt, "%")), check_overlap = T, color = "white", size = 3, position = position_nudge(y = 0.5)) +
  facet_wrap(~scenario, nrow = 1) +
  labs(x = NULL,
       y = NULL,
       fill = NULL,
       subtitle = "PJM_East") +
  scale_fill_manual(values = ls_colors) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3, 
                        nrow = 3, 
                        common.legend = TRUE, 
                        legend = "bottom")

annotate_figure(final_plot, 
                left = text_grob("Energy Generation (×10 million MWh)", 
                                 rot = 90, vjust = 1, size = 12),
                top = text_grob("Energy generation by different sources from the grid", size = 14))

ggsave(filename = "gen.png", path = figs_path, units = "in", height = 8, width = 8, dpi = 300)





#### CAISO ####
gea_example <- "CAISO"
sce_example <- "MidCase"

df_annual %>% 
  filter(gea == gea_example, 
         scenario == sce_example) %>% 
  mutate(type = as.factor("Annual avg."), 
         year = factor(year)) %>% 
  ggplot() +
  geom_col(aes(x = year, y = er, fill = type, alpha = year)) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_fill_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 300, by = 50)) +
  coord_cartesian(ylim = c(0, 300)) +
  labs(x = NULL,
       y = "Emissions rate",
       fill = NULL,
       title = "Annual averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("{emissions}_annual.svg"), path = paste0(figs_path, str_glue("{gea_example}/")), units = "in", height = 8, width = 8, dpi = 300)

df_season %>% 
  filter(gea == gea_example, 
         scenario == sce_example) %>% 
  mutate(type = as.factor("Season avg."), 
         year = as.factor(year), 
         season = factor(season, levels = c("spring", "summer", "fall", "winter"))) %>% 
  ggplot() +
  geom_col(aes(x = season, y = er, fill = type, alpha = year), position = "dodge") +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_fill_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 300, by = 50)) +
  coord_cartesian(ylim = c(0, 300)) +
  labs(x = NULL,
       y = "Emissions rate",
       fill = NULL,
       title = "Season averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("{emissions}_seaons.svg"), path = paste0(figs_path, str_glue("{gea_example}/")), units = "in", height = 8, width = 8, dpi = 300)

df_tod %>% 
  filter(gea == gea_example, 
         scenario == sce_example) %>% 
  mutate(type = as.factor("Time-of-day avg."), 
         year = as.factor(year)) %>% 
  ggplot() +
  geom_line(aes(x = hour, y = er, color = type, alpha = year), linewidth = 3) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 300, by = 50)) +
  coord_cartesian(ylim = c(0, 300)) +
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

ggsave(filename = str_glue("{emissions}_tod.svg"), path = paste0(figs_path, str_glue("{gea_example}/")), units = "in", height = 8, width = 8, dpi = 300)

df_season_hour %>% 
  filter(gea == gea_example, 
         scenario == sce_example) %>% 
  mutate(type = as.factor("Season-hour avg."), 
         year = as.factor(year), 
         season = factor(season, levels = c("spring", "summer", "fall", "winter"))) %>% 
  ggplot() +
  geom_line(aes(x = hour, y = er, color = type, alpha = year), linewidth = 2) +
  facet_wrap(~season) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 400, by = 100)) +
  coord_cartesian(ylim = c(0, 400)) +
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

ggsave(filename = str_glue("{emissions}_season_hour.svg"), path = paste0(figs_path, str_glue("{gea_example}/")), units = "in", height = 8, width = 12, dpi = 300)

df_month_hour %>% 
  filter(gea == gea_example, 
         scenario == sce_example) %>% 
  mutate(type = as.factor("Month-hour avg."), 
         year = as.factor(year), 
         month = factor(month.abb[month],levels = month.abb)) %>% 
  ggplot() +
  geom_line(aes(x = hour, y = er, color = type, alpha = year), linewidth = 1.5) +
  facet_wrap(~month, nrow = 2) +
  scale_alpha_manual(name = "Year", values = alpha_values) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 450, by = 100)) +
  coord_cartesian(ylim = c(0, 450)) +
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

ggsave(filename = str_glue("{emissions}_month_hour.svg"), path = paste0(figs_path, str_glue("{gea_example}/")), units = "in", height = 8, width = 18, dpi = 300)

df_hourly %>% 
  filter(gea == gea_example, 
         scenario == sce_example) %>% 
  mutate(type = as.factor("Hourly avg."), 
         year = as.factor(year)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = er, color = type, group = year), size = 0.1, alpha = 0.2) +
  geom_smooth(aes(x = datetime, y = er, color = type, group = year), linewidth = 0.8) +
  facet_wrap(~year, nrow = 3, scales = "free") +
  scale_color_manual(values = ls_colors) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = breaks_pretty(n = 5)) +
  labs(x = NULL,
       y = "Emissions rate",
       fill = NULL,
       title = "Time-of-day averaged emissions rate") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("{emissions}_hour.svg"), path = paste0(figs_path, str_glue("{gea_example}/")), units = "in", height = 8, width = 8, dpi = 300)

  
# all scenarios at hourly resolution at 2050
p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_smooth(aes(x = datetime, y = er), linewidth = 0.8, alpha = 0.4) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 400, by = 100), 
                     labels = c("0", "100", "200", "300", "400")) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_smooth(aes(x = datetime, y = er), linewidth = 0.8, alpha = 0.4) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 150, by = 50), 
                     labels = c("0", "50", "100", "150")) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_smooth(aes(x = datetime, y = er), linewidth = 0.8, alpha = 0.4) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 150, by = 50), 
                     labels = c("0", "50", "100", "150")) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3, 
                        nrow = 3, 
                        common.legend = T, 
                        legend = "bottom")


annotate_figure(final_plot, 
                left = text_grob("Emissions rate (gCO2e/kWh)", 
                                 rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly emissions rate at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 8, width = 8, dpi = 300)

# Summary
sce_example <- c("MidCase", "LowRECost", "HighNGPrice", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in seq(0, 100, by = 50)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
        mutate(type = "Month-hour avg.")) %>% 
      filter(year %in% c(2025, 2035, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 40, by = 10), 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = c(0, 50)) +
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

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 10, width = 16, dpi = 300)


# avoided
# only annual and season
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate annual error distribution
    annual_err <- abs((er_annual - er_hourly) / er_hourly * 100)
    
    # Calculate season error distribution
    season_err <- abs((er_season - er_hourly) / er_hourly * 100)

    
    # Combine plots for all
    df_error <- bind_rows(
      annual_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "-") %>% 
        mutate(type = "Annual avg."), 
      season_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "-") %>% 
        mutate(type = "Season avg.")) %>% 
      filter(year %in% c(2025, 2035, 2050))
    
    breaks <- if (perc == 25) seq(0, 12000, by = 2500) else seq(0, 350, by = 50)
    ylim <- if (perc == 25) c(0, 12000) else c(0, 350)
      
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = breaks, 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) +
      geom_hline(yintercept = 100, color = "red", lty = "dashed")
    
    
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

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided_sep.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 6, width = 12, dpi = 300)

# without annual and season
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate tod error distribution
    tod_err <- abs((er_tod - er_hourly) / er_hourly * 100)
    
    # Calculate month-hour error distribution
    month_hour_err <- abs((er_month_hour - er_hourly) / er_hourly * 100)
    
    # Calculate season-hour error distribution
    season_hour_err <- abs((er_season_hour - er_hourly) / er_hourly * 100)
    
    # Combine plots for all
    df_error <- bind_rows(
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
        mutate(type = "Month-hour avg.")) %>% 
      filter(year %in% c(2025, 2035, 2050))
    
    breaks <- if (perc == 25) seq(0, 2500, by = 500) else seq(0, 60, by = 20)
    ylim <- if (perc == 25) c(0, 2800) else c(0, 60)
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = breaks, 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) + 
      geom_hline(yintercept = 100, color = "red", lty = "dashed")
    
    
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

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 6, width = 12, dpi = 300)





#### ERCOT ####
gea_example <- "ERCOT"

# all scenarios at hourly resolution at 2050
p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_smooth(aes(x = datetime, y = er), linewidth = 0.8, alpha = 0.4) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 500, by = 100), 
                     labels = c("0", "100", "200", "300", "400", "500")) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_smooth(aes(x = datetime, y = er), linewidth = 0.8, alpha = 0.4) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 300, by = 100), 
                     labels = c("0", "100", "200", "300")) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_smooth(aes(x = datetime, y = er), linewidth = 0.8, alpha = 0.4) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 300, by = 100), 
                     labels = c("0", "100", "200", "300")) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3, 
                        nrow = 3, 
                        common.legend = T, 
                        legend = "bottom")


annotate_figure(final_plot, 
                left = text_grob("Emissions rate (gCO2e/kWh)", 
                                 rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly emissions rate at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 8, width = 8, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost", "HighNGPrice", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in seq(0, 100, by = 50)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
        mutate(type = "Month-hour avg.")) %>% 
      filter(year %in% c(2025, 2035, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 40, by = 10), 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = c(0, 45)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) +
      geom_hline(yintercept = 100, color = "red", lty = "dashed")
    
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

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 10, width = 16, dpi = 300)


# avoided
# only annual and season
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate annual error distribution
    annual_err <- abs((er_annual - er_hourly) / er_hourly * 100)
    
    
    # Calculate season error distribution
    season_err <- abs((er_season - er_hourly) / er_hourly * 100)

    
    # Combine plots for all
    df_error <- bind_rows(
      annual_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "-") %>% 
        mutate(type = "Annual avg."), 
      season_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "-") %>% 
        mutate(type = "Season avg.")) %>% 
      filter(year %in% c(2025, 2035, 2050))
    
    breaks <- if (perc == 25) seq(0, 1000, by = 250) else seq(0, 300, by = 50)
    ylim <- if (perc == 25) c(0, 1200) else c(0, 300)
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = breaks, 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) +
      geom_hline(yintercept = 100, color = "red", lty = "dashed")
    
    
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

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided_sep.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 6, width = 12, dpi = 300)

# without annual and season
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate tod error distribution
    tod_err <- abs((er_tod - er_hourly) / er_hourly * 100)
    
    # Calculate month-hour error distribution
    month_hour_err <- abs((er_month_hour - er_hourly) / er_hourly * 100)
    
    # Calculate season-hour error distribution
    season_hour_err <- abs((er_season_hour - er_hourly) / er_hourly * 100)
    
    # Combine plots for all
    df_error <- bind_rows(
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
        mutate(type = "Month-hour avg.")) %>% 
      filter(year %in% c(2025, 2035, 2050))
    
    breaks <- if (perc == 25) seq(0, 200, by = 50) else seq(0, 15, by = 5)
    ylim <- if (perc == 25) c(0, 220) else c(0, 18)
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = breaks, 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) +
      geom_hline(yintercept = 100, color = "red", lty = "dashed")
    
    
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

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 6, width = 12, dpi = 300)





#### PJM EAST ####
gea_example <- "PJM_East"

# all scenarios at hourly resolution at 2050
p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_smooth(aes(x = datetime, y = er), linewidth = 0.8, alpha = 0.4) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 800, by = 200), 
                     labels = c("0", "200", "400", "600", "800")) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_smooth(aes(x = datetime, y = er), linewidth = 0.8, alpha = 0.4) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 400, by = 100), 
                     labels = c("0", "100", "200", "300", "400")) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_smooth(aes(x = datetime, y = er), linewidth = 0.8, alpha = 0.4) +
  scale_x_datetime(labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 400, by = 100), 
                     labels = c("0", "100", "200", "300", "400")) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3, 
                        nrow = 3, 
                        common.legend = T, 
                        legend = "bottom")


annotate_figure(final_plot, 
                left = text_grob("Emissions rate (gCO2e/kWh)", 
                                 rot = 90, vjust = 1, size = 12),
                top = text_grob("Hourly emissions rate", size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 8, width = 8, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost", "HighNGPrice", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in seq(0, 100, by = 50)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/operational/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
        mutate(type = "Month-hour avg.")) %>% 
      filter(year %in% c(2025, 2035, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 25, by = 5), 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = c(0, 25)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) +
      geom_hline(yintercept = 100, color = "red", lty = "dashed")
    
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

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 10, width = 16, dpi = 300)


# avoided
# only annual and season
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate annual error distribution
    annual_err <- abs((er_annual - er_hourly) / er_hourly * 100)
    
    
    # Calculate season error distribution
    season_err <- abs((er_season - er_hourly) / er_hourly * 100)
    
    
    # Combine plots for all
    df_error <- bind_rows(
      annual_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "-") %>% 
        mutate(type = "Annual avg."), 
      season_err %>% 
        pivot_longer(everything(), names_to = "year", values_to = "error") %>% 
        separate(year, into = c("scenario", "year"), sep = "-") %>% 
        mutate(type = "Season avg.")) %>% 
      filter(year %in% c(2025, 2035, 2050))
    
    breaks <- if (perc == 25) seq(0, 200, by = 50) else seq(0, 75, by = 25)
    ylim <- if (perc == 25) c(0, 200) else c(0, 75)
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = breaks, 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) +
      geom_hline(yintercept = 100, color = "red", lty = "dashed")
    
    
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

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided_sep.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 6, width = 12, dpi = 300)

# without annual and season
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
      select(all_of(colname))
    
    # Calculate tod error distribution
    tod_err <- abs((er_tod - er_hourly) / er_hourly * 100)
    
    # Calculate month-hour error distribution
    month_hour_err <- abs((er_month_hour - er_hourly) / er_hourly * 100)
    
    # Calculate season-hour error distribution
    season_hour_err <- abs((er_season_hour - er_hourly) / er_hourly * 100)
    
    # Combine plots for all
    df_error <- bind_rows(
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
        mutate(type = "Month-hour avg.")) %>% 
      filter(year %in% c(2025, 2035, 2050))
    
    breaks <- if (perc == 25) seq(0, 125, by = 25) else seq(0, 30, by = 10)
    ylim <- if (perc == 25) c(0, 125) else c(0, 35)
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error, fill = type)) +
      geom_lv(alpha = 0.4, k = 4, outlier.size = 0.4) +
      geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = type)) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = breaks, 
                         labels = number_format(suffix = " %")) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) +
      geom_hline(yintercept = 100, color = "red", lty = "dashed")
    
    
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

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = paste0(figs_path, str_glue("{gea_example}/{emissions}")), units = "in", height = 6, width = 12, dpi = 300)
