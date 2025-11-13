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
theme_update(plot.title = element_text(size = 16, colour = "grey20", face = "bold", hjust = 0.5),
             plot.subtitle = element_text(size = 13, colour = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 15)),
             plot.caption = element_text(size = 10, colour = "grey20", face = "italic", hjust = 0.5),
             plot.background = element_rect(fill = "white", colour = NA),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             axis.text = element_text(size = 12),
             legend.text = element_text(size = 12),
             strip.text = element_text(size = 12, color = "grey20", face = "bold"),
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
figs_path <- "../manuscript/figs/"

solar_map <- read_csv(paste0(readfile_path, "solar_map.csv"))

df_energy <- read_rds(paste0(readfile_path, "df_energy.rds")) %>% 
  filter(year(timestamp) == 2016, 
         date(timestamp) != "2016-02-29") %>% 
  select(site, 
         timestamp, 
         type, 
         name, 
         eload)

df_compr <- read_rds(paste0(readfile_path, "df_energy.rds")) %>% 
  filter(year(timestamp) == 2017) %>% 
  select(site, 
         timestamp, 
         type, 
         name, 
         eload)

df_meta <- read_rds(paste0(readfile_path, "df_meta.rds"))

df_weather <- read_rds(paste0(readfile_path, "/weather/", "df_weather.rds")) %>% 
  filter(year(timestamp) == 2016)

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

all_scenario <- df_annual %>% 
  select(scenario) %>% 
  distinct()

all_year <- df_annual %>% 
  select(year) %>% 
  distinct()

all_region <- unique(df_annual %>% .$gea)





#### GENERATION ####
df_gen <- read_rds(paste0(readfile_path, "df_gen.rds")) %>% 
  filter(gea %in% c("CAISO", "ERCOT", "PJM_East"), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice"), 
         year %in% c(2025, 2050)) %>% 
  mutate(solar = rowSums(across(c("distpv_MWh", "upv_MWh"))), 
         solar_rt = round(rowSums(across(c("distpv_MWh", "upv_MWh"))) / generation * 100, digits = 0), 
         wind = rowSums(across(c("wind-ons_MWh", "wind-ofs_MWh"))), 
         wind_rt = round(rowSums(across(c("wind-ons_MWh", "wind-ofs_MWh"))) / generation * 100, digits = 0),
         foss = rowSums(across(c("coal_MWh", "gas-cc_MWh", "gas-ct_MWh"))), 
         foss_rt = round(rowSums(across(c("coal_MWh", "gas-cc_MWh", "gas-ct_MWh"))) / generation * 100, digits = 0),
         other = generation - solar - wind - foss, 
         other_rt = round(other / generation * 100, digits = 0), 
         solar_pos = wind + foss + other + solar / 2, 
         wind_pos = foss + other + wind / 2, 
         other_pos = foss + other / 2,
         foss_pos = foss / 2) %>% 
  mutate(scenario = factor(scenario, levels = c("MidCase", "LowRECost_HighNGPrice")), 
         gea = as.factor(gea), 
         year = as.factor(year))

p1 <- df_gen %>% 
  filter(gea == "CAISO") %>% 
  select(c(gea, year, scenario, other, other_rt, solar, solar_rt, wind, wind_rt, foss, foss_rt, other_pos, solar_pos, wind_pos, foss_pos)) %>% 
  pivot_longer(c(solar, wind, foss, other), names_to = "gen", values_to = "value") %>% 
  mutate(gen = factor(gen, levels = c("solar", "wind", "other", "foss")), 
         gen = recode_factor(gen, 
                             "solar" = "PV", 
                             "wind" = "Wind", 
                             "other" = "Others", 
                             "foss" = "Fossil")) %>% 
  ggplot() +
  geom_col(aes(x = year, y = value / 10000000, fill = gen), position = "stack") +
  geom_text(aes(x = year, y = solar_pos / 10000000, label = ifelse(solar_rt < 10, NA, paste0(solar_rt, "%"))), check_overlap = T, size = 4) +
  geom_text(aes(x = year, y = wind_pos / 10000000, label = ifelse(wind_rt < 10, NA, paste0(wind_rt, "%"))), check_overlap = T, size = 4) +
  geom_text(aes(x = year, y = other_pos / 10000000, label = ifelse(other_rt < 10, NA, paste0(other_rt, "%"))), check_overlap = T, color = "white", size = 4) +
  geom_text(aes(x = year, y = foss_pos / 10000000, label = ifelse(foss_rt < 10, NA, paste0(foss_rt, "%"))), check_overlap = T, color = "white", size = 4) +
  facet_wrap(~scenario, nrow = 1) +
  scale_y_continuous(breaks = seq(0, 100, by = 25)) +
  labs(x = NULL,
       y = NULL,
       fill = NULL,
       subtitle = "CAISO") +
  scale_fill_manual(values = ls_colors) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- df_gen %>% 
  filter(gea == "ERCOT") %>% 
  select(c(gea, year, scenario, other, other_rt, solar, solar_rt, wind, wind_rt, foss, other_pos, foss_rt, solar_pos, wind_pos, foss_pos)) %>% 
  pivot_longer(c(solar, wind, foss, other), names_to = "gen", values_to = "value") %>% 
  mutate(gen = factor(gen, levels = c("solar", "wind", "other", "foss")), 
         gen = recode_factor(gen, 
                             "solar" = "PV", 
                             "wind" = "Wind", 
                             "other" = "Others", 
                             "foss" = "Fossil")) %>% 
  ggplot() +
  geom_col(aes(x = year, y = value / 10000000, fill = gen), position = "stack") +
  geom_text(aes(x = year, y = solar_pos / 10000000, label = ifelse(solar_rt < 10, NA, paste0(solar_rt, "%"))), check_overlap = T, size = 4) +
  geom_text(aes(x = year, y = wind_pos / 10000000, label = ifelse(wind_rt < 10, NA, paste0(wind_rt, "%"))), check_overlap = T, size = 4) +
  geom_text(aes(x = year, y = other_pos / 10000000, label = ifelse(other_rt < 10, NA, paste0(other_rt, "%"))), check_overlap = T, color = "white", size = 4) +
  geom_text(aes(x = year, y = foss_pos / 10000000, label = ifelse(foss_rt < 10, NA, paste0(foss_rt, "%"))), check_overlap = T, color = "white", size = 4) +
  facet_wrap(~scenario, nrow = 1) +
  scale_y_continuous(breaks = seq(0, 100, by = 25)) +
  labs(x = NULL,
       y = NULL,
       fill = NULL,
       subtitle = "ERCOT") +
  scale_fill_manual(values = ls_colors) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        axis.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p3 <- df_gen %>% 
  filter(gea == "PJM_East") %>% 
  select(c(gea, year, scenario, other, other_rt, solar, solar_rt, wind, wind_rt, foss, foss_rt, other_pos, solar_pos, wind_pos, foss_pos)) %>% 
  pivot_longer(c(solar, wind, foss, other), names_to = "gen", values_to = "value") %>% 
  mutate(gen = factor(gen, levels = c("solar", "wind", "other", "foss")), 
         gen = recode_factor(gen, 
                             "solar" = "PV", 
                             "wind" = "Wind", 
                             "other" = "Others", 
                             "foss" = "Fossil")) %>% 
  ggplot() +
  geom_col(aes(x = year, y = value / 10000000, fill = gen), position = "stack") +
  geom_text(aes(x = year, y = solar_pos / 10000000, label = ifelse(solar_rt < 10, NA, paste0(solar_rt, "%"))), check_overlap = T, size = 4) +
  geom_text(aes(x = year, y = wind_pos / 10000000, label = ifelse(wind_rt < 10, NA, paste0(wind_rt, "%"))), check_overlap = T, size = 4) +
  geom_text(aes(x = year, y = other_pos / 10000000, label = ifelse(other_rt < 10, NA, paste0(other_rt, "%"))), check_overlap = T, color = "white", size = 4) +
  geom_text(aes(x = year, y = foss_pos / 10000000, label = ifelse(foss_rt < 10, NA, paste0(foss_rt, "%"))), check_overlap = T, color = "white", size = 4) +
  facet_wrap(~scenario, nrow = 1) +
  scale_y_continuous(breaks = seq(0, 100, by = 25)) +
  labs(x = NULL,
       y = NULL,
       fill = NULL,
       subtitle = "PJM_East") +
  scale_fill_manual(values = ls_colors) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3, 
                        nrow = 3, 
                        heights = c(1, 1, 1.2), 
                        common.legend = TRUE, 
                        legend = "bottom")

y_axis <- expression(Grid ~ electricity ~ generation ~ (10^7 ~ MWh))
annotate_figure(final_plot, 
                left = text_grob(y_axis, 
                                 rot = 90, vjust = 1, size = 12),
                top = text_grob("Energy generation by different sources from the grid", size = 14))

ggsave(filename = "gen.png", path = figs_path, units = "in", height = 8, width = 6, dpi = 300)





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
       y = "Emission rates",
       fill = NULL,
       title = "Annual averaged emission rates") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("{gea_example}_{emissions}_annual.svg"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)

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
       y = "Emission rates",
       fill = NULL,
       title = "Season averaged emission rates") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("{gea_example}_{emissions}_seaons.svg"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)

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
       y = "Emission rates",
       color = NULL,
       title = "Time-of-day averaged emission rates") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("{gea_example}_{emissions}_tod.svg"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)

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
       y = "Emission rates",
       color = NULL,
       title = "Season-hour averaged emission rates") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("{gea_example}_{emissions}_season_hour.svg"), path = figs_path, units = "in", height = 8, width = 12, dpi = 300)

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
       y = "Emission rates",
       color = NULL,
       title = "Month-hour averaged emission rates") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("{gea_example}_{emissions}_month_hour.svg"), path = figs_path, units = "in", height = 8, width = 18, dpi = 300)

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
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = breaks_pretty(n = 5)) +
  labs(x = NULL,
       y = "Emission rates",
       fill = NULL,
       title = "Time-of-day averaged emission rates") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("{gea_example}_{emissions}_hour.svg"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)

  
# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), 
                   labels = date_format("%b"), 
                   date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 400, by = 100), 
                     labels = c("0", "100", "200", "300", "400")) + 
  coord_cartesian(ylim = c(0, 420)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), 
                   labels = date_format("%b"), 
                   date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 150, by = 50), 
                     labels = c("0", "50", "100", "150")) +
  coord_cartesian(ylim = c(0, 160)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        axis.text.x = element_blank(), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.4) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er),
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), 
                   labels = date_format("%b"), 
                   date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 150, by = 50), 
                     labels = c("0", "50", "100", "150")) +
  coord_cartesian(ylim = c(0, 160)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# Summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 38, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 40, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 42)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
        
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
    
      
    } else {
        
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(1, 10, 100, 1000, 10000)
    ylim <- if (perc == 25) c(1, 40000) else c(1, 800)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 20000, 500), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
      
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}
  
ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### ERCOT ####
gea_example <- "ERCOT"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 500, by = 100), 
                     labels = c("0", "100", "200", "300", "400", "500")) +
  coord_cartesian(ylim = c(0, 550)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed",
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), 
                   labels = date_format("%b"), 
                   date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 300, by = 100), 
                     labels = c("0", "100", "200", "300")) +
  coord_cartesian(ylim = c(0, 350)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), 
                   labels = date_format("%b"), 
                   date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 300, by = 100), 
                     labels = c("0", "100", "200", "300")) +
  coord_cartesian(ylim = c(0, 350)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 38, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 40, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 42)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)
# ggsave(filename = str_glue("{gea_example}_{emissions}_operational.svg"), units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000)
    ylim <- if (perc == 25) c(0.1, 4000) else c(0.1, 800)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 2000, 400), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)






#### PJM EAST ####
gea_example <- "PJM_East"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 800, by = 200), 
                     labels = c("0", "200", "400", "600", "800")) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 650)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 400, by = 100), 
                     labels = c("0", "100", "200", "300", "400")) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 450)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 400, by = 100), 
                     labels = c("0", "100", "200", "300", "400")) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 420)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 18, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 20, by = 5), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 21)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100)
    ylim <- if (perc == 25) c(0.1, 300) else c(0.1, 300)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 180, 150), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### FRCC ####
gea_example <- "FRCC"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 600, by = 200)) +
  coord_cartesian(ylim = c(0, 650)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed",
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 600, by = 200)) +
  coord_cartesian(ylim = c(0, 600)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 600, by = 200)) +
  coord_cartesian(ylim = c(0, 600)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 38, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 40, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 42)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000)
    ylim <- if (perc == 25) c(0.1, 9000) else c(0.1, 5000)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 5000, 2000), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)




#### ISONE ####
gea_example <- "ISONE"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 200)) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed",
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 200, by = 100)) +
  coord_cartesian(ylim = c(0, 200)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 200, by = 100)) +
  coord_cartesian(ylim = c(0, 200)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 18, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 20, by = 5), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 21)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100)
    ylim <- if (perc == 25) c(0.1, 500) else c(0.1, 150)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 300, 100), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### MISO Central ####
gea_example <- "MISO_Central"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 800, by = 200)) +
  coord_cartesian(ylim = c(0, 850)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed",
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 300, by = 100)) +
  coord_cartesian(ylim = c(0, 300)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.5, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 200, by = 100)) +
  coord_cartesian(ylim = c(0, 200)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 28, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 30, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 32)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000)
    ylim <- c(0.1, 4000)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 2000, 2000), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)




#### MISO North ####
gea_example <- "MISO_North"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 800, by = 200)) +
  coord_cartesian(ylim = c(0, 800)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed",
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 500, by = 200)) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.5, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 500, by = 200)) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 28, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 30, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 32)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000)
    ylim <- if (perc == 25) c(0.1, 4000) else c(0.1, 400)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 2000, 200), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)




#### MISO South ####
gea_example <- "MISO_South"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 600, by = 200)) +
  coord_cartesian(ylim = c(0, 650)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed",
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 200)) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.5, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 200)) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 28, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 30, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 32)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000)
    ylim <- if (perc == 25) c(0.1, 5000) else c(0.1, 900)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 3000, 600), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### NYISO ####
gea_example <- "NYISO"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 600, by = 200)) +
  coord_cartesian(ylim = c(0, 650)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed",
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 200)) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.5, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 200)) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_color_manual(values = ls_colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 9, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 10, by = 5), 
                         labels = number_format(suffix = " %")) +
      # scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 11)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100)
    ylim <- if (perc == 25) c(0.1, 300) else c(0.1, 80)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 100, 40), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### PJM WEST ####
gea_example <- "PJM_West"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 200)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 450)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.5, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 450)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 420)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 58, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 60, by = 20), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 65)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000)
    ylim <- if (perc == 25) c(0.1, 800) else c(0.1, 300)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 500, 200), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### SERTP ####
gea_example <- "SERTP"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 600, by = 200)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 600)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 300, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 350)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 300, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 350)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 38, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 40, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 42)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000, 10000)
    ylim <- if (perc == 25) c(0.1, 4000) else c(0.1, 800)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 2000, 500), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %", "10000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### SPP North ####
gea_example <- "SPP_North"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 800, by = 200)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 850)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 200, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 100, by = 50)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 150)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 58, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 60, by = 20), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 65)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000, 10000)
    ylim <- if (perc == 25) c(0.1, 4000) else c(0.1, 800)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 2000, 500), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %", "10000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### SPP South ####
gea_example <- "SPP_South"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 800, by = 200)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 800)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 300, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 350)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 200, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 58, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 60, by = 20), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 65)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000, 10000)
    ylim <- if (perc == 25) c(0.1, 4000) else c(0.1, 800)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 2000, 500), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %", "10000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### WestConnect North ####
gea_example <- "WestConnect_North"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 800, by = 200)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 850)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 200, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 100, by = 50)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 150)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 38, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 40, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 42)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000, 10000, 100000)
    ylim <- if (perc == 25) c(0.1, 500000) else c(0.1, 5000)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 250000, 2000), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %", "10000 %", "100000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### WestConnect South ####
gea_example <- "WestConnect_South"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 600, by = 200)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 600)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 300, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 350)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 200, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 58, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 60, by = 20), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 65)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000, 10000)
    ylim <- if (perc == 25) c(0.1, 8000) else c(0.1, 5000)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 5000, 3000), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %", "10000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)




#### NorthernGrid East ####
gea_example <- "NorthernGrid_East"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 800, by = 200)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 850)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 200, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 100, by = 50)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 150)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 38, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 40, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 42)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000, 10000)
    ylim <- if (perc == 25) c(0.1, 8000) else c(0.1, 800)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 4000, 500), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %", "10000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)




#### NorthernGrid South ####
gea_example <- "NorthernGrid_South"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 800, by = 200)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 850)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 400)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 1.4, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 400, by = 100)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 400)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 58, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 60, by = 20), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 65)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(1, 10, 100, 1000, 10000, 10000, 100000)
    ylim <- if (perc == 25) c(1, 800000) else c(1, 8000)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 400000, 4000), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("1 %", "10 %", "100 %", "1000 %", "10000 %", "100000 %", "1000000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)





#### NorthernGrid West ####
gea_example <- "NorthernGrid_West"

# all scenarios at hourly resolution at 2050
avg_annot <- df_annual %>% 
  filter(gea == gea_example, 
         year %in% c(2025, 2050), 
         scenario %in% c("MidCase", "LowRECost_HighNGPrice")) %>% 
  mutate(x_pos = if_else(year == 2025, as.POSIXct("2026-01-15 00:00:00"), as.POSIXct("2051-01-15 00:00:00")))

p1 <- df_hourly %>% 
  filter(year == 2025, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2025, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2025, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 1.2, color = "Annual avg.", label = er), 
            size = 5, 
            show.legend = FALSE) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 200, by = 50)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 200)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2025") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p2 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "MidCase") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "MidCase"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "MidCase"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 100, by = 50)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "MidCase - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))

p3 <- df_hourly %>% 
  filter(year == 2050, 
         gea == gea_example, 
         scenario == "LowRECost_HighNGPrice") %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = er), linewidth = 0.1, alpha = 0.6) +
  geom_hline(data = df_annual %>% 
               filter(year == 2050, 
                      gea == gea_example,
                      scenario == "LowRECost_HighNGPrice"), 
             aes(yintercept = er, color = "Annual avg."), 
             lty = "dashed", 
             linewidth = 1.2) +
  geom_text(data = avg_annot %>% 
              filter(year == 2050, 
                     scenario == "LowRECost_HighNGPrice"), 
            aes(x = x_pos, y = er * 2, color = "Annual avg.", label = er), 
            size = 5) +
  scale_x_datetime(expand = c(0.05, 0), labels = date_format("%b"), date_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0, 100, by = 50)) +
  scale_color_manual(values = ls_colors) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL,
       y = NULL,
       color = NULL,
       subtitle = "LowRECost_HighNGPrice - 2050") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        axis.text = element_text(size = 12), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "mm"))

final_plot <- ggarrange(p1, p2, p3,
                        nrow = 3, 
                        labels = c("a)", "b)", "c)"), 
                        common.legend = T, 
                        legend = "bottom")

y_axis <- expression(Emission ~ rates ~ (gCO[2] ~ e/kWh))

annotate_figure(final_plot, 
                left = text_grob(y_axis, rot = 90, vjust = 1, size = 12),
                top = text_grob(str_glue("Hourly AER at {gea_example}"), size = 14))

ggsave(filename = str_glue("{gea_example}_{emissions}_hourly.png"), path = figs_path, units = "in", height = 8, width = 10, dpi = 300)

# summary
sce_example <- c("MidCase", "LowRECost_HighNGPrice")

# operational
# subfigs_path <- paste0(figs_path, str_glue("operational/{g}/"))
z_index <- 1
plot_list <- list()
for (perc in c(0, 100)){
  
  
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
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
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
      filter(year %in% c(2025, 2050))
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey20", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = 38, label = label), parse = T, size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = seq(0, 40, by = 10), 
                         labels = number_format(suffix = " %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = c(0, 42)) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors) 
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of operational carbon emissions accounting", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_operational.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)


# avoided
z_index <- 1
plot_list <- list()
for (perc in c(25, 100)){
  
  
  for (z in sce_example){
    
    er_hourly <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/hourly.rds"))) %>% 
      select(contains(z))
    
    colname <- colnames(er_hourly)
    
    er_annual <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/annual.rds"))) %>% 
      select(all_of(colname))
    
    er_tod <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/tod.rds"))) %>% 
      select(all_of(colname))
    
    er_month_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/month_hour.rds"))) %>% 
      select(all_of(colname))
    
    er_season <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season.rds"))) %>% 
      select(all_of(colname))
    
    er_season_hour <- read_rds(str_glue(paste0(readfile_path, "/results/{emissions}/avoided/{gea_example}/{perc}/season_hour.rds"))) %>% 
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
      filter(year %in% c(2025, 2050))
    
    breaks <- c(0.1, 1, 10, 100, 1000, 10000)
    ylim <- if (perc == 25) c(0.1, 8000) else c(0.1, 800)
    
    er_hourly_med <- as.data.frame(t(apply(er_hourly, 2, median))) %>% 
      pivot_longer(everything(), names_to = "year", values_to = "er") %>% 
      separate(year, into = c("scenario", "year"), sep = "-") %>% 
      filter(year %in% c(2025, 2050)) %>% 
      mutate(pos = ifelse(perc == 25, 4000, 500), 
             er = ifelse(er < 0.1, "'< 0.1'", round(er, 1)), 
             label = paste0(er, "*~tCO[2]*e"))
    
    error_med <- df_error %>% 
      group_by(year, type, scenario) %>% 
      summarise(med = median(error)) %>% 
      ungroup() %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) 
    
    p <- df_error %>% 
      mutate(scenario = as.factor(scenario), 
             type = factor(type, levels = c("Annual avg.", "Season avg.", "Time-of-day avg.", "Season-hour avg.", "Month-hour avg."))) %>% 
      ggplot(aes(x = year, y = error)) +
      geom_lv(alpha = 0.4, k = 4, outlier.shape = NA, position = position_dodge(width = 0.8), color = "grey60", linewidth = 0.2, aes(fill = type)) +
      geom_text(data = er_hourly_med, aes(x = year, y = pos, label = label), parse = T, size = 5) +
      geom_text(data = error_med, aes(x = year, y = med, label = paste0(round(med, digits = 0), "%"), group = type), position = position_dodge(width = 0.8), size = 5) +
      scale_y_continuous(expand = c(0, 0), 
                         trans = "log10", 
                         breaks = breaks, 
                         labels = c("0.1 %", "1 %", "10 %", "100 %", "1000 %", "10000 %")) +
      scale_x_discrete(expand = c(-0.1, 0)) +
      coord_cartesian(ylim = ylim) +
      scale_fill_manual(values = ls_colors) +
      scale_color_manual(values = ls_colors)
    
    
    if (z_index %% 2 == 1){
      
      
      p <- p + 
        labs(x = NULL, 
             y = "Fractional error distribution", 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
      
    } else {
      
      p <- p + 
        labs(x = NULL, 
             y = NULL, 
             color = NULL, 
             fill = NULL, 
             subtitle = str_glue("{z}\n{perc}% PV offset")) +
        theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.2),
              legend.direction = "horizontal",
              legend.position = "bottom",
              axis.text.y = element_blank(), 
              plot.margin = margin(t = 1, r = 1, b = 15, l = 1, unit = "mm"))
      
    }
    
    plot_list[[z_index]] <- p
    
    z_index <- z_index + 1
    
  }
  
  
}

ggarrange(plotlist = plot_list, 
          ncol = 2, nrow = 2,
          align = "hv",
          common.legend = TRUE,
          labels = c("a)", "b)", "c)", "d)"),
          legend = "bottom") +
  plot_annotation(title = "Error distribution of quantifying avoided carbon emissions\nfrom exported utilities", 
                  subtitle = str_glue("{gea_example}"))

ggsave(filename = str_glue("{gea_example}_{emissions}_avoided.png"), path = figs_path, units = "in", height = 8, width = 14, dpi = 300)
