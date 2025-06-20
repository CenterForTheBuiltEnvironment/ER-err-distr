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




#### READ DATA ####
data_path <- "../../Cambium/"
output_path <- "../readfiles/"
col <- c("generation", 
         "variable_generation", 
         "battery_MWh",
         "biomass_MWh",
         "canada_MWh",
         "coal_MWh",
         "csp_MWh",
         "distpv_MWh",
         "gas-cc_MWh",
         "gas-ct_MWh",
         "geothermal_MWh",
         "hydro_MWh",
         "nuclear_MWh",
         "o-g-s_MWh",
         "phs_MWh",
         "upv_MWh",
         "wind-ons_MWh",
         "wind-ofs_MWh")

df_gen <- list.files(path = data_path, pattern = "annual", full.names = TRUE) %>% 
  map_dfr(~ read_csv(.x, skip = 5)) %>% 
  select(scenario, 
         gea, 
         year = t, 
         col) 





#### OUTPUT ####
write_rds(df_gen, paste0(output_path, "df_gen.rds"), compress = "gz")
