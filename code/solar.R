# AER error distribution analysis
# written by Aoyu



#### SETUP ####

# use pacman
require(pacman)

# load libraries
pacman::p_load(tidyverse, lubridate, scales, slider, cowplot, patchwork, RColorBrewer, broom, ggpmisc, ggpubr, eplusr)

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



#### FUNCTIONS ####
radians <- function(degrees) {
  return(degrees * (pi / 180))
}

degrees <- function(radians) {
  return(radians * (180 / pi))
}

#### READ ####
# epw file
readfile_path <- "../readfiles/weather/"
output_path <- paste0(readfile_path, "./solar/")

# solar_map <- read_csv(paste0(readfile_path, "solar_map.csv"))

epw <- read_epw(paste0(readfile_path, "Oakland_USA.epw"))
df_tmy <- epw$data() %>% 
  mutate(year = 2016,
         datetime = ymd_h(paste(paste(year, month, day, sep = "-"), hour, sep = " "))) %>% 
  select(c(datetime, dry_bulb_temperature, global_horizontal_radiation, diffuse_horizontal_radiation))





#### CALCULATE ####
# Time definition
n <- seq(1, 365, by = 1)
st <- seq(1, 24, by = 1)
L = nrow(df_tmy)

# Read solar radiation data
Rt_h <- df_tmy$global_horizontal_radiation
Rs_h <- df_tmy$diffuse_horizontal_radiation
Rd_h <- Rt_h - Rs_h

# Read ambient temperature data
T <- df_tmy$dry_bulb_temperature

# Photovoltaic cell parameters
Tref <- 25  # Standard temperature for polycrystalline silicon cells in °C
eta_ref <- 0.20  # Standard efficiency of polycrystalline silicon cells
c1 <- -3.75
c2 <- 1.14
c3 <- 0.0175
rou <- 0.0045  # Temperature degradation coefficient
sigma <- 0.1

gamma <- 0 # Azimuth angles in degrees: North 360, East 90, South 180, West 270

beta <- 15
fai <- 30 + 16 / 60

# Initialize variables
cos_theta <- array(dim = L)
sin_a <- array(dim = L)
a <- array(dim = L)
A <- array(dim = L)

# Calculate solar altitude angle
for (i in 1:L) {
  temp <- ceiling(i / 24)
  n_n <- n[temp]
  st_t <- (i %% 24)
  
  omega <- radians(15 * (st_t - 12))
  delta <- radians(23.45 * sin(2 * pi * ((284 + n_n) / 365)))
  sin_a[i] <- sin(radians(fai)) * sin(delta) + cos(radians(fai)) * cos(delta) * cos(omega)
  a[i] <- degrees(asin(sin_a[i]))
  
  A[i] <- degrees(atan(sin(omega) / (cos(radians(fai)) * tan(delta) - sin(radians(fai)) * cos(omega))))
  cos_theta[i] <- (sin(radians(fai)) * cos(radians(beta)) - cos(radians(fai)) * cos(radians(gamma)) * sin(radians(beta))) * sin(delta) + 
    (cos(radians(fai)) * cos(radians(beta)) + sin(radians(fai)) * cos(radians(gamma)) * sin(radians(beta))) * cos(delta) * cos(omega) + 
    sin(radians(gamma)) * sin(radians(beta)) * cos(delta) * sin(omega)
}

# Convert solar radiation on a plane to tilted surface solar radiation
Rd <- array(dim = L)
for (i in 1:L) {
  if (cos_theta[i] / sin_a[i] <= 4 && cos_theta[i] / sin_a[i] >= 0) {
    Rd[i] <- Rd_h[i] * cos_theta[i] / sin_a[i]
  } else {
    Rd[i] <- 0
  }
}

Rr <- 0.15 * Rt_h * (1 - cos(radians(beta))) / 2  # Reflected solar radiation on tilted surface
Rs <- Rs_h * (1 + cos(radians(beta))) / 2  # Scattered solar radiation on tilted surface
Rt <- Rd + Rs + Rr

# Calculate PV power output
Ppv <- array(dim = L)
for (i in 1:L) {
  if (Rt[i] > 0) {
    Ppv[i] <- Rt[i] * eta_ref * (1 - rou * (c1 + c2 * T[i] + c3 * Rt[i] - Tref) + sigma * log10(Rt[i]))
  } else {
    Ppv[i] <- 0
  }
}

Ppv_cal <- 9 * 0.496 * Ppv
df_solar <- data.frame(datetime = df_tmy$datetime, PV = Ppv_cal)




#### PLOT ####
df_solar %>% 
  ggplot(aes(x = datetime, y = PV)) +
  geom_line(color = "#fe9929") +
  scale_x_datetime(date_breaks = "2 months",
                   date_labels = "%b")  +
  scale_y_continuous(expand = c(0, 0),
                     breaks = breaks_pretty(n = 4),
                     labels = number_format(suffix = " W/m^2")) +
  labs(x = NULL,
       y = "PV panel power",
       title = "Calculated solar PV power generation") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"), 
        axis.text = element_text(size = 12))




#### EXPORT ####
