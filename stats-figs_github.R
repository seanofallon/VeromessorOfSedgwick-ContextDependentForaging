###### Data analysis and visualization for "Foraging response to environmental cues
###### differs across contexts in a harvester ant" by Sean O'Fallon and Noa
###### Pinter-Wollman. Inputs are included in GitHub repo and should be saved
###### to your working directory for the script to run.

#### loading ####
# clear env, load packages
rm(list = ls())

# load packages
# basics
library(dplyr)
library(tidyr)
library(stringr)

# space/time
library(lubridate)
library(hms)
library(maptiles)
library(sf)

# ggfamily / viz
library(ggplot2)
library(ggeffects)
library(GGally)
library(ggspatial)
library(patchwork)
library(ggmap)
library(gt)
library(broom.mixed)
library(flextable)
library(officer)
library(purrr)

# stats
library(easystats)
library(lme4)
library(corrplot)
library(glmmTMB)
library(lmerTest)
library(DHARMa)
library(VGAM)
library(car)
library(emmeans)
library(performance)

#### inputs ####
# inputs and how to use them
all_obs23 <- read.csv('all_obs23.csv') # all valid 2023 observations along w/ env and behavioral measurements
# ^ 'tagged' w/ columns 'repeated', 'foraging', and 'steady' for overlapping observation
# qualities: filter as needed
col_loc23 <- read.csv('col_loc.csv') # lon/lat for each colony as well as basic info (field/road, rep/not)
# ^ use for map making

all_obs24 <- read.csv('all_obs24.csv') # all valid 2024 observations along w/ env and behavioral measurements
# ^ 'tagged' w/ columns 'repeated', 'foraging', and 'steady' for overlapping observation
# qualities: filter as needed
col_loc24 <- read.csv('col_loc24.csv') # lon/lat for each colony as well as basic info (field/road, rep/not)
# ^ use for map making

airstrip <- read.csv('airstrip.csv') # weather station data, used for Fig S1

### fix formatting
all_obs23$Colony <- as.factor(all_obs23$Colony)
all_obs24$Colony <- as.factor(all_obs24$Colony)

## correct order morning then evening
# 2023
all_obs23$MorningOrEvening <- factor(
  all_obs23$MorningOrEvening,
  levels = c("morning", "evening")
)
# 2024
all_obs24$MorningOrEvening <- factor(
  all_obs24$MorningOrEvening,
  levels = c("morning", "evening")
)

## get Time in readable format, TimeOfDay to look at general trends
# 2023
all_obs23 <- all_obs23 %>%
  # convert character to POSIXct
  mutate(
    Time = ymd_hms(Time, tz = "America/Los_Angeles"),   # choose local timezone
    TimeOfDay = hour(Time) + minute(Time)/60            # numeric hours since midnight
  )# get Time in readable format, TimeOfDay to look at general trends
# 2024
all_obs24 <- all_obs24 %>%
  # convert character to POSIXct
  mutate(
    Time = ymd_hms(Time, tz = "America/Los_Angeles"),   # choose local timezone
    TimeOfDay = hour(Time) + minute(Time)/60            # numeric hours since midnight
  )

# obs from exp colonies 2024 (a few colonies were selected but didn't get enough obs)
all_obs24 <- all_obs24[all_obs24$exp == TRUE,]

### steady state obs ###
# repeated obs w/ steady state in 2023
ssrep_obs23 <- all_obs23[all_obs23$repeated == TRUE &
                           all_obs23$steady == TRUE,]
# obs w/ steady state foraging in summer 2024
steady_sum24 <- all_obs24[all_obs24$steady == TRUE &
                            all_obs24$Season == 'summer',]

### foraging obs ###
# repeated obs w/ foraging in 2023
forrep_obs23 <- all_obs23[all_obs23$repeated == TRUE &
                            all_obs23$foraging == TRUE,]
# single potentially influential LuxMean observation (92050)
forrep_obs23 <- forrep_obs23[forrep_obs23$LuxMean != 92050, ]


# all obs w/ trail foraging summer 2024
for_sum24 <- all_obs24[all_obs24$trail.foraging == TRUE &
                         all_obs24$Season == 'summer',]



#### stats ####
### foraging flow in all foraging obs
#### 2023 ####
#### air ####
# Scale predictors
forrep_air_data23 <- forrep_obs23 %>%
  mutate(across(
    c(AirTemp, RH, IrradianceAvg),
    ~ as.numeric(scale(.x))
  ))

# Fit the model
ffair23_mod1 <- lmer(
  Foraging.Flow ~ IrradianceAvg * MorningOrEvening +
    AirTemp * MorningOrEvening +
    RH * MorningOrEvening +
    (1 | Colony),
  data = forrep_air_data23
)

# Check assumptions and diagnostics
check_model(ffair23_mod1)
# Check model summary
summary(ffair23_mod1)
# check parameters
model_parameters(ffair23_mod1)

## use emtrends to find whether context-specific slopes
## of sig interactions are significant
# Irradiance
ffair23_trendIrr <- summary(
  emtrends(ffair23_mod1, 
           ~ MorningOrEvening, 
           var = "IrradianceAvg"),
  infer = TRUE
)

#### soil ####
# Scale predictors
forrep_soil_data23 <- forrep_obs23 %>%
  mutate(across(
    c(SoilTemp, SoilMoisture, LuxMean),
    ~ as.numeric(scale(.x))
  ))

# Fit the model
ffsoil23_mod1 <- lmer(
  Foraging.Flow ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    (1 | Colony),
  data = forrep_soil_data23
)

# Check assumptions and diagnostics
check_model(ffsoil23_mod1)
# Check model summary
summary(ffsoil23_mod1)
# check parameters
model_parameters(ffsoil23_mod1)

## use emtrends to find whether context-specific slopes
## of sig interactions are significant
# Lux
ffsoil23_trendLux <- summary(
  emtrends(ffsoil23_mod1, 
           ~ MorningOrEvening, 
           var = "LuxMean"),
  infer = TRUE
)

# SoilTemp
ffsoil23_trendSoilTemp <- summary(
  emtrends(ffsoil23_mod1, 
           ~ MorningOrEvening, 
           var = "SoilTemp"),
  infer = TRUE
)

#### 2024 ####
##### air ####
# Scale predictors
forrep_air_data24 <- for_sum24 %>%
  mutate(across(
    c(AirTemp, RH, IrradianceAvg),
    ~ as.numeric(scale(.x))
  ))

# Fit the model
ffair24_mod1 <- lmer(
  Foraging.Flow ~ IrradianceAvg * MorningOrEvening +
    AirTemp * MorningOrEvening +
    RH * MorningOrEvening +
    (1 | Colony),
  data = forrep_air_data24
)

# Check assumptions and diagnostics
check_model(ffair24_mod1)
# Check model summary
summary(ffair24_mod1)
# check parameters
model_parameters(ffair24_mod1)

## use emtrends to find whether context-specific slopes
## of sig interactions are significant
# Irradiance
ffair24_trendIrr <- summary(
  emtrends(ffair24_mod1, 
           ~ MorningOrEvening, 
           var = "IrradianceAvg"),
  infer = TRUE
)

#### soil (ignore) ####
# Scale predictors
forrep_soil_data24 <- for_sum24 %>%
  mutate(across(
    c(SoilTemp, SoilMoisture, LuxMean),
    ~ as.numeric(scale(.x))
  ))

# Fit the model
ffsoil24_mod1 <- lmer(
  Foraging.Flow ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    (1 | Colony),
  data = forrep_soil_data24
)

# Check assumptions and diagnostics
check_model(ffsoil24_mod1)
# Check model summary
summary(ffsoil24_mod1)
# check parameters
model_parameters(ffsoil24_mod1)




#### all traffic in steady state obs ####
#### 2023 ####
#### air ####
# scale explanatory variables
ssair23_data <- ssrep_obs23 %>%
  mutate(across(
    c(AirTemp, IrradianceAvg, RH),
    ~ as.numeric(scale(.x))
  ))

ssair23_mod1 <- glmmTMB(
  All.Traffic ~ IrradianceAvg * MorningOrEvening +
    AirTemp * MorningOrEvening +
    RH * MorningOrEvening +
    (1 | Colony),
  family = nbinom2,
  data = ssair23_data
)

summary(ssair23_mod1)
check_model(ssair23_mod1)
model_parameters(ssair23_mod1)

#### soil ####
# Scale predictors
sssoil23_data <- ssrep_obs23 %>%
  mutate(across(
    c(SoilTemp, SoilMoisture, LuxMean),
    ~ as.numeric(scale(.x))
  ))

# Fit the model
sssoil23_mod1 <- glmmTMB(
  All.Traffic ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    (1 | Colony),
  family = nbinom2,
  data = sssoil23_data
)

# Check assumptions and diagnostics
check_model(sssoil23_mod1)
# Check model summary
summary(sssoil23_mod1)
# check parameters
model_parameters(sssoil23_mod1)

## use emtrends to find whether each slope is significant
# Lux
sssoil23_trendLux <- summary(
  emtrends(sssoil23_mod1, 
           ~ MorningOrEvening, 
           var = "LuxMean"),
  infer = TRUE
)

# SoilTemp
sssoil23_trendSoilTemp <- summary(
  emtrends(sssoil23_mod1, 
           ~ MorningOrEvening, 
           var = "SoilTemp"),
  infer = TRUE
)

#### 2024 ####
#### air ####
# scale explanatory variables
ssair24_data <- steady_sum24 %>%
  mutate(across(
    c(AirTemp, IrradianceAvg, RH),
    ~ as.numeric(scale(.x))
  ))

ssair24_mod1 <- glmmTMB(
  All.Traffic ~ IrradianceAvg * MorningOrEvening +
    AirTemp * MorningOrEvening +
    RH * MorningOrEvening +
    (1 | Colony),
  family = nbinom2,
  data = ssair24_data
)

summary(ssair24_mod1)
check_model(ssair24_mod1)
model_parameters(ssair24_mod1)

#### soil ####
# Scale predictors
sssoil24_data <- steady_sum24 %>%
  mutate(across(
    c(SoilTemp, SoilMoisture, LuxMean),
    ~ as.numeric(scale(.x))
  ))

# Fit the model
sssoil24_mod1 <- glmmTMB(
  All.Traffic ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    (1 | Colony),
  family = nbinom2,
  data = sssoil24_data
)

# Check assumptions and diagnostics
check_model(sssoil24_mod1)
# Check model summary
summary(sssoil24_mod1)
# check parameters
model_parameters(sssoil24_mod1)

# SoilTemp
sssoil23_trendSoilTemp <- summary(
  emtrends(sssoil23_mod1, 
           ~ MorningOrEvening, 
           var = "SoilTemp"),
  infer = TRUE
)

#### 2024 treatment effects ####
#### flow air (no effect, worse AIC, ignore) ####
# flow air model w/ TREATMENT
ffair24_modT <- lmer(
  Foraging.Flow ~ IrradianceAvg * MorningOrEvening +
    AirTemp * MorningOrEvening +
    RH * MorningOrEvening + 
    Treatment * MorningOrEvening +
    (1 | Colony),
  data = forrep_air_data24
)

# Check assumptions and diagnostics
check_model(ffair24_modT)
# Check model summary
summary(ffair24_modT)
# check parameters
model_parameters(ffair24_modT)

AIC(ffair24_modT, ffair24_mod1)

#### flow soil (treatment effect, =AIC, use) ####
# flow soil model w/ TREATMENT
ffsoil24_modT <- lmer(
  Foraging.Flow ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    Treatment * MorningOrEvening +
    (1 | Colony),
  data = forrep_soil_data24
)

# Check assumptions and diagnostics
check_model(ffsoil24_modT)
# Check model summary
summary(ffsoil24_modT)
# check parameters
model_parameters(ffsoil24_modT)

# confirm that treatment improves model - so keep
AIC(ffsoil24_modT, ffsoil24_mod1)

## use emtrends to find whether context-specific slopes
## of sig interactions are significant
# Lux
ffsoil24_trendLux <- summary(
  emtrends(ffsoil24_modT, 
           ~ MorningOrEvening, 
           var = "LuxMean"),
  infer = TRUE
)

# SoilTemp
ffsoil24_trendSoilTemp <- summary(
  emtrends(ffsoil24_modT, 
           ~ MorningOrEvening, 
           var = "SoilTemp"),
  infer = TRUE
)

## emmeans to find what drives relationship
# Treatment
emm_t <- emmeans(ffsoil24_modT, ~ Treatment)

emm_t

pairs(emm_t, adjust = "tukey")

# Treatment x Context
emm_txC <- emmeans(ffsoil24_modT, ~ Treatment | MorningOrEvening)

emm_txC

pairs(emm_txC, adjust = "tukey")


# Context x Treatment
emm_cXt <- emmeans(ffsoil24_modT, ~ MorningOrEvening | Treatment)

pairs(emm_cXt)


#### traffic air (no effect, worse AIC, ignore) ####
# traffic air model w/ TREATMENT
ssair24_modT <- glmmTMB(
  All.Traffic ~ IrradianceAvg * MorningOrEvening +
    AirTemp * MorningOrEvening +
    RH * MorningOrEvening +
    Treatment * MorningOrEvening +
    (1 | Colony),
  family = nbinom2,
  data = ssair24_data
)

summary(ssair24_modT)
check_model(ssair24_modT)
model_parameters(ssair24_modT)

AIC(ssair24_modT, ssair24_mod1)

#### traffic soil (no effect, worse AIC, ignore) ####
# traffic soil model w/ TREATMENT
sssoil24_modT <- glmmTMB(
  All.Traffic ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    Treatment * MorningOrEvening +
    (1 | Colony),
  family = nbinom2,
  data = sssoil24_data
)


summary(sssoil24_modT)
check_model(sssoil24_modT)
model_parameters(sssoil24_modT)

AIC(sssoil24_modT, sssoil24_mod1)

#### batch colony effect ####
model_list <- list(
  ffair23_mod1 = ffair23_mod1,
  ssair23_mod1 = ssair23_mod1,
  ffsoil23_mod1 = ffsoil23_mod1,
  sssoil23_mod1 = sssoil23_mod1,
  ffair24_mod1 = ffair24_mod1,
  ssair24_mod1 = ssair24_mod1,
  ffsoil24_modT = ffsoil24_modT,
  sssoil24_mod1 = sssoil24_mod1
)

r2_summary <- lapply(model_list, function(mod) {
  r2 <- r2_nakagawa(mod)
  data.frame(
    marginal = r2$R2_marginal,
    conditional = r2$R2_conditional,
    random_effect = r2$R2_conditional - r2$R2_marginal
  )
})

do.call(rbind, r2_summary)

#### maps ####
## w/ sf and maptiles
# convert lon/lat to sf format
col_loc23_sf <- col_loc23 %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# look at bounding box
st_bbox(col_loc23_sf)

# Shrink or expand as needed by modifying the bounds manually
bbox23 <- st_bbox(col_loc23_sf)
bbox23["xmin"] <- bbox23["xmin"] - 0.0001  # West
bbox23["xmax"] <- bbox23["xmax"] + 0.0001  # East
bbox23["ymin"] <- bbox23["ymin"] - 0.0001  # South
bbox23["ymax"] <- bbox23["ymax"] + 0.0001  # North

# use bbox to get_tiles()
bg_map23 <- get_tiles(x = bbox23, provider = "Esri.WorldImagery", zoom = 18)

# STEP 3: Plot
ggplot() +
  layer_spatial(bg_map23) +
  geom_sf(data = col_loc23_sf, 
          aes(color = repOrNot, shape = repOrNot),
          size = 4) +  
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 1.2,
    line_width = 1.2,
    bar_cols = c("white", "black")
  ) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  coord_sf(xlim = c(bbox23["xmin"], bbox23["xmax"]),
           ylim = c(bbox23["ymin"], bbox23["ymax"]),
           expand = FALSE) +  
  theme_minimal() +
  labs(title = "2023 Colony Locations") +
  scale_color_manual(values = c("FALSE" = "#000000", "TRUE" = "#00BFFF"),
                     labels = c("Surveyed", "Assayed"),
                     name = NULL) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("Surveyed", "Assayed"),
                     name = NULL) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_blank(),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 13),
    plot.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black"),  # white box with border
    legend.position = c(0.85, 0.9)
  )



### 2024
## w/ sf and maptiles
# STEP 1: Convert to sf
col_loc_sf24 <- col_loc24 %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)

# look at bounding box
st_bbox(col_loc_sf24)

# Shrink or expand as needed by modifying the bounds manually
bbox24 <- st_bbox(col_loc_sf24)
bbox24["xmin"] <- bbox24["xmin"] - 0.0001  # West
bbox24["xmax"] <- bbox24["xmax"] + 0.0001  # East
bbox24["ymin"] <- bbox24["ymin"] - 0.0001  # South
bbox24["ymax"] <- bbox24["ymax"] + 0.0001  # North

# Now use this in get_tiles()
bg_map24 <- get_tiles(x = bbox24, provider = "Esri.WorldImagery", zoom = 18)

# STEP 3: Plot
ggplot() +
  layer_spatial(bg_map24) +
  geom_sf(data = col_loc_sf24, 
          aes(color = Treatment, shape = Treatment),
          size = 4) + 
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 1.2,
    line_width = 1.2,
    bar_cols = c("white", "black")
  ) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_fancy_orienteering()) +
  coord_sf(xlim = c(bbox24["xmin"], bbox24["xmax"]),
           ylim = c(bbox24["ymin"], bbox24["ymax"]),
           expand = FALSE) +  
  theme_minimal() +
  labs(title = "2024 Colony Locations") +
  scale_color_manual(values = c("shade" = "magenta", "water" = "green", "control" = "yellow"),
                     labels = c("control" = "Control","water" = "Water","shade" = "Shade"),
                     name = NULL) +
  scale_shape_manual(values = c("shade" = 16, "water" = 17, "control" = 18),
                     labels = c("control" = "Control","water" = "Water","shade" = "Shade"),
                     name = NULL) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_blank(),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 13),
    plot.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black"),  # white box with border
    legend.position = c(0.9, 0.9)  # inside map, adjust x/y as needed
  )


#### figures ####
# Assign colors explicitly
colors <- c("morning" = "#ffd859", "evening" = "#7570b3")
#### fig 1 - conceptual ####
# FORAGING FLOW VS. ALL TRAFFIC
all_obs23 <- all_obs23 %>%
  mutate(steady_dir = case_when(
    steady ~ "steady",
    Foraging.Flow >= 1.089 ~ "above",
    Foraging.Flow <= -1.089 ~ "below",
    Foraging.Flow > -1.089 & Foraging.Flow < 1.089 ~ "steady"
  ))

ggplot(data = all_obs23, aes(x = All.Traffic, y = Foraging.Flow, color = steady_dir)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "below" = "#dc5809",
      "steady" = "black",
      "above" = "#074761"
    )) +
  labs(
    x = "Total Traffic",
    y = "Foraging Flow"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 14),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "none",
    legend.justification = c("right", "bottom")
  )


#### fig 2 - morning vs. evening ####

# traffic 2023
mVe_traffic23 <- ggplot(ssrep_obs23, aes(x = MorningOrEvening, y = All.Traffic, fill = MorningOrEvening)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.9) +
  scale_y_continuous(limits = c(0, 260)) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Total Traffic") +
  annotate("text", x = -Inf, y = Inf, label = "C",
           hjust = -0.1, vjust = 1.5, size = 10, fontface = "bold") +
  annotate("text", x = 1.5, y = 250, label = "*", size = 12) +
  theme_minimal(base_size = 20) +
  theme(
    axis.title = element_text(face = "bold", size = 22),
    axis.text = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# flow 2023
mVe_flow23 <- ggplot(forrep_obs23, aes(x = MorningOrEvening, y = Foraging.Flow, fill = MorningOrEvening)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.9) +
  scale_y_continuous(limits = c(-6, 5)) +
  scale_fill_manual(values = colors) +
  labs(title = "Observational Assays (2023)", x = NULL, y = "Foraging Flow") +
  annotate("text", x = -Inf, y = Inf, label = "A",
           hjust = -0.1, vjust = 1.5, size = 10, fontface = "bold") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 22),
    axis.text = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# traffic 2024
mVe_traffic24 <- ggplot(steady_sum24, aes(x = MorningOrEvening, y = All.Traffic, fill = MorningOrEvening)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.9) +
  scale_y_continuous(limits = c(0, 260)) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = NULL) +
  annotate("text", x = -Inf, y = Inf, label = "D",
           hjust = -0.1, vjust = 1.5, size = 10, fontface = "bold") +
  annotate("text", x = 1.5, y = 250, label = "*", size = 12) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold", size = 22),
    axis.text = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# flow 2024
mVe_flow24 <- ggplot(for_sum24, aes(x = MorningOrEvening, y = Foraging.Flow, fill = MorningOrEvening)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.9) +
  scale_y_continuous(limits = c(-6, 5)) +
  scale_fill_manual(values = colors) +
  labs(title = "Experimental Trials (2024)", x = NULL, y = NULL) +
  annotate("text", x = -Inf, y = Inf, label = "B",
           hjust = -0.1, vjust = 1.5, size = 10, fontface = "bold") +
  annotate("text", x = 1.5, y = 4.5, label = "*", size = 12) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    axis.title.y = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
    axis.text = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# plot together
fig2 <- (mVe_flow23 | mVe_flow24) /
  (mVe_traffic23 | mVe_traffic24)

#### fig 3 - flow vs. treatment x context ####
# just look at flow - traffic not affected
fig3 <- ggplot(for_sum24, 
               aes(x = Treatment, 
                   y = Foraging.Flow, 
                   fill = MorningOrEvening)) +
  geom_boxplot(position = position_dodge(width = 0.75),
               outlier.shape = 21,
               outlier.size = 2,
               alpha = 0.9) +
  scale_x_discrete(labels = c("control" = "Control",
                              "shade"   = "Shade",
                              "water"   = "Water")) +
  scale_fill_manual(values = colors,
                    name = "Context") +
  labs(x = "Treatment",
       y = "Foraging Flow") +
  annotate("text", x = 2, y = 2.1, label = "*", size = 8) +
  theme_minimal(base_size = 18) +   
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(size = 16),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 16),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

#### fig 4 - flow vs. environment ####
# use forrep_obs23 for unscaled predictor values - model is qualitatively the same but for context due to non-0's
ffsoil23_modU <- lmer(
  Foraging.Flow ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    (1 | Colony),
  data = forrep_obs23
)

ff23_lux_preds <- ggpredict(ffsoil23_modU, terms = c("LuxMean", "MorningOrEvening"))

ff23_lux_preds$MorningOrEvening <- ff23_lux_preds$group

ff23_lux_plot <- ggplot(ff23_lux_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(data = forrep_obs23, aes(x = LuxMean, y = Foraging.Flow, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(-10, 81000), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-6, 5)) +
  labs(
    x = "Light Intensity (Lux)",
    y = "Foraging Flow (2023)"
  ) +
  annotate("text", x = 4000, y = Inf, label = "B",
           hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )

# FF23 - Soil Temp
ff23_soilTemp_preds <- ggpredict(ffsoil23_modU, terms = c("SoilTemp", "MorningOrEvening"))

ff23_soilTemp_preds$MorningOrEvening <- ff23_soilTemp_preds$group

ff23_soilTemp_plot <- ggplot(ff23_soilTemp_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2, linetype = "dashed") +
  geom_point(data = forrep_obs23, aes(x = SoilTemp, y = Foraging.Flow, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(23, 33), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-6, 5), expand = c(0, 0)) +
  labs(
    x = "Soil Temperature (°C)",
    y = "Foraging Flow (2023)"
  ) +
  annotate("text", x = 23.5, y = Inf, label = "A",
           hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )

## 2024 flow
ffsoil24_modU <- lmer(
  Foraging.Flow ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    Treatment * MorningOrEvening +
    (1 | Colony),
  data = for_sum24
)

ff24_lux_preds <- ggpredict(ffsoil24_modU, terms = c("LuxMean", "MorningOrEvening"))

ff24_lux_preds$MorningOrEvening <- ff24_lux_preds$group

ff24_lux_plot <- ggplot(ff24_lux_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(data = for_sum24, aes(x = LuxMean, y = Foraging.Flow, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(-10, 81000), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-6, 5)) +
  labs(
    x = "Light Intensity (Lux)",
    y = "Foraging Flow (2024)"
  ) +
  annotate("text", x = 4000, y = Inf, label = "E",
           hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )



# FF24 - Soil Temp
ff24_soilTemp_preds <- ggpredict(ffsoil24_modU, terms = c("SoilTemp", "MorningOrEvening"))

ff24_soilTemp_preds$MorningOrEvening <- ff24_soilTemp_preds$group

ff24_soilTemp_plot <- ggplot(ff24_soilTemp_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(aes(linetype = MorningOrEvening), linewidth = 1.2) +
  geom_point(data = for_sum24, aes(x = SoilTemp, y = Foraging.Flow, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_linetype_manual(values = c("morning" = "solid", "evening" = "dashed")) +
  scale_color_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(23, 33), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-6, 5), expand = c(0, 0)) +
  labs(
    y = "Foraging Flow (2024)",
    x = "Soil Temperature (°C)") +
  annotate("text", x = 23.5, y = Inf, label = "D",
           hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )

# FF23 - Irradiance
## 2023 air flow
ffair23_modU <- lmer(
  Foraging.Flow ~ IrradianceAvg * MorningOrEvening +
    AirTemp * MorningOrEvening +
    RH * MorningOrEvening +
    (1 | Colony),
  data = ssrep_obs23
)

ff23_irr_preds <- ggpredict(ffair23_modU, terms = c("IrradianceAvg", "MorningOrEvening"))

ff23_irr_preds$MorningOrEvening <- ff23_irr_preds$group

ff23_irradiance_plot <- ggplot(ff23_irr_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(data = for_sum24, aes(x = IrradianceAvg, y = Foraging.Flow, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_color_manual(name = "Context",
                     values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(name = "Context",
                    values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(0, 500), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-6, 5)) +
  labs(x = expression(Irradiance~(W/m^2))) +
  annotate("text", x = 25, y = Inf, label = "C",
           hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = c(0.85, 0.9),
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )

# FF24 - Irradiance
## 2024 air flow
ffair24_modU <- lmer(
  Foraging.Flow ~ IrradianceAvg * MorningOrEvening +
    AirTemp * MorningOrEvening +
    RH * MorningOrEvening +
    (1 | Colony),
  data = for_sum24
)

ff24_irr_preds <- ggpredict(ffair24_modU, terms = c("IrradianceAvg", "MorningOrEvening"))

ff24_irr_preds$MorningOrEvening <- ff24_irr_preds$group

ff24_irradiance_plot <- ggplot(ff24_irr_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(data = for_sum24, aes(x = IrradianceAvg, y = Foraging.Flow, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(0, 500), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-6, 5)) +
  labs(x = expression(Irradiance~(W/m^2))) +
  annotate("text", x = 25, y = Inf, label = "F",
           hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 24, margin = margin(t = -5)),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )


# plot together
fig4 <- (ff23_soilTemp_plot | ff23_lux_plot | ff23_irradiance_plot) /
  (ff24_soilTemp_plot | ff24_lux_plot | ff24_irradiance_plot)


#### fig 5 - AT vs. environment ####
sssoil23_modU <- glmmTMB(
  All.Traffic ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    (1 | Colony),
  family = nbinom2,
  data = ssrep_obs23
)

ss23_lux_preds <- ggpredict(sssoil23_modU, terms = c("LuxMean", "MorningOrEvening"))

ss23_lux_preds$MorningOrEvening <- ss23_lux_preds$group

ss23_lux_plot <- ggplot(ss23_lux_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2, linetype = "dashed") +
  geom_point(data = ssrep_obs23, aes(x = LuxMean, y = All.Traffic, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(0, 55000), expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0, 250)) +
  labs(
    x = "Light Intensity (Lux)",
    y = "Total Traffic (2023)"
  ) +
  annotate("text", x = 3000, y = Inf, label = "A",
         hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )

# FF23 - Soil Temp
ss23_soilTemp_preds <- ggpredict(sssoil23_modU, terms = c("SoilTemp", "MorningOrEvening"))

ss23_soilTemp_preds$MorningOrEvening <- ss23_soilTemp_preds$group

ss23_soilTemp_plot <- ggplot(ss23_soilTemp_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(aes(linetype = MorningOrEvening), linewidth = 1.2) +
  geom_point(data = ssrep_obs23, aes(x = SoilTemp, y = All.Traffic, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_linetype_manual(name = "Context",
                        values = c("morning" = "solid", "evening" = "dashed")) +
  scale_color_manual(name = "Context",
                     values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(name = "Context",
                    values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(22.5, 32.6), expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0, 250)) +
  labs(
    x = "Soil Temperature (°C)",
    y = "Total Traffic (2023)"
  ) +
  annotate("text", x = 23, y = Inf, label = "B",
           hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = c(0.25,0.85),
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )


## SS 24
sssoil24_modU <- glmmTMB(
  All.Traffic ~ LuxMean * MorningOrEvening +
    SoilTemp * MorningOrEvening +
    SoilMoisture * MorningOrEvening +
    (1 | Colony),
  family = nbinom2,
  data = steady_sum24
)

ss24_lux_preds <- ggpredict(sssoil24_modU, terms = c("LuxMean", "MorningOrEvening"))

ss24_lux_preds$MorningOrEvening <- ss24_lux_preds$group

ss24_lux_plot <- ggplot(ss24_lux_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2, linetype = "dashed") +
  geom_point(data = steady_sum24, aes(x = LuxMean, y = All.Traffic, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(0, 55000), expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0, 250)) +
  labs(
    x = "Light Intensity (Lux)",
    y = "Total Traffic (2024)"
  ) +
  annotate("text", x = 3000, y = Inf, label = "C",
           hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 24),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )


ss24_soilTemp_preds <- ggpredict(sssoil24_modU, terms = c("SoilTemp", "MorningOrEvening"))

ss24_soilTemp_preds$MorningOrEvening <- ss24_soilTemp_preds$group

ss24_soilTemp_plot <- ggplot(ss24_soilTemp_preds, aes(x = x, y = predicted, color = MorningOrEvening, fill = MorningOrEvening)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2, linetype = "dashed") +
  geom_point(data = steady_sum24, aes(x = SoilTemp, y = All.Traffic, color = MorningOrEvening), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_fill_manual(values = c("morning" = "#ffd859", "evening" = "#7570b3")) +
  scale_x_continuous(limits = c(22.5, 32.6), expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0, 250)) +
  labs(
    x = "Soil Temperature (°C)",
    y = "Total Traffic (2024)"
  ) +
  annotate("text", x = 23, y = Inf, label = "D",
           hjust = -0.2, vjust = 1.5, size = 10, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 24),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    legend.background = element_rect(fill = alpha("transparent", 0))  # make legend background transparent
  )

fig5 <- (ss23_lux_plot | ss23_soilTemp_plot) /
  (ss24_lux_plot | ss24_soilTemp_plot)

#### Fig S1 - weather over the day ####

airstrip_summer <- airstrip %>%
  
  # parse timestamps
  mutate(
    datetime = with_tz(
      ymd_hms(TIMESTAMP, tz = "UTC"),
      tzone = "America/Los_Angeles"
    )
  ) %>%
  
  # keep June + July
  filter(month(datetime) %in% c(6, 7)) %>%
  
  # retain focal variables
  transmute(
    datetime,
    
    AirTemp      = as.numeric(AirTC_1_Avg),
    RH           = as.numeric(RH_avg_1),
    Irradiance   = SlrW_Avg,
    SoilTemp     = SoilTC_2_Avg,
    SoilMoisture = SMwfv_2_Avg,
    
    TimeOfDay = as_hms(datetime)
  ) %>%
  
  # remove rows where all vars are NA
  filter(
    !(is.na(AirTemp) &
        is.na(RH) &
        is.na(Irradiance) &
        is.na(SoilTemp) &
        is.na(SoilMoisture))
  )

airstrip_summer <- airstrip_summer %>%
  mutate(
    TimeOfDay = as_hms(datetime)
  )

# LONG FORMAT FOR PLOTTING

airstrip_long <- airstrip_summer %>%
  pivot_longer(
    cols = c(AirTemp, RH, Irradiance, SoilTemp, SoilMoisture),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "Irradiance",
        "AirTemp",
        "SoilTemp",
        "RH",
        "SoilMoisture"
      )
    )
  )

airstrip_long <- airstrip_long %>%
  mutate(Date = as.Date(datetime))

airstrip_z <- airstrip_long %>%
  group_by(Date, Variable) %>%
  mutate(
    Z = (Value - mean(Value, na.rm = TRUE)) /
      sd(Value, na.rm = TRUE)
  ) %>%
  ungroup()

daily_profiles <- airstrip_z %>%
  distinct(Date, TimeOfDay, Variable, Z)

airstrip_summary <- daily_profiles %>%
  group_by(Variable, TimeOfDay) %>%
  summarize(
    MeanZ = mean(Z, na.rm = TRUE),
    SD    = sd(Z, na.rm = TRUE),
    N     = n_distinct(Date),
    SE    = SD / sqrt(N),
    LowerSE = MeanZ - SE,
    UpperSE = MeanZ + SE,
    LowerSD = MeanZ - SD,
    UpperSD = MeanZ + SD,
    .groups = "drop"
  )

fig_se <- ggplot(
  airstrip_summary,
  aes(
    x = TimeOfDay,
    y = MeanZ,
    color = Variable,
    fill = Variable,
    group = Variable
  )
) +
  geom_ribbon(
    aes(ymin = LowerSE, ymax = UpperSE),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Time of day",
    y = "Standardized value (z-score)",
    color = "Variable:",
    fill = "Variable:"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    legend.position = "bottom"
  )

# DISPLAY
fig_se

#### Stats Tables S1-S4 ####
#### Table S1 - Air model summaries ####
air_results <- bind_rows(
  
  extract_fixed_effects(
    ffair23_mod1,
    "Foraging flow 2023"
  ),
  
  extract_fixed_effects(
    ffair24_mod1,
    "Foraging flow 2024"
  ),
  
  extract_fixed_effects(
    ssair23_mod1,
    "Surface traffic 2023"
  ),
  
  extract_fixed_effects(
    ssair24_mod1,
    "Surface traffic 2024"
  )
)

term_labels <- c(
  
  "(Intercept)" =
    "Intercept",
  
  "IrradianceAvg" =
    "Irradiance",
  
  "MorningOrEveningevening" =
    "Evening context",
  
  "AirTemp" =
    "Air temperature",
  
  "RH" =
    "Relative humidity",
  
  "IrradianceAvg:MorningOrEveningevening" =
    "Irradiance × Evening",
  
  "MorningOrEveningevening:AirTemp" =
    "Evening × Air temperature",
  
  "MorningOrEveningevening:RH" =
    "Evening × RH"
)


air_results <- air_results %>%
  
  mutate(
    term = dplyr::recode(term, !!!term_labels)
  )

air_results <- air_results %>%
  
  mutate(
    
    estimate =
      round(estimate, 3),
    
    statistic =
      round(statistic, 2),
    
    p.value =
      case_when(
        p.value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p.value)
      ),
    
    result =
      paste0(
        "β = ", estimate,
        ", ",
        
        ifelse(
          str_detect(model, "Surface"),
          "z = ",
          "t = "
        ),
        
        statistic,
        ", p = ",
        p.value
      )
  )



air_table <- air_results %>%
  
  dplyr::select(
    term,
    model,
    result
  ) %>%
  
  pivot_wider(
    names_from = model,
    values_from = result
  )


ft_air <- flextable(air_table)


ft_air <- ft_air %>%
  
  autofit() %>%
  
  theme_booktabs() %>%
  
  fontsize(size = 10) %>%
  
  bold(part = "header") %>%
  
  align(
    align = "left",
    part = "all"
  ) %>%
  
  set_caption(
    caption =
      "Table S1. Fixed effects from air-environment mixed models."
  )



save_as_docx(
  
  "Table S1" = ft_air,
  
  path = "Table_S1_air_models.docx"
)



#### Table S2 - Soil model summaries ####
soil_results <- bind_rows(
  
  extract_fixed_effects(
    ffsoil23_mod1,
    "Foraging flow 2023"
  ),
  
  extract_fixed_effects(
    ffsoil24_modT,
    "Foraging flow 2024"
  ),
  
  extract_fixed_effects(
    sssoil23_mod1,
    "Surface traffic 2023"
  ),
  
  extract_fixed_effects(
    sssoil24_mod1,
    "Surface traffic 2024"
  )
)

soil_term_labels <- c(
  
  "(Intercept)" =
    "Intercept",
  
  "SoilTemp" =
    "Soil temperature",
  
  "SoilMoisture" =
    "Soil moisture",
  
  "MorningOrEveningevening" =
    "Evening context",
  
  "SoilTemp:MorningOrEveningevening" =
    "Soil temperature × Evening",
  
  "MorningOrEveningevening:SoilMoisture" =
    "Evening × Soil moisture"
)

soil_results <- soil_results %>%
  
  mutate(
    term = dplyr::recode(term, !!!soil_term_labels)
  )

# keep only shared terms
shared_terms <- soil_results %>%
  dplyr::count(term) %>%
  dplyr::filter(n == length(unique(soil_results$model))) %>%
  dplyr::pull(term)

soil_results <- soil_results %>%
  dplyr::filter(term %in% shared_terms)

# format
soil_results <- soil_results %>%
  
  mutate(
    
    estimate = round(estimate, 3),
    
    statistic = round(statistic, 2),
    
    p.value = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    ),
    
    result = paste0(
      "β = ", estimate,
      ", ",
      
      ifelse(
        str_detect(model, "Surface"),
        "z = ",
        "t = "
      ),
      
      statistic,
      ", p = ",
      p.value
    )
  )

soil_table <- soil_results %>%
  
  dplyr::select(
    term,
    model,
    result
  ) %>%
  
  pivot_wider(
    names_from = model,
    values_from = result
  )

ft_soil <- flextable(soil_table)

ft_soil <- ft_soil %>%
  
  autofit() %>%
  
  theme_booktabs() %>%
  
  fontsize(size = 10) %>%
  
  bold(part = "header") %>%
  
  align(
    align = "left",
    part = "all"
  ) %>%
  
  set_caption(
    caption =
      "Table S2. Fixed effects from soil-environment mixed models."
  )

save_as_docx(
  
  "Table S2" = ft_soil,
  
  path = "Table_S2_soil_models.docx"
)

#### Table S3 - post-hoc trends ####
library(dplyr)
library(tidyr)
library(flextable)

extract_trends <- function(x, response, predictor){
  
  df <- as.data.frame(x)
  
  trend_col <- grep("\\.trend$", names(df), value = TRUE)
  
  stat_col <- if("t.ratio" %in% names(df)) "t.ratio" else "z.ratio"
  
  df %>%
    rename(
      Estimate = all_of(trend_col),
      Statistic = all_of(stat_col)
    ) %>%
    mutate(
      Response = response,
      Predictor = predictor,
      Test = ifelse(stat_col == "t.ratio", "t", "z")
    )
}

trend_table <- bind_rows(
  extract_trends(ffair23_trendIrr,       "Foraging flow (2023)", "Irradiance"),
  extract_trends(ffsoil23_trendLux,      "Foraging flow (2023)", "Lux"),
  extract_trends(ffsoil23_trendSoilTemp, "Foraging flow (2023)", "Soil temperature"),
  
  extract_trends(ffair24_trendIrr,       "Foraging flow (2024)", "Irradiance"),
  extract_trends(ffsoil24_trendLux,      "Foraging flow (2024)", "Lux"),
  extract_trends(ffsoil24_trendSoilTemp, "Foraging flow (2024)", "Soil temperature"),
  
  extract_trends(sssoil23_trendLux,      "Steady-state traffic (2023)", "Lux"),
  extract_trends(sssoil23_trendSoilTemp, "Steady-state traffic (2023)", "Soil temperature")
) %>%
  mutate(
    stat_label = ifelse(Test == "t", "t", "z"),
    
    result = paste0(
      "β = ", round(Estimate, 3),
      ", ", stat_label, " = ", round(Statistic, 2),
      ", p = ",
      ifelse(
        p.value < 0.001,
        "<0.001",
        sprintf("%.3f", p.value)
      )
    )
  ) %>%
  select(Response, Predictor, MorningOrEvening, result) %>%
  pivot_wider(
    names_from = MorningOrEvening,
    values_from = result
  ) %>%
  rename(
    Term = Predictor
  )

ft_trends <- flextable(trend_table)

ft_trends <- ft_trends %>%
  
  autofit() %>%
  
  theme_booktabs() %>%
  
  fontsize(size = 10) %>%
  
  bold(part = "header") %>%
  
  align(
    align = "left",
    part = "all"
  ) %>%
  
  set_caption(
    caption =
      "Table S3. Post-hoc marginal trends (simple slopes) for significant environmental predictor × context interactions in foraging flow and steady state foraging models."
  )

save_as_docx(
  
  "Table S3" = ft_trends,
  
  path = "Table_S3_posthoc_trends.docx"
)


#### Table S4 - AIC Comparison Table ####
aic_table <- tibble(
  
  model = c(
    "Foraging flow 2024 (air)",
    "Foraging flow 2024 (soil)",
    "Steady state foraging 2024 (air)",
    "Steady state foraging (soil)"
  ),
  
  AIC_base = c(
    AIC(ffair24_mod1),
    AIC(ffsoil24_mod1),
    AIC(ssair24_mod1),
    AIC(sssoil24_mod1)
  ),
  
  AIC_treatment = c(
    AIC(ffair24_modT),
    AIC(ffsoil24_modT),
    AIC(ssair24_modT),
    AIC(sssoil24_modT)
  )
) %>%
  
  mutate(
    
    AIC_base =
      round(AIC_base, 2),
    
    AIC_treatment =
      round(AIC_treatment, 2),
    
    delta_AIC =
      round(
        AIC_treatment - AIC_base,
        2
      ),
    
    preferred_model =
      ifelse(
        delta_AIC < 0,
        "Treatment",
        "Base"
      )
  )

ft_aic <- flextable(aic_table)

ft_aic <- ft_aic %>%
  
  autofit() %>%
  
  theme_booktabs() %>%
  
  fontsize(size = 10) %>%
  
  bold(part = "header") %>%
  
  align(
    align = "left",
    part = "all"
  ) %>%
  
  set_caption(
    caption =
      "Table S4. Comparison of models with and without treatment effects. Negative ΔAIC values indicate improved model fit when treatment was included."
  )

save_as_docx(
  
  "Table S4" = ft_aic,
  
  path = "Table_S4_AIC_comparison.docx"
)

mutate(
  preferred_model =
    ifelse(
      delta_AIC < 0,
      "Treatment",
      "Base"
    )
)

#### ####
