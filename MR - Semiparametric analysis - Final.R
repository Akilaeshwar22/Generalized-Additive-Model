# Final - Rubella 

library(mgcv)
library(dplyr)
library(ggplot2)
library(readxl)

# --- Step 1: Read data ---
df_r <- read_excel("D:\\CMC-WTRL\\MR\\RUBELLA\\For analysis\\Semiparametric analysis - Rubella.xlsx", sheet = "Sheet1")

df_r <- df_r %>%
  mutate(Seropositive = ifelse(Seropositive == "Seropositive", 1, 0))

# --- Step 2: Create age groups 
df_r <- df_r %>%
  mutate(AgeGroup = cut(
    Age,
    breaks = seq(0, 80, by = 5),   # every 5 years
    right = TRUE, include.lowest = TRUE,
    labels = paste0(seq(0, 75, by = 5), "-", seq(5, 80, by = 5))
  ))


# --- Step 3: Fit GAM model ---
gam_fit_r <- gam(Seropositive ~ s(Age, by = Year) + Year,
               data = df_r, family = binomial)

# --- Step 4: Predictions over continuous age ---
newdata_r <- expand.grid(
  Age = seq(0, 80, by = 1),
  Year = unique(df_r$Year)
)
pred_r <- predict(gam_fit_r, newdata_r, type = "link", se.fit = TRUE)
newdata_r$fit <- plogis(pred_r$fit)
newdata_r$lwr <- plogis(pred_r$fit - 1.96 * pred_r$se.fit)
newdata_r$upr <- plogis(pred_r$fit + 1.96 * pred_r$se.fit)

# --- Step 5: Calculate observed prevalence per AgeGroup ---
obs_summary_r <- df_r %>%
  group_by(Year, AgeGroup) %>%
  summarise(
    n = n(),
    prev = mean(Seropositive, na.rm = TRUE),
    .groups = "drop"
  )

# --- Step 6: Calculate observed prevalence per age group & year
df_grouped_r <- df_r %>%
  group_by(Year, AgeGroup) %>%
  summarise(prevalence = mean(Seropositive), .groups = "drop")

# --- Step 7: Plot ---
p1 <- ggplot() +
  # Smoothed GAM fit
  geom_line(data = newdata_r, aes(x = Age, y = fit, colour = factor(Year)), size = 1) +
  geom_ribbon(data = newdata_r, aes(x = Age, ymin = lwr, ymax = upr, fill = factor(Year)),
              alpha = 0.2, colour = NA) +
  # Observed prevalence points (midpoint of each group)
  geom_point(data = df_grouped_r,
             aes(x = as.numeric(sub("-.*", "", AgeGroup)) + 
                   (as.numeric(sub(".*-", "", AgeGroup)) - 
                      as.numeric(sub("-.*", "", AgeGroup)))/2,
                 y = prevalence, colour = factor(Year)),
             size = 0.6, shape = 21, stroke = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Age group (years)",
    y = "Proportion of participants\n seropositive for rubella (%)",  # extra \n for spacing
    colour = "Year",
    fill = "Year"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13, face = "bold", colour = "black"),
    axis.text = element_text(size = 12, face = "bold", colour = "black"),
    legend.title = element_text(size = 13, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, face = "bold", colour = "black")
  )

# --- Step 7: Save as TIFF ---
tiff("D:\\CMC-WTRL\\MR\\MR graphs\\Rubella_GAM.tiff", width = 8, height = 6, units = "in", res = 300)
print(p1)
dev.off()

#######################################################################################################################################


# Final - Measles 

library(mgcv)
library(dplyr)
library(ggplot2)
library(readxl)

# --- Step 1: Read data ---
df_m <- read_excel("D:\\CMC-WTRL\\MR\\MEASLES\\For analysis\\Semiparametric analysis - Measles.xlsx", sheet = "Sheet1")

df_m <- df_m %>%
  mutate(Seropositive = ifelse(Seropositive == "Seropositive", 1, 0))

# --- Step 2: Create age groups 
df_m <- df_m %>%
  mutate(AgeGroup = cut(
    Age,
    breaks = seq(0, 80, by = 5),   # every 5 years
    right = TRUE, include.lowest = TRUE,
    labels = paste0(seq(0, 75, by = 5), "-", seq(5, 80, by = 5))
  ))


# --- Step 3: Fit GAM model ---
gam_fit_m <- gam(Seropositive ~ s(Age, by = Year) + Year,
               data = df_m, family = binomial)

# --- Step 4: Predictions over continuous age ---
newdata_m <- expand.grid(
  Age = seq(0, 80, by = 1),
  Year = unique(df_m$Year)
)
pred_m <- predict(gam_fit_m, newdata_m, type = "link", se.fit = TRUE)
newdata_m$fit <- plogis(pred_m$fit)
newdata_m$lwr <- plogis(pred_m$fit - 1.96 * pred_m$se.fit)
newdata_m$upr <- plogis(pred_m$fit + 1.96 * pred_m$se.fit)

# --- Step 5: Calculate observed prevalence per AgeGroup ---
obs_summary_m <- df_m %>%
  group_by(Year, AgeGroup) %>%
  summarise(
    n = n(),
    prev = mean(Seropositive, na.rm = TRUE),
    .groups = "drop"
  )

# --- Step 6: Calculate observed prevalence per age group & year
df_grouped_m <- df_m %>%
  group_by(Year, AgeGroup) %>%
  summarise(prevalence = mean(Seropositive), .groups = "drop")

# --- Step 7: Plot ---
p2 <- ggplot() +
  # Smoothed GAM fit
  geom_line(data = newdata_m, aes(x = Age, y = fit, colour = factor(Year)), size = 1) +
  geom_ribbon(data = newdata_m, aes(x = Age, ymin = lwr, ymax = upr, fill = factor(Year)),
              alpha = 0.2, colour = NA) +
  # Observed prevalence points (midpoint of each group)
  geom_point(data = df_grouped_m,
             aes(x = as.numeric(sub("-.*", "", AgeGroup)) + 
                   (as.numeric(sub(".*-", "", AgeGroup)) - 
                      as.numeric(sub("-.*", "", AgeGroup)))/2,
                 y = prevalence, colour = factor(Year)),
             size = 0.6, shape = 21, stroke = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Age group (years)",
    y = "Proportion of participants\n seropositive for measles (%)",  # extra \n for spacing
    colour = "Year",
    fill = "Year"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none" ,
    axis.title = element_text(size = 13, face = "bold", colour = "black"),
    axis.text = element_text(size = 12, face = "bold", colour = "black"),
    legend.title = element_text(size = 13, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, face = "bold", colour = "black")
  )

# --- Step 7: Save as TIFF ---
tiff("D:\\CMC-WTRL\\MR\\MR graphs\\Measles_GAM.tiff", width = 8, height = 6, units = "in", res = 300)
print(p2)
dev.off()

summary(gam_fit_m)



library (patchwork)

# Combine plots with shared legend at the bottom
p_combined <- p1+ p2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Save as TIFF
tiff("D:\\CMC-WTRL\\MR\\MR graphs\\Rubella and measles_GAM.tiff",
     width = 14, height = 6, units = "in", res = 300)
print(p_combined)
dev.off()


####################################################################################################################

# Sample size is reflected in the size of the dots - Rubella 

library(mgcv)
library(dplyr)
library(ggplot2)
library(readxl)

# --- Step 1: Read data ---
df_r <- read_excel("D:\\CMC-WTRL\\MR\\RUBELLA\\For analysis\\Semiparametric analysis - Rubella.xlsx", sheet = "Sheet1")

df_r <- df_r %>%
  mutate(Seropositive = ifelse(Seropositive == "Seropositive", 1, 0))

# --- Step 2: Create age groups 
df_r <- df_r %>%
  mutate(AgeGroup = cut(
    Age,
    breaks = seq(0, 80, by = 5),   # every 5 years
    right = TRUE, include.lowest = TRUE,
    labels = paste0(seq(0, 75, by = 5), "-", seq(5, 80, by = 5))
  ))


# --- Step 3: Fit GAM model ---
gam_fit_r <- gam(Seropositive ~ s(Age, by = Year) + Year,
                 data = df_r, family = binomial)

# --- Step 4: Predictions over continuous age ---
newdata_r <- expand.grid(
  Age = seq(0, 80, by = 1),
  Year = unique(df_r$Year)
)
pred_r <- predict(gam_fit_r, newdata_r, type = "link", se.fit = TRUE)
newdata_r$fit <- plogis(pred_r$fit)
newdata_r$lwr <- plogis(pred_r$fit - 1.96 * pred_r$se.fit)
newdata_r$upr <- plogis(pred_r$fit + 1.96 * pred_r$se.fit)

# --- Step 5: Calculate observed prevalence per AgeGroup ---
obs_summary_r <- df_r %>%
  group_by(Year, AgeGroup) %>%
  summarise(
    n = n(),
    prev = mean(Seropositive, na.rm = TRUE),
    .groups = "drop"
  )

# --- Step 6: Calculate observed prevalence per age group & year
df_grouped_r <- df_r %>%
  group_by(Year, AgeGroup) %>%
  summarise(
    prevalence = mean(Seropositive),
    n = n(),   # sample size
    .groups = "drop"
  )


# --- Step 7: Plot ---
p1 <- ggplot() +
  # Smoothed GAM fit
  geom_line(data = newdata_r, aes(x = Age, y = fit, colour = factor(Year)), size = 1) +
  geom_ribbon(data = newdata_r, aes(x = Age, ymin = lwr, ymax = upr, fill = factor(Year)),
              alpha = 0.2, colour = NA) +
  # Observed prevalence points (midpoint of each group, size = n)
  geom_point(
    data = df_grouped_r,
    aes(
      x = as.numeric(sub("-.*", "", AgeGroup)) +
        (as.numeric(sub(".*-", "",AgeGroup)) -
           as.numeric(sub("-.*","",AgeGroup)))/2,
      y = prevalence,
      colour = factor(Year),   # this controls the interior fill for 15–20
      size = n
    ),
    shape = 16   # solid filled circle
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_size(range = c(1, 3)) +  # adjust size range
  labs(
    x = "Age group (years)",
    y = "Proportion of participants\nseropositive for rubella (%)",  
    colour = "Year",
    fill = "Year"
  ) +
  guides(size = "none") +  # remove size legend
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 13, face = "bold", colour = "black"),
    axis.text = element_text(size = 12, face = "bold", colour = "black"),
    legend.title = element_text(size = 13, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, face = "bold", colour = "black")
  )

# --- Step 7: Save as TIFF ---
tiff("D:\\CMC-WTRL\\MR\\MR graphs\\Rubella_GAM-2.tiff", width = 8, height = 6, units = "in", res = 300)
print(p1)
dev.off()

##############################################################################################################
# Sample size is reflected in the size of the dots - Measles 

library(mgcv)
library(dplyr)
library(ggplot2)
library(readxl)

# --- Step 1: Read data ---
df_m <- read_excel("D:\\CMC-WTRL\\MR\\MEASLES\\For analysis\\Semiparametric analysis - Measles.xlsx", sheet = "Sheet1")

df_m <- df_m %>%
  mutate(Seropositive = ifelse(Seropositive == "Seropositive", 1, 0))

# --- Step 2: Create age groups 
df_m <- df_m %>%
  mutate(AgeGroup = cut(
    Age,
    breaks = seq(0, 80, by = 5),   # every 5 years
    right = TRUE, include.lowest = TRUE,
    labels = paste0(seq(0, 75, by = 5), "-", seq(5, 80, by = 5))
  ))


# --- Step 3: Fit GAM model ---
gam_fit_m <- gam(Seropositive ~ s(Age, by = Year) + Year,
                 data = df_m, family = binomial)

# --- Step 4: Predictions over continuous age ---
newdata_m <- expand.grid(
  Age = seq(0, 80, by = 1),
  Year = unique(df_m$Year)
)
pred_m <- predict(gam_fit_m, newdata_m, type = "link", se.fit = TRUE)
newdata_m$fit <- plogis(pred_m$fit)
newdata_m$lwr <- plogis(pred_m$fit - 1.96 * pred_m$se.fit)
newdata_m$upr <- plogis(pred_m$fit + 1.96 * pred_m$se.fit)

# --- Step 5: Calculate observed prevalence per AgeGroup ---
obs_summary_m <- df_m %>%
  group_by(Year, AgeGroup) %>%
  summarise(
    n = n(),
    prev = mean(Seropositive, na.rm = TRUE),
    .groups = "drop"
  )

# --- Step 6: Calculate observed prevalence per age group & year
df_grouped_m <- df_m %>%
  group_by(Year, AgeGroup) %>%
  summarise(
    prevalence = mean(Seropositive),
    n = n(),   # sample size
    .groups = "drop"
  )


# --- Step 7: Plot ---
p2 <- ggplot() +
  # Smoothed GAM fit
  geom_line(data = newdata_m, aes(x = Age, y = fit, colour = factor(Year)), size = 1) +
  geom_ribbon(data = newdata_m, aes(x = Age, ymin = lwr, ymax = upr, fill = factor(Year)),
              alpha = 0.2, colour = NA) +
  # Observed prevalence points (midpoint of each group, size = n)
  geom_point(
    data = df_grouped_m,
    aes(
      x = as.numeric(sub("-.*", "", AgeGroup)) +
        (as.numeric(sub(".*-", "",AgeGroup)) -
           as.numeric(sub("-.*","",AgeGroup)))/2,
      y = prevalence,
      colour = factor(Year),   # this controls the interior fill for 15–20
      size = n
    ),
    shape = 16   # solid filled circle
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_size(range = c(1, 3)) +  # adjust size range
  labs(
    x = "Age group (years)",
    y = "Proportion of participants\nseropositive for rubella (%)",  
    colour = "Year",
    fill = "Year"
  ) +
  guides(size = "none") +  # remove size legend
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 13, face = "bold", colour = "black"),
    axis.text = element_text(size = 12, face = "bold", colour = "black"),
    legend.title = element_text(size = 13, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, face = "bold", colour = "black")
  )

# --- Step 7: Save as TIFF ---
tiff("D:\\CMC-WTRL\\MR\\MR graphs\\Measles_GAM-2.tiff", width = 8, height = 6, units = "in", res = 300)
print(p2)
dev.off()

summary(gam_fit_m)



library (patchwork)

# Combine plots with shared legend at the bottom
p_combined <- p1+ p2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Save as TIFF
tiff("D:\\CMC-WTRL\\MR\\MR graphs\\Rubella and measles_GAM-2.tiff",
     width = 14, height = 6, units = "in", res = 300)
print(p_combined)
dev.off()


####################################################################################################################
