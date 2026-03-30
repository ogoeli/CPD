library(sf)
library(sp)
library(tidyverse)
library(raster)
library(dplyr)
library(data.table)
library(terra)
library(jsonlite)
library(bfast)
library(Rbeast)
library(DBEST)


#data <- read.csv('/scratch/ope4/Downloads/Band_Values_By_Class.csv', quote = "\"")

data <- boundary



# Step 1: Extract the date part (first 8 characters from 'system.index')
data$timestamp <- substr(data$system.index, 1, 8)

# Step 2: Convert the extracted date to Date class in R
data$timestamp <- as.Date(data$timestamp, format = "%Y%m%d")

# Step 3: Extract the year and month from timestamp
data$year_month_1 <- format(data$timestamp, "%Y-%m")

# Extract year and month from the system_index column
data <- data %>%
  mutate(
    year = substr(system.index, 1, 4),
    month = substr(system.index, 5, 6),
    year_month = paste(year, month, sep = "-")
  )


# Calculate indices
data$NDVI <- (data$B8 - data$B4) / (data$B8 + data$B4)
data$NDWI <- (data$B8 - data$B11) / (data$B8 + data$B11)
data$RENDVI <- (data$B8 - data$B5) / (data$B8 + data$B5)
data$CCCI <- data$RENDVI/data$NDVI
data$CCCI = ifelse(data$NDVI == 0, NA, data$RENDVI / data$NDVI) 


df <- data |>
  group_by(.geo, year_month_1) |>
  # distinct(.geo, year_month_1, .keep_all = TRUE) |>
  filter(cover == 4) # | cover == 3




# Assuming 'geometry' identifies the location, filter by the first value of 'geometry'
df_2 <- df[df$.geo == df$.geo[1], ]

str(df_2)

# Step 2: Convert the extracted date to Date class in R
df_2$year_month <- as.Date(df_2$year_month, format = "%Y%m")



############# DBEST MODEL ###################

time_series <- ts(df_2$NDWI, start=c(2019, 1), end = c(2024, 12), frequency = 12)  # assuming monthly data from 2021

# Interpolate missing values
time_series_imputed <- na.approx(time_series)



DBEST.bettle <- DBEST(data=time_series_imputed, data.type="cyclical", 
                     algorithm="change detection", 
                     breakpoints.no=3, first.level.shift=0.1, 
                     second.level.shift=0.2, duration=148, 
                     distance.threshold="default", alpha=0.05, plot="off")

print(DBEST.bettle)
str(DBEST.bettle)

# Convert your time series to a data frame
df <- data.frame(
  Time = time(DBEST.bettle$Data),
  NDWI = as.numeric(DBEST.bettle$Data),
  Fit = as.numeric(DBEST.bettle$Fit)
)

# Extract breakpoints (use Start or End indices)
breakpoints <- DBEST.bettle$End  # or DBEST.bettle$Start if you prefer

# Plot
ggplot(df, aes(x = Time, y = NDWI)) +
  geom_line(color = "forestgreen", size = 1) +
  geom_line(aes(y = Fit), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(xintercept = df$Time[breakpoints], linetype = "dashed", color = "red") +
  labs(
    title = "Detected Breakpoints in Beetle Time Series Using DBEST",
    x = "Time (Years)",
    y = "Normalized Difference Water Index (NDWI)"
  ) +
  theme_minimal(base_size = 13) +
 theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face="bold", color = "black")
  ) +  
  annotate(
    "text", x = min(df$Time)+2.2, y = max(df$NDWI)*1.0,
    label = "Green: NDWI (Observed)\nBlue dashed: Trend (DBEST Fit)\nRed dashed: Change points",
    hjust = 0, vjust = 1, size = 3.8, color = "black", face = "bold"
  ) 



######################################################


############################### BEAST MODEL ####################################



time_series <- ts(df_2$NDVI, start=c(2019, 1), end = c(2024, 12), frequency = 12)  # assuming monthly data from 2021

# Interpolate missing values
time_series_imputed <- na.approx(time_series)


BEAST_BETTLE <- beast(time_series_imputed, # start  = c(2019,01,01),deltat = 1/365,  # period = 365      
          season  = 'harmonic')  
print(BEAST_BETTLE)
str(BEAST_BETTLE)


# Convert BEAST output to data frame
df_beast <- data.frame(
  Time = BEAST_BETTLE$time,
  NDVI = BEAST_BETTLE$data,
  Trend = BEAST_BETTLE$trend$Y
)

# Extract all trend changepoints and remove NAs
trend_cps_all <- sort(BEAST_BETTLE$trend$cp[!is.na(BEAST_BETTLE$trend$cp)])

# Select last 3
trend_cps <- tail(trend_cps_all, 3)
trend_cps

# Create the plot
ggplot(df_beast, aes(x = Time, y = NDVI)) +
  geom_line(color = "forestgreen", size = 1) +  # Observed NDWI
  geom_line(aes(y = Trend), color = "blue", linetype = "dashed", size = 1) +  # Trend fit
  geom_vline(xintercept = trend_cps, linetype = "dashed", color = "red") +  # Changepoints
  annotate(
    "text", x = min(df_beast$Time) + 0.3, y = max(df_beast$NDVI)*1.0,
    label = "Green: NDVI (Observed)\nBlue dashed: Trend (BEAST Fit)\nRed dashed: Change points",
    hjust = 0, vjust = 1, size = 3.8, color = "black", face = "bold"
  ) +
  labs(
    title = "Detected Breakpoints in Beetle Time Series Using BEAST",
    x = "Time (Years)",
    y = "Normalized Difference Vegetation Index (NDVI)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face="bold", color = "black")
  )


############################################################################################



#################################### BFAST MODE ###########################


time_series <- ts(df_2$CCCI, start=c(2019, 1), end = c(2024, 12), frequency = 12)  # assuming monthly data from 2021

# Interpolate missing values
time_series_imputed <- na.approx(time_series)


# ratio of distance between breaks (time steps) and length of the time series
BFAST_BETTLE <- bfast(time_series_imputed, season = "harmonic")
BFAST_BETTLE
str(BFAST_BETTLE)

# Extract observed and trend series
obs <- as.numeric(BFAST_BETTLE$Yt)
trend <- as.numeric(BFAST_BETTLE$output[[1]]$Tt)  # trend component from first segment

# Time vector
time <- as.numeric(time(time_series_imputed))


# Breakpoints: last 3 if needed
trend_cps <- tail(BFAST_BETTLE$output[[1]]$Vt.bp, 3)

# Prepare data in long format
df_plot <- data.frame(
  Time = as.numeric(time(time_series_imputed)),
  CCCI = as.numeric(time_series_imputed),
  Trend = as.numeric(BFAST_BETTLE$output[[1]]$Tt)
)

bp_indices <- unlist(BFAST_BETTLE$output[[1]]$Vt.bp)
bp_times   <- df_plot$Time[bp_indices]  # convert indices to actual time


# Create the plot
ggplot(df_plot, aes(x = Time, y = CCCI)) +
  geom_line(color = "forestgreen", size = 1) +  # Observed NDWI
  geom_line(aes(y = Trend), color = "blue", linetype = "dashed", size = 1) +  # Trend fit
  geom_vline(xintercept = bp_times, linetype = "dashed", color = "red") +  # Changepoints
  annotate(
   "text", x = min(df_plot$Time) + 0, y = max(df_plot$CCCI)*1.0,
    label = "Green: CCCI (Observed)\nBlue dashed: Trend (BFAST Fit)\nRed dashed: Change points",
    hjust = 0, vjust = 1, size = 3.8, color = "black", face = "bold"
  ) +
  labs(
    title = "Detected Breakpoints in Beetle Time Series Using BFAST",
    x = "Time (Years)",
    y = "Canopy Chlorophyll Content Index (CCCI)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face="bold", color = "black")
  )














################################################################################
#BFAST FOR CCCI WHEN IT IS DIFFICULT TO GET BREAKPOINTS

###############################################
# 1. CLEAN + PREPARE DF
###############################################
df_clean <- df %>%
  filter(cover == 4) %>%                  # target cover class
  filter(.geo == .geo[1]) %>%             # use one location
  mutate(
    year_month = paste0(year_month_1, "-01"),
    year_month = as.Date(year_month)
  ) %>%
  group_by(year_month) %>%
  summarise(
    CCCI = mean(CCCI, na.rm = TRUE)
  ) %>%
  arrange(year_month)

###############################################
# 2. CREATE MONTHLY TIME SERIES STARTING 2019-01
###############################################
full_months <- seq.Date(from = as.Date("2019-01-01"),
                        to   = max(df_clean$year_month),
                        by   = "month")

df_full <- data.frame(year_month = full_months) %>%
  left_join(df_clean, by = "year_month")

time_series <- ts(df_full$CCCI,
                  start = c(2019, 1),
                  frequency = 12)

# Interpolate missing months
time_series <- na.approx(time_series)

###############################################
# 3. RUN BFAST
###############################################
BFAST_BETTLE <- bfast(
  time_series,
  season = "harmonic",
  h = 0.05,
  max.iter = 10
)

###############################################
# 4. EXTRACT BREAKPOINTS + TREND
###############################################
trend <- as.numeric(BFAST_BETTLE$output[[1]]$Tt)
bp_index <- BFAST_BETTLE$output[[1]]$Vt.bp
bp_index <- bp_index[bp_index > 0]  # keep positive indices
bp_dates <- full_months[bp_index]

# Select last 3 break points
bp_dates <- bp_dates[3:5]
bp_dates


df_plot <- data.frame(
  Time = full_months,
  CCCI = as.numeric(time_series),
  Trend = trend
)

###############################################
# 5. PLOT WITH TREND, BREAKPOINTS, AND ANNOTATIONS
###############################################
ggplot(df_plot, aes(x = Time, y = CCCI)) +
  geom_line(color = "forestgreen", size = 1) +                       # Observed
  geom_line(aes(y = Trend), color = "blue", linetype = "dashed", size = 1) +  # Trend
  geom_vline(xintercept = bp_dates, linetype = "dashed", color = "red") +    # Breakpoints
  annotate(
    "text", x = min(df_plot$Time), y = max(df_plot$CCCI)*1.05,
    label = "Green: CCCI (Observed)\nBlue dashed: Trend (BFAST Fit)\nRed dashed: change points",
    hjust = 0, vjust = 1, size = 3.8, color = "black"
  ) +
  labs(
    title = "Detected Breakpoints in Beetle Time Series Using BFAST",
    x = "Time (Years)",
    y = "Canopy Chlorophyll Content Index (CCCI)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(face = "bold", color = "black")
  ) + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") 
