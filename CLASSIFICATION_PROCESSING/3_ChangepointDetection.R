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

#use data from step 1_VegTrend
data <- data

##filter just the dead trees which are class 4 and group by pixel per month
df <- data |>
  group_by(.geo, year_month_1) 

# Assuming 'geometry' identifies the location, filter by the first value of 'geometry'
df_2 <- df[df$.geo == df$.geo[1], ]

#Convert the extracted date to Date class in R
df_2$year_month <- as.Date(df_2$year_month, format = "%Y%m")



############# DBEST MODEL ###################

time_series <- ts(df_2$CCCI, start=c(2019, 1), end = c(2024, 12), frequency = 12)  # assuming monthly data from 2021

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
  CCCI = as.numeric(DBEST.bettle$Data),
  Fit = as.numeric(DBEST.bettle$Fit)
)

# Extract breakpoints (use Start or End indices)
breakpoints <- DBEST.bettle$End  # or DBEST.bettle$Start if you prefer

# Plot
ggplot(df, aes(x = Time, y = CCCI)) +
  geom_line(color = "forestgreen", size = 1) +
  geom_line(aes(y = Fit), color = "blue", linetype = "dashed", size = 1) +
  geom_vline(xintercept = df$Time[breakpoints], linetype = "dashed", color = "red") +
  labs(
    title = "Detected Breakpoints in Beetle Time Series Using DBEST",
    x = "Time (Years)",
    y = "Canopy Chlorophyll Content Index (CCCI)"
  ) +
  theme_minimal(base_size = 13) +
 theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face="bold", color = "black")
  ) +  
  annotate(
    "text", x = min(df$Time)+0, y = max(df$CCCI)*1.2,
    label = "Green: CCCI (Observed)\nBlue dashed: DBEST Fit\nRed dashed: Change points",
    hjust = 0, vjust = 1, size = 3.8, color = "black", face = "bold"
  ) 



######################################################


############################### BEAST MODEL ####################################

time_series <- ts(df_2$CCCI, start=c(2019, 1), end = c(2024, 12), frequency = 12)  # assuming monthly data from 2021

# Interpolate missing values
time_series_imputed <- na.approx(time_series)


BEAST_BETTLE <- beast(time_series_imputed, # start  = c(2019,01,01),deltat = 1/365,  # period = 365      
          season  = 'harmonic')  
print(BEAST_BETTLE)
str(BEAST_BETTLE)


# Convert BEAST output to data frame
df_beast <- data.frame(
  Time = BEAST_BETTLE$time,
  CCCI = BEAST_BETTLE$data,
  Trend = BEAST_BETTLE$trend$Y
)

# Extract all trend changepoints and remove NAs
trend_cps_all <- sort(BEAST_BETTLE$trend$cp[!is.na(BEAST_BETTLE$trend$cp)])

# Select last 3
trend_cps <- tail(trend_cps_all, 3)
trend_cps

# Create the plot
ggplot(df_beast, aes(x = Time, y = CCCI)) +
  geom_line(color = "forestgreen", size = 1) +  # Observed NDWI
  geom_line(aes(y = Trend), color = "blue", linetype = "dashed", size = 1) +  # Trend fit
  geom_vline(xintercept = trend_cps, linetype = "dashed", color = "red") +  # Changepoints
  annotate(
    "text", x = min(df_beast$Time) + 1.4, y = max(df_beast$CCCI)*1.1,
    label = "Green: CCCI (Observed)\nBlue dashed: BEAST Fit\nRed dashed: Change points",
    hjust = 0, vjust = 1, size = 3.8, color = "black", face = "bold"
  ) +
  labs(
    title = "Detected Breakpoints in Beetle Time Series Using BEAST",
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


# Prepare data in long format
df_plot <- data.frame(
  Time = as.numeric(time(time_series_imputed)),
  CCCI = as.numeric(time_series_imputed),
  Trend = as.numeric(BFAST_BETTLE$output[[1]]$Tt)
)

#plot the last 3 breakpoints
bp_indices <- tail(BFAST_BETTLE$output[[1]]$Vt.bp, 3)
bp_times   <- df_plot$Time[bp_indices] # convert indices to actual time

# Create the plot
ggplot(df_plot, aes(x = Time, y = CCCI)) +
  geom_line(color = "forestgreen", size = 1) +  # Observed NDWI
  geom_line(aes(y = Trend), color = "blue", linetype = "dashed", size = 1) +  # Trend fit
  geom_vline(xintercept = bp_times, linetype = "dashed", color = "red") +  # Changepoints
  annotate(
   "text", x = min(df_plot$Time) + 0, y = max(df_plot$CCCI)*0.9,
    label = "Green: CCCI (Observed)\nBlue dashed: BFAST Fit\nRed dashed: Change points",
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



