# Load necessary libraries
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


#load in the ground truth file
data <- read.csv('/scratch/ope4/CPD_PAPER/Band_Values_POLY_Monica.csv') ##do not apply tree mask

# Extract coordinates and split into longitude and latitude
coords <- t(sapply(data$.geo, function(x) {
  coords <- fromJSON(x)$coordinates
  return(coords)
}))

# Assign longitude and latitude to separate columns
data$longitude <- coords[, 1]
data$latitude <- coords[, 2]

# View the new columns
head(data[, c("longitude", "latitude")])

#get a shapefile to find boundary
boundary <- read_sf("/scratch/ope4/POLY_2/POLY_2/POLY_2.shp")

# Ensure your 'data' dataframe has 'longitude' and 'latitude' columns
data_sf <- st_as_sf(data, coords = c("longitude", "latitude"), crs = st_crs(boundary))

# Assign the CRS of the boundary to data_sf
data_sf <- st_set_crs(data_sf, st_crs(boundary)) 

# Reproject data_sf to match boundary's CRS
data_sf <- st_transform(data_sf, st_crs(boundary))

#apply tree mask from GEE
tree_mask <- brick("/scratch/ope4/MERGE/TREE_MASK/TREE_MASK.TIF")

# Check CRS of points and raster
crs_points <- st_crs(data_sf)  
crs_raster <- crs(tree_mask)

if (!identical(crs_points, crs_raster)) {
  # Transform points to match raster CRS if they differ
  points <- st_transform(data_sf, crs_raster)  
}

# Use nearest neighbor interpolation to extract tree mask values at each Sentinel 2 pixel
extracted_values <- raster::extract(tree_mask, points)  

#number of points for tree mask and actual sentinel 2 pixels
nrow(extracted_values)
nrow(data_sf)

# Ensure extracted_values is a data frame
extracted_values_2 <- as.data.frame(extracted_values)

#bind both sentinel 2 pixels to tree mask values such that each pixel have tree height
extracted_values_15 <- cbind(extracted_values_2, data_sf)

# Now, convert the dataframe to an sf object and assign CRS, 
boundary <- extracted_values_15 

# Step 1: Extract the date part (first 8 characters from 'system.index')
boundary$timestamp <- substr(boundary$system.index, 1, 8)

# Step 2: Convert the extracted date to Date class in R
boundary$timestamp <- as.Date(boundary$timestamp, format = "%Y%m%d")

# Step 3: Extract the year and month from timestamp
boundary$year_month_1 <- format(boundary$timestamp, "%Y-%m")

# Extract year and month from the system_index column
data <- boundary %>%
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
summary(data)


#------------------------------------ Decomposed Indices Plotting 

# Group by year_month and cover, and calculate the mean for each band
grouped_data <- data %>%
  group_by(year_month, cover) %>%
  summarise(across(matches('^(N|CC)'), mean, na.rm = TRUE)) 

##select year
grouped_data$year_month <- as.Date(paste0(grouped_data$year_month, "-01"))


# Reshape the data to long format focusing only on the indices
indices_data <- grouped_data %>%
  pivot_longer(cols = c(NDWI, NDVI, CCCI),  
               names_to = 'Index', 
               values_to = 'Value') 

indices_data <- indices_data %>%
  mutate(year_month = as.Date(paste0(year_month, "-01")))


# Create a complete monthly sequence per Index
indices_data_filled <- indices_data %>%
  group_by(Index) %>%
  complete(year_month = seq(min(year_month), max(year_month), by = "month")) %>%
  arrange(Index, year_month) %>%
  mutate(
    Value = na.approx(Value, na.rm = FALSE),
    Value = ifelse(is.infinite(Value) | is.nan(Value), NA, Value)  
  ) %>%
  mutate(Value = na.approx(Value, na.rm = FALSE)) %>%  
  ungroup()

## decomposed the original dataset
decomposed_results <- indices_data_filled %>%
  split(.$Index) %>%
  map(~ {
    x <- .x %>% arrange(year_month)
    ts_data <- ts(x$Value, start = c(year(x$year_month[1]), month(x$year_month[1])), frequency = 12)
    decompose(ts_data, type = "additive")
  })

#from the decomposed dataset, select an index to plot
decomposed_results_1 <- decomposed_results$CCCI

# Convert the index time series trend to a data frame
df_trend <- data.frame(
  Time = as.numeric(time(decomposed_results_1$trend)),
  Trend = as.numeric(decomposed_results_1$trend)
)

# Convert time to Date
df_trend$Date <- as.Date(as.yearmon(time(decomposed_results_1$trend)))

# Attack dates
first_date <- as.Date("2023-06-17")
last_date  <- as.Date("2024-08-15")

# ggplot
ggplot(df_trend, aes(x = Date, y = Trend)) +
  geom_line(color = "black", size = 1.2) +
  labs(
    title = "Trend Component of CCCI",
    x = "Year",
    y = "Trend Value"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face="bold", color = "black")
  ) +
  geom_vline(xintercept = first_date, linetype = "dashed", color = "red") +
  geom_vline(xintercept = last_date,  linetype = "dashed", color = "red") +
  scale_x_date(limits = as.Date(c("2019-01-01", "2024-12-31")),
               date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0))






