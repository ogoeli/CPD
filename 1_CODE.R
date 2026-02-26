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



data_1 <- read.csv('/scratch/ope4/CPD_PAPER/Band_Values_POLY_1.csv')

data_3 <- read.csv('/scratch/ope4/CPD_PAPER/Band_Values_POLY_3.csv')

str(data)

colnames(data_1)[colnames(data_1) == 'classification'] <- 'cover'

colnames(data_3)[colnames(data_3) == 'classification'] <- 'cover'



poly_1 <- data_1 |>
  filter (cover == 4) #3 for healthy and 4 for dead

poly_3 <- data_3 |>
  filter (cover == 4)

data <- rbind(poly_1, poly_3)

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

boundary <- read_sf("/scratch/ope4/POLY_2/POLY_2/POLY_2.shp")


# Ensure your 'data' dataframe has 'longitude' and 'latitude' columns
data_sf <- st_as_sf(data, coords = c("longitude", "latitude"), crs = st_crs(boundary))

# Assign the CRS of the boundary to data_sf
data_sf <- st_set_crs(data_sf, st_crs(boundary)) #|>
  #dplyr::select(-.geo)

# Reproject data_sf to match boundary's CRS
data_sf <- st_transform(data_sf, st_crs(boundary))



#apply tree mask
tree_mask <- brick("/scratch/ope4/MERGE/TREE_MASK/TREE_MASK.TIF")



# Check CRS of points and raster
crs_points <- st_crs(data_sf)  # or crs(points) if using sp package
crs_raster <- crs(tree_mask)


if (!identical(crs_points, crs_raster)) {
  # Transform points to match raster CRS if they differ
  points <- st_transform(data_sf, crs_raster)  # or spTransform for sp objects
}


# Use nearest neighbor interpolation to extract values
extracted_values <- raster::extract(tree_mask, points)  # or 'simple' if preferred

nrow(extracted_values)
nrow(data_sf)

# Ensure extracted_values is a data frame
extracted_values_2 <- as.data.frame(extracted_values)


extracted_values_15 <- cbind(extracted_values_2, data_sf)

table(extracted_values_15$TREE_MASK)
table(extracted_values_15$cover_code)

# Now, convert the dataframe to an sf object and assign CRS, 
##Do not filter for monica points
boundary <- extracted_values_15 |>
  dplyr::filter(TREE_MASK >= 1) ##height greater than or equal to 5

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

# Group by year_month and cover, and calculate the mean for each band
grouped_data <- data %>%
  group_by(year_month, cover) %>%
  summarise(across(starts_with('B'), mean, na.rm = TRUE)) #%>%
  #summarise(across(matches('^(N|CC)'), mean, na.rm = TRUE)) 

grouped_data$year_month <- as.Date(paste0(grouped_data$year_month, "-01"))

# Reshape the data for plotting
long_data <- grouped_data %>%
  pivot_longer(cols = starts_with('B'), names_to = 'Band', values_to = 'Mean_Value') 



#------------------------------------ INDICES

# Group by year_month and cover, and calculate the mean for each band
grouped_data <- data %>%
  group_by(year_month, cover) %>%
  #summarise(across(everything(), mean, na.rm = TRUE)) |>
  summarise(across(matches('^(N|CC|E)'), mean, na.rm = TRUE)) 

##select year
grouped_data$year_month <- as.Date(paste0(grouped_data$year_month, "-01"))


# Reshape the data to long format focusing only on the indices
indices_data <- grouped_data %>%
  pivot_longer(cols = c(NDWI, NDVI, CCCI),  #CIG,
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
    Value = ifelse(is.infinite(Value) | is.nan(Value), NA, Value)  # convert -Inf/NaN to NA
  ) %>%
  mutate(Value = na.approx(Value, na.rm = FALSE)) %>%  # interpolate again
  ungroup()



##we decomposed the original dataset
decomposed_results <- indices_data_filled %>%
  #filter(cover == 4) %>%
  split(.$Index) %>%
  map(~ {
    x <- .x %>% arrange(year_month)
    ts_data <- ts(x$Value, start = c(year(x$year_month[1]), month(x$year_month[1])), frequency = 12)
    decompose(ts_data, type = "additive")
  })

decomposed_results_1 <- decomposed_results$CCCI

# Convert the time series trend to a data frame
df_trend <- data.frame(
  Time = as.numeric(time(decomposed_results_1$trend)),
  Trend = as.numeric(decomposed_results_1$trend)
)

# ggplot version
ggplot(df_trend, aes(x = Time, y = Trend)) +
  geom_line(color = "black", size = 1.2) +
  labs(
    title = "Trend Component of CCCI",
    x = "Date",
    y = "Trend Value"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face="bold", color = "black")
    ) 


##Monica
# Convert time to Date
df_trend$Date <- as.Date(as.yearmon(time(decomposed_results_1$trend)))

# Attack dates
first_date <- as.Date("2023-06-17")
last_date  <- as.Date("2024-08-15")

# ggplot
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
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")


#Table 2. Ranges of vegetation index values for each canopy condition class---------------------------------------
##to get min and max value for each month
summary(na.omit(grouped_data))

###for monica points: 
##to get min and max value for each month
summary(grouped_data)




#### Explan the variation in CPD------------------PRECIP------------------------------------------
# 1. Read your CSV (replace with your actual file path)
precip <- read_csv("/scratch/ope4/CPD_PAPER/Monica_precipitation_time_series.csv")

# 2. Convert date column to Date format
precip_2 <- precip %>%
  filter(precipitation > 0) |>
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

# Aggregate monthly averages
monthly_avg <- precip_2 %>%
  mutate(year_month = floor_date(date, "month")) %>%
  group_by(year_month) %>%
  summarise(avg_precip = mean(precipitation, na.rm = TRUE)) %>%
  ungroup()

# Convert to time series
precip_ts <- ts(monthly_avg$avg_precip,
                start = c(year(min(monthly_avg$year_month)),
                          month(min(monthly_avg$year_month))),
                frequency = 12)  # monthly data

##decompse the trend
precip_decomp <- decompose(precip_ts, type = "additive")  # or type = "multiplicative"

# Plot the decomposition
plot(precip_decomp)

# 3. Plot time series
trend <- precip_decomp$trend

# Assuming you have a data frame with your monthly averages and trend
df_precip <- data.frame(
  Month = monthly_avg$year_month,  # Date or numeric
  Trend = trend                    # Trend values
)

# Attack dates
first_date <- as.Date("2023-06-17")
last_date  <- as.Date("2024-08-15")


# ggplot
ggplot(df_precip, aes(x = Month, y = Trend)) +
  geom_line(color = "black", size = 1) +
  labs(
    title = "Trend Component of Monthly Precipitation",
    x = "Time (Years)",
    y = "Trend in Precipitation (mm)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(face = "bold", color = "black")
  ) +
  geom_vline(xintercept = first_date, linetype = "dashed", color = "red") +
  geom_vline(xintercept = last_date,  linetype = "dashed", color = "red") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  
#
#
  
