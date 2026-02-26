library(sf)
library(sp)
library(randomForest)
library(GWmodel)
library(terra)
library(ranger)
library(caret)
library(blockCV)
library(tidyverse)
library(cowplot)
library(automap)
library(spdep)


data <- st_read('/scratch/ope4/MERGE/USDA/PIPOsplatNewTrees24.gpx') 


data <- dplyr::select(data, ele, time, name, geometry, sym)

data$cover <- 1

st_write(data, "/scratch/ope4/MERGE/USDA/SHP/monica.shp")


plot(data$ele)

data_1 <- read_csv('/scratch/ope4/MERGE/USDA/For Ogonna.csv')

data_1 <- data_1 %>% dplyr::rename(name = `Tree  #`)

data_1 <- data_1 %>% dplyr::rename(attack_period = `When attacks began`)

data_1 <- data_1 %>% dplyr::rename(attack_year = Year)


data_1 <- data_1 %>% mutate(name = as.character(name))

class(data)
class(data_1)

new <- data_1 %>% 
  inner_join(data, by = join_by(name))



head(new)

#new <- dplyr::select(new, attack_period, attack_year, name)

library(sf)

library(sf)

new <- dplyr::select(new, attack_period, attack_year, name, geometry)

# Convert to sf object (if not already)
damage_points_sf <- new %>%
  st_as_sf() %>%  # Convert to sf object
  mutate(
    longitude = st_coordinates(.)[, 1],  # Extract X (longitude)
    latitude = st_coordinates(.)[, 2]   # Extract Y (latitude)
  )


class(damage_points_sf)


damage_points_sf <- damage_points_sf |>
  select(longitude, latitude, `When attacks began`, Year, geometry )




write_csv(new, "/scratch/ope4/MERGE/USDA/damage_points_sf.csv", append=FALSE)

