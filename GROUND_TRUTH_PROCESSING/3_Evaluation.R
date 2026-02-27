#Rename the column .geo in df from step 2 to geom
#so it can be used consistently as a join key.
colnames(df)[colnames(df) == '.geo'] <- 'geom'
df <- df[!duplicated(df[c("geom")]), ]

#Removes duplicate rows in df based on the geom column, 
#keeping the first occurrence of each unique geometry.
inner_join_BEAST <- inner_join(final_results_BEAST, df, by = "geom")
inner_join_BFAST <- inner_join(final_results_BFAST, df, by = "geom")
inner_join_DBEST <- inner_join(final_results_DBEST, df, by = "geom")

#import actual infected points, similair to the points from the Sentinel 2 SR time series
#same points from step 1 and 2
USDA <- read.csv('/scratch/ope4/CPD_PAPER/damage_points_sf.csv')

USDA <- USDA %>%
  mutate(name = as.character(name))

##join the result from step 2 to the actual ground truth from field visit
#this retains the period between last and next visit when a tree was marked as infected
dataset_BEAST <- inner_join_BEAST %>% 
  inner_join(USDA, by = join_by(name)) |>
  dplyr::select(band, RMSE, R2, first_cp, second_cp, third_cp, period_1, period_2)

dataset_BFAST <- inner_join_BFAST %>% 
  inner_join(USDA, by = join_by(name)) |>
  dplyr::select(band, RMSE, R2, break1, break2, break3, period_1, period_2)

dataset_DBEST <- inner_join_DBEST %>% 
  inner_join(USDA, by = join_by(name)) |>
  dplyr::select(band, RMSE, R2, cp_date_1, cp_date_2, cp_date_3, period_1, period_2)

##let all dataset have same col names as date1..
#where date1,2,3 signifys the date of 3 changepoint detected 
# Standardize DBEST names
dataset_DBEST <- dataset_DBEST %>%
  rename(
    date1 = cp_date_1,
    date2 = cp_date_2,
    date3 = cp_date_3
  )
# Standardize BEAST names
dataset_BEAST <- dataset_BEAST %>%
  rename(
    date1 = first_cp,
    date2 = second_cp,
    date3 = third_cp
  )
# Standardize BFAST names
dataset_BFAST <- dataset_BFAST %>%
  rename(
    date1 = break1,
    date2 = break2,
    date3 = break3
  )

##GET MEAN DIFF IN DAYS
#Table 4. Variable performance. ----------------------------------------

#function to convert decimal years to actual date
#some algorithms outputs CP in decimal
decimal_year_to_date <- function(decimal_year) {
  year <- floor(decimal_year)
  fraction <- decimal_year - year
  days_in_year <- 365
  days_to_add <- round(fraction * days_in_year)
  as.Date(paste0(year, "-01-01")) + days_to_add
}
#change the date str to make Date format in R
dataset_BFAST <- dataset_BFAST %>%
  mutate(
    date1  = as.numeric(date1),
    date2 = as.numeric(date2),
    date3  = as.numeric(date3)
  )


##DOES NOT APPLY TO BFAST AND DBEST
#only BEAST outputs their date in decimal
dataset_BFAST <- dataset_BFAST %>%
  mutate(
    date1  = decimal_year_to_date(date1),
    date2 = decimal_year_to_date(date2),
    date3  = decimal_year_to_date(date3)
  )
# Ensure dates are in Date format
dataset_BFAST <- dataset_BFAST %>%
  mutate(
    period_1   = as.Date(period_1, format = "%m/%d/%Y"),
    period_2   = as.Date(period_2, format = "%m/%d/%Y"),
    date1  = as.Date(date1),
    date2  = as.Date(date2),
    date3  = as.Date(date3)
  )

# Compute the date of the actual period
#used the last date when a tree was marked as infected
dataset_BFAST <- dataset_BFAST %>%
  mutate(actual_mid = period_2)
# Convert to Date if needed
dataset_BFAST <- dataset_BFAST %>%
  mutate(across(c(date1, date2, date3, period_1, period_2), as.Date))

# Calculate absolute differences for each predicted-actual pair
#using the last date with observed infections
dataset_diff <- dataset_BFAST %>%
  rowwise() %>%
  mutate(
    diff_1_2 = abs(as.numeric(date1 - period_2)),
    diff_2_2 = abs(as.numeric(date2 - period_2)),
    diff_3_2 = abs(as.numeric(date3 - period_2))
  ) %>%
  ungroup()

##get the mean difference 
#average all 11 points for each algorithm across 3 index
performance_dates_summary <- dataset_diff %>%
  group_by(band) %>%
  summarise(
    mean_diff_cp1 = mean(diff_1_2, na.rm = TRUE),
    mean_diff_cp2 = mean(diff_2_2, na.rm = TRUE),
    mean_diff_cp3 = mean(diff_3_2, na.rm = TRUE),
    
    median_diff_cp1 = median(diff_1_2, na.rm = TRUE),
    median_diff_cp2 = median(diff_2_2, na.rm = TRUE),
    median_diff_cp3 = median(diff_3_2, na.rm = TRUE),
    
    n = n()
  ) %>%
  arrange(mean_diff_cp1)




##TABLE 3: MEAN R2 AND RMSE---------------------------------
dataset_DBEST %>%
  group_by(band) %>%
  summarise(
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    mean_R2 = mean(R2, na.rm = TRUE)
  )


#Table 2. Ranges of vegetation index values for each canopy condition class---------------------------------------
##to get min and max value for each month
summary(grouped_data)
