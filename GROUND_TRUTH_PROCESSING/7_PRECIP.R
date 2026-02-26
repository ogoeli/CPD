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

