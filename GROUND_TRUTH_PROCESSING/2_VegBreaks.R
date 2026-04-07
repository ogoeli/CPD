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


#load in the ground truth dataset with tree mask applied, index calculated
data <- data

#group the dataset by pixel location and time
df <- data |>
  group_by(.geo, year_month_1) 


#---------------------------------------------------------------CPD algorithm
#loop through each pixel time series and detect change point by 3 different algorithms


#BEAST Algorithm--------------------
#A loop to go through each ground truth pixel time series per index
bands <- c("NDVI", "NDWI", "CCCI")

run_all_algorithms <- function(df, band_name) {
  results <- list()
  
  for (geom in unique(df$.geo)) {
    df_2 <- df[df$.geo == geom, ]
    #ts_data <- ts(df_2[[band_name]], start = c(2019, 1), end = c(2024, 12), frequency = 12)
    #ts_filled <- na.approx(ts_data)
    df_2$timestamp <- as.Date(df_2$timestamp)  # adjust if it's named differently
    # Extract year-month and aggregate
    df_2$ym <- format(df_2$timestamp, "%Y-%m")
    monthly_df <- aggregate(df_2[[band_name]], by = list(df_2$ym), FUN = mean)
    names(monthly_df) <- c("ym", band_name)
    # Sort by year-month
    monthly_df <- monthly_df[order(monthly_df$ym), ]
    # Create time series from cleaned data
    time_series <- ts(monthly_df[[band_name]], start=c(2019, 1), frequency=12)
    
    ts_filled <- na.approx(time_series, na.rm = FALSE)
    
    
    
    ### --- BEAST ---
    beast_out <- tryCatch(beast(ts_filled, tcp.minmax = c(0, 3), season = 'harmonic'), error = function(e) NULL)
    
    beast_result <- if (!is.null(beast_out)) {
      change_points <- beast_out$trend$cp
      cp_probabilities <- beast_out$trend$cpPr
      sorted_indices <- order(change_points)
      cp_sorted <- change_points[sorted_indices]
      pr_sorted <- cp_probabilities[sorted_indices]
      
      list(
        method = "BEAST",
        geom = geom,
        band = band_name,
        RMSE = beast_out$RMSE,
        R2 = beast_out$R2,
        ncp = beast_out$trend$ncp,
        first_cp = ifelse(length(cp_sorted) >= 1, cp_sorted[1], NA),
        second_cp = ifelse(length(cp_sorted) >= 2, cp_sorted[2], NA),
        third_cp = ifelse(length(cp_sorted) >= 3, cp_sorted[3], NA),
        first_pr = ifelse(length(pr_sorted) >= 1, pr_sorted[1], NA),
        second_pr = ifelse(length(pr_sorted) >= 2, pr_sorted[2], NA),
        third_pr = ifelse(length(pr_sorted) >= 3, pr_sorted[3], NA)
      )
    } else NULL
    
    
    # Append all available results
    for (res in list(beast_result)) {
      if (!is.null(res)) {
        results[[length(results) + 1]] <- as.data.frame(res)
      }
    }
  }
  
  return(do.call(rbind, results))
}




all_results <- lapply(bands, function(b) run_all_algorithms(df, b))
final_results_BEAST <- do.call(rbind, all_results)


#-------------------------------------------------------------------------------------------------------------------

#

#BFAST Algorithm------------

results_list <- list()
bands <- c("NDVI", "NDWI", "CCCI")

# Build a complete monthly ts from 2019-01 to 2024-12
make_monthly_ts <- function(df_band, band, start_ym = "2019-01", end_ym = "2024-12") {
  df_band$timestamp <- as.Date(df_band$timestamp)
  df_band$ym <- as.yearmon(df_band$timestamp)
  monthly <- aggregate(df_band[[band]], by = list(df_band$ym), FUN = mean, na.rm = TRUE)
  names(monthly) <- c("ym", band)
  monthly <- monthly[order(monthly$ym), ]
  
  ym_seq <- seq(as.yearmon(start_ym), as.yearmon(end_ym), by = 1/12)
  full <- data.frame(ym = ym_seq)
  out <- merge(full, monthly, by = "ym", all.x = TRUE)
  
  start_year  <- floor(as.numeric(ym_seq[1]))
  start_month <- round((as.numeric(ym_seq[1]) - start_year) * 12) + 1
  ts(out[[band]], start = c(start_year, start_month), frequency = 12)
}

# Fractional ts time -> Date (assumes monthly freq=12)
ts_time_to_date <- function(frac_time, freq = 12) {
  year  <- floor(frac_time)
  month <- round((frac_time - year) * freq) + 1
  as.Date(sprintf("%04d-%02d-01", year, month))
}

for (geom in unique(df$.geo)) {
  for (band in bands) {
    
    df_2 <- df[df$.geo == geom, ]
    
    ts_y <- make_monthly_ts(df_2, band, start_ym = "2019-01", end_ym = "2024-12")
    ts_y_imp <- na.approx(ts_y, na.rm = FALSE)
    
    # Skip degenerate series
    if (all(is.na(ts_y_imp)) || length(unique(na.omit(as.numeric(ts_y_imp)))) <= 1) {
      results_list[[paste(geom, band, sep = "_")]] <- data.frame(
        geom, band, RMSE = NA_real_, R2 = NA_real_,
        break1 = as.Date(NA), break2 = as.Date(NA), break3 = as.Date(NA)
      )
      next
    }
    
    out <- tryCatch(
      bfast(ts_y_imp, breaks = 3, season = "harmonic"),
      error = function(e) NULL
    )
    
    if (is.null(out) || is.null(out$output) || length(out$output) == 0) {
      results_list[[paste(geom, band, sep = "_")]] <- data.frame(
        geom, band, RMSE = NA_real_, R2 = NA_real_,
        break1 = as.Date(NA), break2 = as.Date(NA), break3 = as.Date(NA)
      )
      next
    }
    
    # Final iteration
    k <- length(out$output)
    Yt <- out$Yt
    Tt <- out$output[[k]]$Tt
    St <- out$output[[k]]$St
    fit <- Tt + St
    
    # Metrics on common finite indices
    Y <- as.numeric(Yt); F <- as.numeric(fit)
    ok <- is.finite(Y) & is.finite(F)
    if (!any(ok)) {
      rmse <- NA_real_; r2 <- NA_real_
    } else {
      Y_ok <- Y[ok]; F_ok <- F[ok]
      rmse <- sqrt(mean((Y_ok - F_ok)^2))
      sse  <- sum((Y_ok - F_ok)^2)
      sst  <- sum((Y_ok - mean(Y_ok))^2)
      r2   <- if (sst > 0) 1 - sse/sst else NA_real_
    }
    
    # Breakpoints -> Dates (use Vt time, drop NAs first)
    bpV <- out$output[[k]]$bp.Vt
    if (is.null(bpV) || !inherits(bpV, "breakpoints") || length(bpV$breakpoints) == 0) {
      br_dates <- rep(as.Date(NA), 3)
    } else {
      Vt <- out$output[[k]]$Vt
      Vn <- as.numeric(Vt)
      ok_v <- !is.na(Vn)
      t_ok <- time(Vt)[ok_v]           # fractional time values
      idx  <- bpV$breakpoints
      idx  <- idx[idx >= 1 & idx <= length(t_ok)]
      brt  <- t_ok[idx]
      br_dates <- ts_time_to_date(brt)
      if (length(br_dates) < 3) br_dates <- c(br_dates, rep(as.Date(NA), 3 - length(br_dates)))
      if (length(br_dates) > 3) br_dates <- br_dates[1:3]
    }
    
    results_list[[paste(geom, band, sep = "_")]] <- data.frame(
      geom = geom, band = band, RMSE = rmse, R2 = r2,
      break1 = br_dates[1], break2 = br_dates[2], break3 = br_dates[3]
    )
  }
}

final_results_BFAST <- do.call(rbind, results_list) %>%
  distinct(geom, band, .keep_all = TRUE)

head(final_results_BFAST)

#-------------------------------------------------------------------------------


#DBEST Algorithm-----------------
# Build a complete monthly ts from 2019-01 to 2024-12
make_monthly_ts <- function(df_band, band, start_ym="2019-01", end_ym="2024-12") {
  df_band$timestamp <- as.Date(df_band$timestamp)
  df_band$ym <- as.yearmon(df_band$timestamp)
  monthly <- aggregate(df_band[[band]], by=list(df_band$ym), mean, na.rm=TRUE)
  names(monthly) <- c("ym", band)
  monthly <- monthly[order(monthly$ym), ]
  ym_seq <- seq(as.yearmon(start_ym), as.yearmon(end_ym), by = 1/12)
  full <- data.frame(ym = ym_seq)
  out <- merge(full, monthly, by="ym", all.x=TRUE)
  
  start_year  <- floor(as.numeric(ym_seq[1]))
  start_month <- round((as.numeric(ym_seq[1]) - start_year) * 12) + 1
  ts(out[[band]], start=c(start_year, start_month), frequency=12)
}

idx_to_date <- function(idx, start_year=2019, start_month=1) {
  if (is.na(idx)) return(NA)
  year  <- start_year + ((start_month - 1 + idx - 1) %/% 12)
  month <- ((start_month - 1 + idx - 1) %% 12) + 1
  as.Date(sprintf("%04d-%02d-01", year, month))
}

results_list_DBEST <- list()
bands <- c("NDVI", "NDWI", "CCCI")

for (geom in unique(df$.geo)) {
  for (band in bands) {
    df_2 <- df[df$.geo == geom, ]
    
    ts_y <- make_monthly_ts(df_2, band, start_ym="2019-01", end_ym="2024-12")
    ts_filled <- na.approx(ts_y, na.rm=FALSE)
    
    out <- DBEST(ts_filled, data.type="cyclical", algorithm="change detection",
                 breakpoints.no=3, first.level.shift=0.1, second.level.shift=0.2,
                 duration=148, distance.threshold="default", alpha=0.05, plot="off")
    
    # metrics on common finite indices just in case
    ok <- is.finite(out$Data) & is.finite(out$Fit)
    rmse <- if (any(ok)) sqrt(mean((out$Data[ok] - out$Fit[ok])^2)) else NA_real_
    r2   <- if (any(ok)) {
      sse <- sum((out$Data[ok] - out$Fit[ok])^2)
      sst <- sum((out$Data[ok] - mean(out$Data[ok]))^2)
      if (sst > 0) 1 - sse/sst else NA_real_
    } else NA_real_
    
    cp_idx <- out$Start
    sig    <- out$Significance
    length(cp_idx) <- 3; length(sig) <- 3
    valid <- !is.na(cp_idx)
    if (sum(valid) > 1) {
      ord <- order(cp_idx[valid])
      cp_idx[valid] <- cp_idx[valid][ord]
      sig[valid]    <- sig[valid][ord]
    }
    
    results_list_DBEST[[paste(geom, band, sep="_")]] <- data.frame(
      geom = geom, band = band, RMSE = rmse, R2 = r2,
      ncp = sum(!is.na(cp_idx)),
      cp_idx_1 = cp_idx[1], cp_idx_2 = cp_idx[2], cp_idx_3 = cp_idx[3],
      cp_date_1 = idx_to_date(cp_idx[1]),
      cp_date_2 = idx_to_date(cp_idx[2]),
      cp_date_3 = idx_to_date(cp_idx[3]),
      cp_sig_1 = sig[1], cp_sig_2 = sig[2], cp_sig_3 = sig[3]
    )
  }
}

final_results_DBEST <- do.call(rbind, results_list_DBEST)

summary(final_results_DBEST)

#------------------------------------------------------------------------------------

