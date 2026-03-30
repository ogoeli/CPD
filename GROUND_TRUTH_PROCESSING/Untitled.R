library(Rbeast)
library(zoo)

run_beast_and_save_plots <- function(df, band_name, output_dir = "BEAST_Plots") {
  
  # Create main output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Create band-specific folder
  band_dir <- file.path(output_dir, band_name)
  if (!dir.exists(band_dir)) {
    dir.create(band_dir)
  }
  
  for (geom in unique(df$.geo)) {
    
    df_2 <- df[df$.geo == geom, ]
    
    # Convert timestamp
    df_2$timestamp <- as.Date(df_2$timestamp)
    df_2$ym <- format(df_2$timestamp, "%Y-%m")
    
    # Monthly aggregation
    monthly_df <- aggregate(df_2[[band_name]], 
                            by = list(df_2$ym), 
                            FUN = mean)
    
    names(monthly_df) <- c("ym", band_name)
    monthly_df <- monthly_df[order(monthly_df$ym), ]
    
    if (nrow(monthly_df) < 12) next  # skip short time series
    
    # Extract correct start year/month
    start_year  <- as.numeric(substr(monthly_df$ym[1], 1, 4))
    start_month <- as.numeric(substr(monthly_df$ym[1], 6, 7))
    
    # Create time series
    ts_data <- ts(monthly_df[[band_name]],
                  start = c(start_year, start_month),
                  frequency = 12)
    
    ts_filled <- na.approx(ts_data, na.rm = FALSE)
    
    if (all(is.na(ts_filled))) next
    
    # Run BEAST safely
    beast_out <- tryCatch(
      beast(ts_filled, tcp.minmax = c(0, 3), season = 'harmonic'),
      error = function(e) NULL
    )
    
    if (!is.null(beast_out)) {
      
      # Save plot
      file_name <- paste0("BEAST_", band_name, "_", gsub("[^0-9]", "", geom), ".png")
      file_path <- file.path(band_dir, file_name)
      
      png(file_path, width = 1200, height = 800, res = 150)
      plot(beast_out)
      dev.off()
    }
  }
}


##plot all index for all pixel
bands <- c("NDVI", "NDWI", "CCCI")

for (b in bands) {
  run_beast_and_save_plots(df, b)
}
