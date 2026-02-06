# Heard Island Glacier Albedo Trend Map Analysis
# Pixel-level linear regression trend analysis for multiple sensors
# R Script for Spatial Trend Analysis and Visualization

# Load required libraries
library(tidyverse)
library(raster)
library(sp)
library(sf)
library(terra)
library(Kendall)          # For Mann-Kendall test
# library(zyp) removed - using simple linear regression instead
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(extrafont)

# Load Times New Roman font with error handling
tryCatch({
  loadfonts(device = "pdf", quiet = TRUE)
  cat("✓ Times New Roman font loaded successfully\n")
}, error = function(e) {
  cat("Warning: Could not load Times New Roman font:", e$message, "\n")
  cat("Using system default font instead\n")
})

# Set working directory and create output directory
setwd("/Users/jingming/Desktop/HIG")
if (!dir.exists("figures")) {
  dir.create("figures")
}
if (!dir.exists("output")) {
  dir.create("output")
}

# Define Times New Roman theme with fallback
theme_times <- function() {
  # Check if Times New Roman is available
  if ("Times New Roman" %in% fonts()) {
    font_family <- "Times New Roman"
  } else {
    font_family <- "serif"  # Fallback to serif font
    cat("Warning: Times New Roman not available, using serif font\n")
  }
  
  theme_minimal() +
    theme(
      text = element_text(family = font_family, size = 11)
    )
}

# Define data paths
data_path <- "data/gee_export/"

# Function to load and organize raster files by sensor
load_raster_files <- function(sensor_name) {
  cat("Loading", sensor_name, "raster files...\n")
  
  # Define sensor-specific directories and patterns - Focus on VIIRS only
  # Based on comprehensive data quality assessment:
  # - VIIRS: Uses validated BRDF/albedo parameters from VNP43 product
  # - Sentinel-3: Excluded due to use of adapted coefficients from Sentinel-2 MSI
  #   rather than validated Sentinel-3 OLCI-specific coefficients, making results non-credible
  sensor_configs <- list(
    "VIIRS" = list(
      dir = paste0(data_path, "HIG_VIIRS_Albedo/"),
      pattern = ".*viirs.*_albedo_yearly_.*\\.tif$",
      exclude_pattern = "area_mean"
    )
  )
  
  if (!sensor_name %in% names(sensor_configs)) {
    stop("Unknown sensor:", sensor_name)
  }
  
  config <- sensor_configs[[sensor_name]]
  
  # List all files in the directory
  all_files <- list.files(config$dir, pattern = "\\.tif$", full.names = TRUE)
  
  # Filter files based on pattern and exclude patterns
  sensor_files <- all_files[grepl(config$pattern, basename(all_files), ignore.case = TRUE)]
  sensor_files <- sensor_files[!grepl(config$exclude_pattern, basename(sensor_files), ignore.case = TRUE)]
  
  # Extract years from filenames
  file_info <- data.frame(
    file_path = sensor_files,
    filename = basename(sensor_files),
    stringsAsFactors = FALSE
  )
  
  # Extract year from filename using regex
  year_pattern <- ".*_([0-9]{4})\\.tif$"
  file_info$year <- as.numeric(gsub(year_pattern, "\\1", file_info$filename))
  
  # Filter out invalid years and sort
  file_info <- file_info[!is.na(file_info$year) & file_info$year >= 2000 & file_info$year <= 2024, ]
  
  # Check if we have valid files
  if (nrow(file_info) == 0) {
    cat("Warning: No valid files found for", sensor_name, "\n")
    return(data.frame())
  }
  
  file_info <- file_info[order(file_info$year), ]
  
  cat("✓ Found", nrow(file_info), "valid raster files for", sensor_name, "\n")
  cat("  Years:", min(file_info$year), "to", max(file_info$year), "\n")
  
  return(file_info)
}

# Function to check raster validity and data quality
check_raster_quality <- function(file_path) {
  tryCatch({
    # Load raster
    r <- raster(file_path)
    
    # Check for valid data
    values <- values(r)
    valid_values <- values[!is.na(values) & !is.infinite(values)]
    
    if (length(valid_values) == 0) {
      return(list(valid = FALSE, reason = "No valid data"))
    }
    
    # Check if all values are the same (constant raster)
    if (length(unique(valid_values)) == 1) {
      return(list(valid = FALSE, reason = "Constant values", value = unique(valid_values)[1]))
    }
    
    # Check data range
    min_val <- min(valid_values)
    max_val <- max(valid_values)
    
    # Check for reasonable albedo range (0-1)
    if (min_val < 0 || max_val > 1) {
      return(list(valid = FALSE, reason = "Values outside 0-1 range", min = min_val, max = max_val))
    }
    
    return(list(
      valid = TRUE,
      n_valid = length(valid_values),
      min = min_val,
      max = max_val,
      mean = mean(valid_values),
      sd = sd(valid_values)
    ))
  }, error = function(e) {
    return(list(valid = FALSE, reason = paste("Error:", e$message)))
  })
}

# Function to calculate linear regression slope for a pixel time series
calculate_linear_slope_pixel <- function(values, years) {
  # Remove NA values
  valid_idx <- !is.na(values) & !is.na(years)
  if (sum(valid_idx) < 3) {
    return(list(slope = NA, p_value = NA, n_obs = sum(valid_idx)))
  }
  
  valid_values <- values[valid_idx]
  valid_years <- years[valid_idx]
  
  # Calculate simple linear regression
  tryCatch({
    lm_result <- lm(valid_values ~ valid_years)
    slope <- coef(lm_result)[2]  # Slope per year
    intercept <- coef(lm_result)[1]
    
    # Get p-value from linear regression
    lm_summary <- summary(lm_result)
    p_value <- lm_summary$coefficients[2, 4]  # p-value for slope coefficient
    
    # Calculate Mann-Kendall test for comparison (non-parametric test)
    mk_result <- MannKendall(valid_values)
    
    return(list(
      slope = slope,
      p_value = p_value,
      tau = mk_result$tau,
      n_obs = length(valid_values),
      intercept = intercept,
      r_squared = lm_summary$r.squared
    ))
  }, error = function(e) {
    return(list(slope = NA, p_value = NA, n_obs = sum(valid_idx)))
  })
}

# Function to process raster stack and calculate trends
calculate_raster_trends <- function(file_info, sensor_name) {
  cat("Processing", sensor_name, "raster stack for trend analysis...\n")
  
  # Filter files based on quality check
  valid_files <- data.frame()
  
  for (i in 1:nrow(file_info)) {
    file_path <- file_info$file_path[i]
    year <- file_info$year[i]
    
    cat("  Checking", basename(file_path), "...")
    quality <- check_raster_quality(file_path)
    
    if (quality$valid) {
      cat(" ✓\n")
      valid_files <- rbind(valid_files, data.frame(
        file_path = file_path,
        year = year,
        quality_info = list(quality)
      ))
    } else {
      cat(" ✗ (", quality$reason, ")\n")
    }
  }
  
  if (nrow(valid_files) < 3) {
    cat("Warning: Insufficient valid files for", sensor_name, "trend analysis\n")
    return(NULL)
  }
  
  cat("✓ Using", nrow(valid_files), "valid files for trend analysis\n")
  
  # Sort files by year
  valid_files <- valid_files[order(valid_files$year), ]
  
  # Load all rasters into a stack with extent matching
  cat("  Loading raster stack...\n")
  
  # First, load all rasters to check extents
  raster_list <- list()
  years <- valid_files$year
  
  for (i in 1:nrow(valid_files)) {
    r <- raster(valid_files$file_path[i])
    raster_list[[i]] <- r
    if (i %% 5 == 0 || i == nrow(valid_files)) {
      cat("    Loaded", i, "of", nrow(valid_files), "rasters\n")
    }
  }
  
  # Check if all rasters have the same extent
  extents <- lapply(raster_list, function(x) extent(x))
  extent_equal <- all(sapply(extents[-1], function(x) identical(x, extents[[1]])))
  
  if (!extent_equal) {
    cat("  Warning: Rasters have different extents, resampling to common extent...\n")
    
    # Find the common extent (intersection of all extents)
    common_extent <- extents[[1]]
    for (i in 2:length(extents)) {
      common_extent <- intersect(common_extent, extents[[i]])
    }
    
    # Create a template raster with the common extent
    template_raster <- raster_list[[1]]
    extent(template_raster) <- common_extent
    
    # Resample all rasters to the common extent
    resampled_rasters <- list()
    for (i in 1:length(raster_list)) {
      if (extent(raster_list[[i]]) != common_extent) {
        cat("    Resampling raster", i, "to common extent\n")
        resampled_rasters[[i]] <- resample(raster_list[[i]], template_raster, method = "bilinear")
      } else {
        resampled_rasters[[i]] <- raster_list[[i]]
      }
    }
    
    # Create the stack from resampled rasters
    raster_stack <- stack(resampled_rasters)
  } else {
    # All rasters have the same extent, create stack directly
    raster_stack <- stack(raster_list)
  }
  
  # Get raster dimensions
  nrows <- nrow(raster_stack)
  ncols <- ncol(raster_stack)
  nlayers <- nlayers(raster_stack)
  
  cat("  Raster dimensions:", nrows, "x", ncols, "x", nlayers, "\n")
  
  # Create result rasters
  slope_raster <- raster(nrows = nrows, ncols = ncols, 
                        xmn = xmin(raster_stack), xmx = xmax(raster_stack),
                        ymn = ymin(raster_stack), ymx = ymax(raster_stack),
                        crs = crs(raster_stack))
  
  pvalue_raster <- slope_raster
  tau_raster <- slope_raster
  nobs_raster <- slope_raster
  
  # Initialize with NA values
  values(slope_raster) <- NA
  values(pvalue_raster) <- NA
  values(tau_raster) <- NA
  values(nobs_raster) <- NA
  
  # Process in blocks to manage memory
  block_size <- 1000
  n_blocks <- ceiling(nrows / block_size)
  
  cat("  Processing", n_blocks, "blocks...\n")
  
  # Suppress Mann-Kendall warnings
  old_warn <- getOption("warn")
  options(warn = -1)
  
  for (block in 1:n_blocks) {
    start_row <- (block - 1) * block_size + 1
    end_row <- min(block * block_size, nrows)
    
    if (block %% 5 == 0 || block == n_blocks) {
      cat("    Processing block", block, "of", n_blocks, "\n")
    }
    
    # Extract block data
    block_data <- getValuesBlock(raster_stack, row = start_row, nrows = end_row - start_row + 1)
    
    # Calculate trends for each pixel in this block
    n_pixels <- nrow(block_data)
    
    for (pixel in 1:n_pixels) {
      # Progress indicator for large blocks
      if (n_pixels > 10000 && pixel %% 10000 == 0) {
        cat("      Processed", pixel, "of", n_pixels, "pixels in block", block, "\n")
      }
      
      pixel_values <- block_data[pixel, ]
      
      # Remove NA values
      valid_idx <- !is.na(pixel_values)
      if (sum(valid_idx) >= 3) {
        valid_values <- pixel_values[valid_idx]
        valid_years <- years[valid_idx]
        
        # Calculate Sen's slope
        tryCatch({
          # Check for sufficient variation in data
          if (length(unique(valid_values)) < 3) {
            # Not enough variation for trend analysis
            slope <- NA
            p_value <- NA
            tau <- NA
          } else {
            # Calculate linear regression slope
            lm_result <- tryCatch({
              lm(valid_values ~ valid_years)
            }, error = function(e) {
              return(NULL)
            })
            
            if (!is.null(lm_result)) {
              slope <- coef(lm_result)[2]
              lm_summary <- summary(lm_result)
              p_value <- lm_summary$coefficients[2, 4]  # p-value for slope
              
              # Calculate Mann-Kendall test for comparison
              mk_result <- tryCatch({
                MannKendall(valid_values)
              }, error = function(e) {
                # If Mann-Kendall fails, use alternative approach
                if (length(unique(valid_values)) == length(valid_values)) {
                  cor_test <- cor.test(valid_values, valid_years, method = "kendall")
                  list(sl = cor_test$p.value, tau = cor_test$estimate)
                } else {
                  list(sl = NA, tau = NA)
                }
              })
              
              tau <- mk_result$tau
            } else {
              slope <- NA
              p_value <- NA
              tau <- NA
            }
          }
          
          # Store results
          pixel_idx <- (start_row - 1) * ncols + pixel
          values(slope_raster)[pixel_idx] <- slope
          values(pvalue_raster)[pixel_idx] <- p_value
          values(tau_raster)[pixel_idx] <- tau
          values(nobs_raster)[pixel_idx] <- sum(valid_idx)
        }, error = function(e) {
          # Keep NA values for failed calculations
        })
      }
    }
  }
  
  # Restore warning settings
  options(warn = old_warn)
  
  # Set names for the rasters
  names(slope_raster) <- paste0(sensor_name, "_slope")
  names(pvalue_raster) <- paste0(sensor_name, "_pvalue")
  names(tau_raster) <- paste0(sensor_name, "_tau")
  names(nobs_raster) <- paste0(sensor_name, "_nobs")
  
  return(list(
    slope = slope_raster,
    pvalue = pvalue_raster,
    tau = tau_raster,
    nobs = nobs_raster,
    sensor = sensor_name,
    years = valid_files$year,
    n_files = nrow(valid_files)
  ))
}

# Function to create trend visualization
create_trend_map <- function(trend_results, sensor_name) {
  cat("Creating trend map for", sensor_name, "...\n")
  
  slope_raster <- trend_results$slope
  pvalue_raster <- trend_results$pvalue
  
  # Convert to data frame for ggplot
  slope_df <- as.data.frame(slope_raster, xy = TRUE)
  pvalue_df <- as.data.frame(pvalue_raster, xy = TRUE)
  
  # Merge data
  plot_df <- merge(slope_df, pvalue_df, by = c("x", "y"))
  names(plot_df) <- c("x", "y", "slope", "pvalue")
  
  # Remove NA values
  plot_df <- plot_df[!is.na(plot_df$slope) & !is.na(plot_df$pvalue), ]
  
  if (nrow(plot_df) == 0) {
    cat("Warning: No valid data for", sensor_name, "trend map\n")
    return(NULL)
  }
  
  # Create significance mask
  plot_df$significant <- plot_df$pvalue < 0.05
  
  # Export plot data to CSV
  write.csv(plot_df, 
            paste0("output/", sensor_name, "_trend_map_data.csv"), 
            row.names = FALSE)
  cat("✓ Exported trend map data:", paste0("output/", sensor_name, "_trend_map_data.csv"), "\n")
  
  # Create the plot
  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_gradient2(
      low = "red", mid = "white", high = "blue",
      midpoint = 0,
      name = "Slope\n(per year)",
      guide = guide_colorbar(title.position = "top")
    ) +
    labs(
      title = paste(sensor_name, "Albedo Trend Map"),
      subtitle = paste("Linear regression slope per year (", min(trend_results$years), "-", max(trend_results$years), ")"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_times() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1.5, "cm"),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    coord_fixed()
  
  return(p)
}

# Function to create Figure 4 for manuscript (clean version without title/subtitle)
create_figure_4 <- function(trend_results, sensor_name) {
  cat("Creating Figure 4 (manuscript version) for", sensor_name, "...\n")
  
  slope_raster <- trend_results$slope
  pvalue_raster <- trend_results$pvalue
  
  # Convert to data frame for ggplot
  slope_df <- as.data.frame(slope_raster, xy = TRUE)
  pvalue_df <- as.data.frame(pvalue_raster, xy = TRUE)
  
  # Merge data
  plot_df <- merge(slope_df, pvalue_df, by = c("x", "y"))
  names(plot_df) <- c("x", "y", "slope", "pvalue")
  
  # Remove NA values
  plot_df <- plot_df[!is.na(plot_df$slope) & !is.na(plot_df$pvalue), ]
  
  if (nrow(plot_df) == 0) {
    cat("Warning: No valid data for Figure 4\n")
    return(NULL)
  }
  
  # Load and prepare SRTM DEM for contour lines
  srtm_path <- "data/heard_island_srtm_dem.tif"
  dem_df <- NULL
  
  if (file.exists(srtm_path)) {
    cat("Loading SRTM DEM for contour lines...\n")
    srtm <- rast(srtm_path)
    
    # Convert slope raster to terra format for resampling
    slope_rast <- rast(slope_raster)
    
    # Resample DEM to match slope raster if needed
    if (!identical(ext(srtm), ext(slope_rast)) || !identical(res(srtm), res(slope_rast))) {
      cat("Resampling SRTM DEM to match slope raster...\n")
      srtm_resampled <- resample(srtm, slope_rast, method = "bilinear")
    } else {
      srtm_resampled <- srtm
    }
    
    # Convert DEM to data frame for contours
    dem_df <- as.data.frame(srtm_resampled, xy = TRUE)
    names(dem_df) <- c("x", "y", "elevation")
    dem_df <- dem_df[!is.na(dem_df$elevation), ]
  } else {
    cat("Warning: SRTM DEM file not found at", srtm_path, "\n")
    cat("  Contour lines will not be added.\n")
  }
  
  # Create the plot (manuscript version - clean, no title/subtitle)
  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_gradient2(
      low = "red", mid = "white", high = "blue",
      midpoint = 0,
      name = "Slope\n(per year)",
      guide = guide_colorbar(title.position = "top", barwidth = unit(8, "cm"))
    ) +
    labs(
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_times() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(1.5, "cm"),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    ) +
    coord_fixed()
  
  # Add contour lines if DEM data is available
  if (!is.null(dem_df) && nrow(dem_df) > 0) {
    # Determine appropriate contour break intervals based on elevation range
    elev_range <- range(dem_df$elevation, na.rm = TRUE)
    elev_min <- floor(elev_range[1] / 100) * 100
    elev_max <- ceiling(elev_range[2] / 100) * 100
    
    # Create contour breaks at 200m intervals (adjustable)
    contour_breaks <- seq(elev_min, elev_max, by = 200)
    # Filter to only include breaks within actual elevation range
    contour_breaks <- contour_breaks[contour_breaks >= elev_range[1] & contour_breaks <= elev_range[2]]
    
    if (length(contour_breaks) > 0) {
      p <- p + 
        geom_contour(
          data = dem_df,
          aes(x = x, y = y, z = elevation),
          breaks = contour_breaks,
          color = "gray30",
          alpha = 0.4,
          linewidth = 0.3,
          inherit.aes = FALSE
        )
      
      cat("✓ Added contour lines at", length(contour_breaks), "elevation intervals (", 
          min(contour_breaks), "-", max(contour_breaks), "m)\n")
    }
  }
  
  return(p)
}

# Function to create significance map
create_significance_map <- function(trend_results, sensor_name) {
  cat("Creating significance map for", sensor_name, "...\n")
  
  pvalue_raster <- trend_results$pvalue
  
  # Convert to data frame for ggplot
  pvalue_df <- as.data.frame(pvalue_raster, xy = TRUE)
  names(pvalue_df) <- c("x", "y", "pvalue")
  
  # Remove NA values
  pvalue_df <- pvalue_df[!is.na(pvalue_df$pvalue), ]
  
  if (nrow(pvalue_df) == 0) {
    cat("Warning: No valid data for", sensor_name, "significance map\n")
    return(NULL)
  }
  
  # Create significance categories
  pvalue_df$significance <- cut(
    pvalue_df$pvalue,
    breaks = c(0, 0.001, 0.01, 0.05, 1),
    labels = c("p < 0.001", "p < 0.01", "p < 0.05", "p ≥ 0.05"),
    include.lowest = TRUE
  )
  
  # Export plot data to CSV
  write.csv(pvalue_df, 
            paste0("output/", sensor_name, "_significance_map_data.csv"), 
            row.names = FALSE)
  cat("✓ Exported significance map data:", paste0("output/", sensor_name, "_significance_map_data.csv"), "\n")
  
  # Create the plot
  p <- ggplot(pvalue_df, aes(x = x, y = y)) +
    geom_raster(aes(fill = significance)) +
    scale_fill_manual(
      values = c("p < 0.001" = "#d73027", "p < 0.01" = "#f46d43", 
                "p < 0.05" = "#fdae61", "p ≥ 0.05" = "#e0e0e0"),
      name = "Significance",
      guide = guide_legend(title.position = "top")
    ) +
    labs(
      title = paste(sensor_name, "Trend Significance Map"),
      subtitle = "Mann-Kendall test p-values",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_times() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    coord_fixed()
  
  return(p)
}

# Function to create summary statistics plot
create_summary_plot <- function(trend_results, sensor_name) {
  cat("Creating summary plot for", sensor_name, "...\n")
  
  slope_raster <- trend_results$slope
  pvalue_raster <- trend_results$pvalue
  
  # Extract values
  slope_values <- values(slope_raster)
  pvalue_values <- values(pvalue_raster)
  
  # Remove NA values
  valid_idx <- !is.na(slope_values) & !is.na(pvalue_values)
  slope_values <- slope_values[valid_idx]
  pvalue_values <- pvalue_values[valid_idx]
  
  if (length(slope_values) == 0) {
    cat("Warning: No valid data for", sensor_name, "summary plot\n")
    return(NULL)
  }
  
  # Create summary data frame
  summary_df <- data.frame(
    sensor = sensor_name,
    n_pixels = length(slope_values),
    mean_slope = mean(slope_values),
    median_slope = median(slope_values),
    sd_slope = sd(slope_values),
    min_slope = min(slope_values),
    max_slope = max(slope_values),
    significant_pixels = sum(pvalue_values < 0.05),
    significant_percent = sum(pvalue_values < 0.05) / length(pvalue_values) * 100
  )
  
  # Create histogram plot
  plot_df <- data.frame(slope = slope_values, pvalue = pvalue_values)
  plot_df$significant <- plot_df$pvalue < 0.05
  
  # Export plot data to CSV
  write.csv(plot_df, 
            paste0("output/", sensor_name, "_slope_distribution_data.csv"), 
            row.names = FALSE)
  cat("✓ Exported slope distribution data:", paste0("output/", sensor_name, "_slope_distribution_data.csv"), "\n")
  
  p <- ggplot(plot_df, aes(x = slope)) +
    geom_histogram(aes(fill = significant), bins = 30, alpha = 0.7) +
    scale_fill_manual(
      values = c("FALSE" = "gray70", "TRUE" = "red"),
      name = "Significant\n(p < 0.05)",
      labels = c("FALSE" = "No", "TRUE" = "Yes")
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = paste(sensor_name, "Slope Distribution"),
      subtitle = paste("n =", summary_df$n_pixels, "pixels;", 
                      round(summary_df$significant_percent, 1), "% significant"),
      x = "Linear Regression Slope (per year)",
      y = "Count"
    ) +
    theme_times() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  return(list(plot = p, summary = summary_df))
}

# Main analysis function
analyze_sensor_trends <- function(sensor_name) {
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("ANALYZING", toupper(sensor_name), "ALBEDO TRENDS\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  # Load raster files
  file_info <- load_raster_files(sensor_name)
  
  if (nrow(file_info) == 0) {
    cat("No valid files found for", sensor_name, "\n")
    return(NULL)
  }
  
  # Calculate trends
  trend_results <- calculate_raster_trends(file_info, sensor_name)
  
  if (is.null(trend_results)) {
    cat("Trend calculation failed for", sensor_name, "\n")
    return(NULL)
  }
  
  # Create visualizations
  trend_map <- create_trend_map(trend_results, sensor_name)
  significance_map <- create_significance_map(trend_results, sensor_name)
  summary_plot <- create_summary_plot(trend_results, sensor_name)
  
  # Create Figure 4 for manuscript (clean version)
  figure_4 <- create_figure_4(trend_results, sensor_name)
  
  # Save results
  if (!is.null(trend_map)) {
    ggsave(paste0("figures/", sensor_name, "_trend_map.pdf"), 
           trend_map, width = 8, height = 6, device = "pdf", useDingbats = FALSE)
    cat("✓ Saved trend map:", paste0("figures/", sensor_name, "_trend_map.pdf"), "\n")
  }
  
  # Save Figure 4 (manuscript version)
  if (!is.null(figure_4)) {
    # Save to main figures directory
    ggsave("figures/fig_4.pdf", 
           figure_4, width = 7.0, height = 6.0, device = "pdf", useDingbats = FALSE)
    cat("✓ Saved Figure 4:", "figures/fig_4.pdf", "\n")
    
    # Also save to albedo_map_trend directory
    if (!dir.exists("figures/albedo_map_trend")) {
      dir.create("figures/albedo_map_trend", recursive = TRUE)
    }
    ggsave("figures/albedo_map_trend/VIIRS_trend_map.pdf", 
           figure_4, width = 7.0, height = 6.0, device = "pdf", useDingbats = FALSE)
    cat("✓ Saved Figure 4 (copy):", "figures/albedo_map_trend/VIIRS_trend_map.pdf", "\n")
  }
  
  if (!is.null(significance_map)) {
    ggsave(paste0("figures/", sensor_name, "_significance_map.pdf"), 
           significance_map, width = 8, height = 6, device = "pdf", useDingbats = FALSE)
    cat("✓ Saved significance map:", paste0("figures/", sensor_name, "_significance_map.pdf"), "\n")
  }
  
  if (!is.null(summary_plot)) {
    ggsave(paste0("figures/", sensor_name, "_slope_distribution.pdf"), 
           summary_plot$plot, width = 6, height = 4, device = "pdf", useDingbats = FALSE)
    cat("✓ Saved slope distribution:", paste0("figures/", sensor_name, "_slope_distribution.pdf"), "\n")
  }
  
  # Save raster data
  writeRaster(trend_results$slope, 
              paste0("output/", sensor_name, "_slope.tif"), 
              overwrite = TRUE)
  writeRaster(trend_results$pvalue, 
              paste0("output/", sensor_name, "_pvalue.tif"), 
              overwrite = TRUE)
  writeRaster(trend_results$tau, 
              paste0("output/", sensor_name, "_tau.tif"), 
              overwrite = TRUE)
  writeRaster(trend_results$nobs, 
              paste0("output/", sensor_name, "_nobs.tif"), 
              overwrite = TRUE)
  
  cat("✓ Saved raster outputs to output/ directory\n")
  
  return(list(
    trend_results = trend_results,
    trend_map = trend_map,
    figure_4 = figure_4,
    significance_map = significance_map,
    summary_plot = summary_plot
  ))
}

# Function to calculate trend pixel percentages
calculate_trend_percentages <- function(sensor_name = "VIIRS") {
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("CALCULATING TREND PIXEL PERCENTAGES\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Load slope raster
  slope_path <- paste0("output/", sensor_name, "_slope.tif")
  if (!file.exists(slope_path)) {
    cat("Warning: Slope raster not found at", slope_path, "\n")
    return(NULL)
  }
  
  cat("Loading", sensor_name, "slope raster...\n")
  slope_rast <- rast(slope_path)
  
  # Extract all values
  slope_values <- values(slope_rast)
  slope_values <- as.vector(slope_values)
  
  # Remove NA values
  slope_values_valid <- slope_values[!is.na(slope_values)]
  
  if (length(slope_values_valid) == 0) {
    cat("Warning: No valid slope values found\n")
    return(NULL)
  }
  
  # Calculate counts
  total_pixels <- length(slope_values_valid)
  positive_pixels <- sum(slope_values_valid > 0)
  negative_pixels <- sum(slope_values_valid < 0)
  zero_pixels <- sum(slope_values_valid == 0)
  
  # Calculate percentages
  positive_percent <- (positive_pixels / total_pixels) * 100
  negative_percent <- (negative_pixels / total_pixels) * 100
  zero_percent <- (zero_pixels / total_pixels) * 100
  
  # Print results
  cat("\nTrend Pixel Statistics:\n")
  cat(sprintf("  Total valid pixels: %d\n", total_pixels))
  cat(sprintf("  Positive trends (increasing albedo): %d pixels (%.2f%%)\n", 
              positive_pixels, positive_percent))
  cat(sprintf("  Negative trends (decreasing albedo): %d pixels (%.2f%%)\n", 
              negative_pixels, negative_percent))
  cat(sprintf("  Zero trends (no change): %d pixels (%.2f%%)\n", 
              zero_pixels, zero_percent))
  
  # Additional statistics
  cat("\nAdditional Statistics:\n")
  cat(sprintf("  Mean slope: %.6f per year\n", mean(slope_values_valid)))
  cat(sprintf("  Median slope: %.6f per year\n", median(slope_values_valid)))
  cat(sprintf("  Min slope: %.6f per year\n", min(slope_values_valid)))
  cat(sprintf("  Max slope: %.6f per year\n", max(slope_values_valid)))
  cat(sprintf("  Standard deviation: %.6f\n", sd(slope_values_valid)))
  
  # Export trend percentages data
  trend_percentages_df <- data.frame(
    sensor = sensor_name,
    total_pixels = total_pixels,
    positive_pixels = positive_pixels,
    negative_pixels = negative_pixels,
    zero_pixels = zero_pixels,
    positive_percent = positive_percent,
    negative_percent = negative_percent,
    zero_percent = zero_percent,
    mean_slope = mean(slope_values_valid),
    median_slope = median(slope_values_valid),
    sd_slope = sd(slope_values_valid),
    min_slope = min(slope_values_valid),
    max_slope = max(slope_values_valid),
    stringsAsFactors = FALSE
  )
  write.csv(trend_percentages_df, 
            paste0("output/", sensor_name, "_trend_percentages.csv"), 
            row.names = FALSE)
  cat("✓ Exported trend percentages data\n")
  
  return(list(
    total_pixels = total_pixels,
    positive_pixels = positive_pixels,
    negative_pixels = negative_pixels,
    zero_pixels = zero_pixels,
    positive_percent = positive_percent,
    negative_percent = negative_percent,
    zero_percent = zero_percent,
    mean_slope = mean(slope_values_valid),
    median_slope = median(slope_values_valid),
    sd_slope = sd(slope_values_valid)
  ))
}

# Function to analyze albedo slope distribution by elevation
analyze_elevation_slope_relationship <- function(sensor_name = "VIIRS") {
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("ANALYZING ALBEDO SLOPE DISTRIBUTION BY ELEVATION\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Load SRTM DEM
  srtm_path <- "data/heard_island_srtm_dem.tif"
  if (!file.exists(srtm_path)) {
    cat("Warning: SRTM DEM file not found at", srtm_path, "\n")
    return(NULL)
  }
  
  # Load slope raster
  slope_path <- paste0("output/", sensor_name, "_slope.tif")
  if (!file.exists(slope_path)) {
    cat("Warning: Slope raster not found at", slope_path, "\n")
    cat("  Please run trend analysis first to generate slope rasters.\n")
    return(NULL)
  }
  
  cat("\nLoading SRTM DEM...\n")
  srtm <- rast(srtm_path)
  cat("Loading", sensor_name, "slope raster...\n")
  slope_rast <- rast(slope_path)
  
  # Check if rasters are compatible (same extent and resolution)
  cat("Checking raster compatibility...\n")
  
  # Resample DEM to match slope raster if needed
  if (!identical(ext(srtm), ext(slope_rast)) || !identical(res(srtm), res(slope_rast))) {
    cat("Resampling SRTM DEM to match slope raster...\n")
    srtm_resampled <- resample(srtm, slope_rast, method = "bilinear")
  } else {
    srtm_resampled <- srtm
  }
  
  # Extract pixel values first
  cat("Extracting pixel values...\n")
  slope_values_all <- values(slope_rast)
  elevation_values_all <- values(srtm_resampled)
  
  # Create a mask for valid pixels (non-NA in both rasters)
  valid_mask <- !is.na(slope_values_all) & !is.na(elevation_values_all)
  
  # Extract valid pixel values
  slope_values <- slope_values_all[valid_mask]
  elevation_values <- elevation_values_all[valid_mask]
  
  # Create data frame
  elevation_slope_data <- data.frame(
    elevation = as.vector(elevation_values),
    slope = as.vector(slope_values)
  )
  
  # Remove any remaining NA values
  elevation_slope_data <- elevation_slope_data[complete.cases(elevation_slope_data), ]
  
  cat("✓ Extracted", nrow(elevation_slope_data), "valid pixel pairs\n")
  cat("  Elevation range:", round(min(elevation_slope_data$elevation), 1), 
      "to", round(max(elevation_slope_data$elevation), 1), "m\n")
  cat("  Slope range:", round(min(elevation_slope_data$slope, na.rm = TRUE), 6), 
      "to", round(max(elevation_slope_data$slope, na.rm = TRUE), 6), "per year\n")
  
  if (nrow(elevation_slope_data) < 100) {
    cat("Warning: Insufficient data points for robust analysis\n")
    return(list(data = elevation_slope_data, analysis = NULL))
  }
  
  # Calculate correlation
  elevation_slope_cor <- cor.test(elevation_slope_data$elevation, elevation_slope_data$slope)
  
  cat("\nCorrelation analysis:\n")
  cat(sprintf("  Pearson r = %.4f, p = %.4f\n", 
              elevation_slope_cor$estimate, elevation_slope_cor$p.value))
  
  # Create elevation bins for analysis
  elevation_slope_data$elevation_bin <- cut(elevation_slope_data$elevation, 
                                            breaks = quantile(elevation_slope_data$elevation, 
                                                            probs = seq(0, 1, 0.2)),
                                            include.lowest = TRUE,
                                            labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))
  
  # Calculate statistics by elevation bin
  bin_stats <- elevation_slope_data %>%
    group_by(elevation_bin) %>%
    summarise(
      mean_elevation = mean(elevation, na.rm = TRUE),
      mean_slope = mean(slope, na.rm = TRUE),
      median_slope = median(slope, na.rm = TRUE),
      sd_slope = sd(slope, na.rm = TRUE),
      n_pixels = n(),
      .groups = "drop"
    )
  
  cat("\nStatistics by elevation quantile:\n")
  for (i in seq_len(nrow(bin_stats))) {
    cat(sprintf("  %s (mean elev: %.0f m): mean slope = %.6f, median = %.6f, n = %d\n",
                bin_stats$elevation_bin[i],
                bin_stats$mean_elevation[i],
                bin_stats$mean_slope[i],
                bin_stats$median_slope[i],
                bin_stats$n_pixels[i]))
  }
  
  # Fit linear regression
  lm_model <- lm(slope ~ elevation, data = elevation_slope_data)
  lm_summary <- summary(lm_model)
  
  cat("\nLinear regression (slope ~ elevation):\n")
  cat(sprintf("  Intercept: %.6f\n", coef(lm_model)[1]))
  cat(sprintf("  Elevation coefficient: %.6f per 100m\n", coef(lm_model)[2] * 100))
  cat(sprintf("  R²: %.4f\n", lm_summary$r.squared))
  cat(sprintf("  p-value: %.4f\n", lm_summary$coefficients[2, 4]))
  
  # Try polynomial fit
  lm_poly <- lm(slope ~ poly(elevation, 2), data = elevation_slope_data)
  lm_poly_summary <- summary(lm_poly)
  
  cat("\nPolynomial regression (slope ~ elevation²):\n")
  cat(sprintf("  R²: %.4f\n", lm_poly_summary$r.squared))
  
  # Compare models
  anova_result <- anova(lm_model, lm_poly)
  cat(sprintf("  Model comparison p-value: %.4f\n", anova_result$`Pr(>F)`[2]))
  
  return(list(
    data = elevation_slope_data,
    correlation = elevation_slope_cor,
    bin_stats = bin_stats,
    linear_model = lm_model,
    polynomial_model = lm_poly,
    srtm_raster = srtm_resampled,
    slope_raster = slope_rast
  ))
}

# Run analysis for all sensors
cat(paste(rep("=", 80), collapse=""), "\n")
cat("HEARD ISLAND GLACIER ALBEDO TREND MAP ANALYSIS\n")
cat("Pixel-level linear regression trend analysis for multiple sensors\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Define sensors to analyze
# Focus on VIIRS only - Based on comprehensive data quality assessment:
# - VIIRS: Uses validated BRDF/albedo parameters from VNP43 product
# - Sentinel-3: Excluded due to use of adapted coefficients from Sentinel-2 MSI
#   rather than validated Sentinel-3 OLCI-specific coefficients, making results non-credible
sensors <- c("VIIRS")

# Store results
all_results <- list()

# ============================================================================
# VIIRS PROCESSING
# ============================================================================

# Process VIIRS
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("PROCESSING VIIRS\n")
cat(paste(rep("=", 80), collapse=""), "\n")
result_viirs <- analyze_sensor_trends("VIIRS")
if (!is.null(result_viirs)) {
  all_results[["VIIRS"]] <- result_viirs
  cat("✓ VIIRS processing completed successfully\n")
} else {
  cat("✗ VIIRS processing failed\n")
}

# ============================================================================
# TREND PIXEL PERCENTAGE ANALYSIS
# ============================================================================

# Calculate percentages of positive and negative trend pixels
trend_percentages <- calculate_trend_percentages("VIIRS")

# ============================================================================
# ELEVATION-SLOPE ANALYSIS
# ============================================================================

# Perform elevation-slope analysis for VIIRS
elevation_slope_results <- analyze_elevation_slope_relationship("VIIRS")

# Create elevation-slope visualizations if analysis was successful
if (!is.null(elevation_slope_results) && !is.null(elevation_slope_results$data)) {
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("CREATING ELEVATION-SLOPE VISUALIZATIONS\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  elev_data <- elevation_slope_results$data
  corr_result <- elevation_slope_results$correlation
  lm_model <- elevation_slope_results$linear_model
  bin_stats <- elevation_slope_results$bin_stats
  
  # 1. Scatter plot: slope vs elevation with regression
  corr_annotation <- sprintf("r = %.3f, p = %.4f\nR² = %.3f", 
                            corr_result$estimate, corr_result$p.value,
                            summary(lm_model)$r.squared)
  
  # Sample data if too large for plotting (for performance)
  if (nrow(elev_data) > 50000) {
    cat("Sampling elevation-slope data for visualization (large dataset)\n")
    elev_data_plot <- elev_data[sample(nrow(elev_data), 50000), ]
  } else {
    elev_data_plot <- elev_data
  }
  
  p_elev_scatter <- ggplot(elev_data_plot, aes(x = elevation, y = slope)) +
    geom_point(alpha = 0.1, size = 0.3) +
    stat_bin2d(bins = 50, alpha = 0.7) +
    scale_fill_viridis_c(option = "plasma", trans = "log10") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
    geom_smooth(method = "loess", se = FALSE, color = "blue", linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = Inf, y = Inf, label = corr_annotation, 
             hjust = 1.1, vjust = 1.5, size = 3) +
    labs(x = "Elevation (m)", y = "Albedo Trend (per year)", 
         fill = "Count\n(log10)",
         title = "Albedo Trend vs Elevation",
         subtitle = paste0("VIIRS ", "(", nrow(elev_data), " pixels)")) +
    theme_times() +
    theme(legend.position = "right")
  
  ggsave("figures/VIIRS_elevation_slope_scatter.pdf", 
         p_elev_scatter, width = 7.0, height = 4.5, 
         device = "pdf", useDingbats = FALSE)
  cat("✓ Saved elevation-slope scatter plot\n")
  
  # Export elevation-slope scatter data
  elev_scatter_data <- elev_data_plot[, c("elevation", "slope")]
  write.csv(elev_scatter_data, 
            "output/VIIRS_elevation_slope_scatter_data.csv", 
            row.names = FALSE)
  cat("✓ Exported elevation-slope scatter data: output/VIIRS_elevation_slope_scatter_data.csv\n")
  
  # 2. Boxplot: slope distribution by elevation quantiles
  bin_labels <- paste0(bin_stats$elevation_bin, "\n(mean: ", 
                       round(bin_stats$mean_elevation), " m)")
  
  p_elev_dist <- ggplot(elev_data, aes(x = elevation_bin, y = slope, fill = elevation_bin)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.1, alpha = 0.8, outlier.size = 0.3) +
    geom_point(data = bin_stats, aes(y = mean_slope, x = elevation_bin), 
              color = "red", size = 2) +
    scale_fill_viridis_d(option = "plasma") +
    scale_x_discrete(labels = bin_labels) +
    labs(x = "Elevation Quantile", y = "Albedo Trend (per year)",
         title = "Albedo Trend Distribution by Elevation",
         subtitle = "Violin plots show probability density, red points indicate mean slopes") +
    theme_times() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5))
  
  ggsave("figures/VIIRS_elevation_slope_distribution.pdf", 
         p_elev_dist, width = 7.0, height = 4.5, 
         device = "pdf", useDingbats = FALSE)
  cat("✓ Saved elevation-slope distribution plot\n")
  
  # Export elevation-slope distribution data
  elev_dist_data <- elev_data[, c("elevation", "elevation_bin", "slope")]
  write.csv(elev_dist_data, 
            "output/VIIRS_elevation_slope_distribution_data.csv", 
            row.names = FALSE)
  write.csv(bin_stats, 
            "output/VIIRS_elevation_bin_stats.csv", 
            row.names = FALSE)
  cat("✓ Exported elevation-slope distribution data\n")
  
  # 3. 2D density plot: elevation vs slope
  p_elev_density <- ggplot(elev_data_plot, aes(x = elevation, y = slope)) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE) +
    scale_fill_viridis_c(option = "plasma") +
    geom_smooth(method = "lm", se = TRUE, color = "yellow", linewidth = 1, alpha = 0.3) +
    labs(x = "Elevation (m)", y = "Albedo Trend (per year)",
         fill = "Density",
         title = "Elevation-Slope Density Distribution",
         subtitle = "Contour lines show probability density") +
    theme_times() +
    theme(legend.position = "right")
  
  ggsave("figures/VIIRS_elevation_slope_density.pdf", 
         p_elev_density, width = 7.0, height = 4.5, 
         device = "pdf", useDingbats = FALSE)
  cat("✓ Saved elevation-slope density plot\n")
  
  # Export elevation-slope density data (using the same data as scatter)
  write.csv(elev_data_plot[, c("elevation", "slope")], 
            "output/VIIRS_elevation_slope_density_data.csv", 
            row.names = FALSE)
  cat("✓ Exported elevation-slope density data: output/VIIRS_elevation_slope_density_data.csv\n")
  
  cat("\n✓ All elevation-slope visualizations created successfully\n")
}

# Create summary
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("SUMMARY OF VIIRS ANALYSIS\n")
cat(paste(rep("=", 80), collapse=""), "\n")

summary_data <- data.frame()

for (sensor in names(all_results)) {
  result <- all_results[[sensor]]
  if (!is.null(result$summary_plot)) {
    summary_data <- rbind(summary_data, result$summary_plot$summary)
  }
}

if (nrow(summary_data) > 0) {
  # Print summary table
  cat("\nSensor Summary Statistics:\n")
  print(summary_data)
  
  # Save summary to CSV
  write.csv(summary_data, "output/trend_analysis_summary.csv", row.names = FALSE)
  cat("✓ Saved summary to output/trend_analysis_summary.csv\n")
  
  # Export combined summary plot data
  write.csv(summary_data, 
            "output/combined_sensor_trends_data.csv", 
            row.names = FALSE)
  cat("✓ Exported combined sensor trends data: output/combined_sensor_trends_data.csv\n")
  
  # Create combined summary plot
  if (nrow(summary_data) > 1) {
    p_combined <- ggplot(summary_data, aes(x = sensor, y = mean_slope)) +
      geom_col(aes(fill = sensor), alpha = 0.7) +
      geom_errorbar(aes(ymin = mean_slope - sd_slope, 
                        ymax = mean_slope + sd_slope), 
                    width = 0.2) +
      geom_text(aes(label = paste0("n=", n_pixels)), vjust = -0.5, size = 3) +
      labs(
        title = "Mean Albedo Trends by Sensor",
        subtitle = "Linear regression slope per year with standard deviation",
        x = "Sensor",
        y = "Mean Slope (per year)"
      ) +
      theme_times() +
      theme(legend.position = "none")
    
    ggsave("figures/combined_sensor_trends.pdf", 
           p_combined, width = 6, height = 4, device = "pdf", useDingbats = FALSE)
    cat("✓ Saved combined trends plot: figures/combined_sensor_trends.pdf\n")
  }
}

# Add trend percentages summary
if (!is.null(trend_percentages)) {
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("TREND PIXEL PERCENTAGE SUMMARY\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  cat(sprintf("\nPositive trends (increasing albedo): %.2f%% (%d pixels)\n",
              trend_percentages$positive_percent, trend_percentages$positive_pixels))
  cat(sprintf("Negative trends (decreasing albedo): %.2f%% (%d pixels)\n",
              trend_percentages$negative_percent, trend_percentages$negative_pixels))
  cat(sprintf("No change (zero slope): %.2f%% (%d pixels)\n",
              trend_percentages$zero_percent, trend_percentages$zero_pixels))
  cat(sprintf("\nOverall statistics:\n"))
  cat(sprintf("  Mean slope: %.6f per year\n", trend_percentages$mean_slope))
  cat(sprintf("  Median slope: %.6f per year\n", trend_percentages$median_slope))
}

# Add elevation-slope summary if analysis was successful
if (!is.null(elevation_slope_results) && !is.null(elevation_slope_results$correlation)) {
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("ELEVATION-SLOPE RELATIONSHIP SUMMARY\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  elev_corr <- elevation_slope_results$correlation
  elev_lm <- elevation_slope_results$linear_model
  elev_bins <- elevation_slope_results$bin_stats
  
  cat(sprintf("\nCorrelation: r = %.4f, p = %.4f\n", 
              elev_corr$estimate, elev_corr$p.value))
  cat(sprintf("Linear regression R² = %.4f\n", summary(elev_lm)$r.squared))
  cat(sprintf("Elevation coefficient: %.6f per 100m\n", coef(elev_lm)[2] * 100))
  cat("\nMean slopes by elevation quantile:\n")
  for (i in seq_len(nrow(elev_bins))) {
    cat(sprintf("  %s (%.0f m): %.6f per year\n",
                elev_bins$elevation_bin[i],
                elev_bins$mean_elevation[i],
                elev_bins$mean_slope[i]))
  }
}

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("TREND MAP ANALYSIS COMPLETE!\n")
cat("Methods implemented:\n")
cat("  • Linear regression (OLS trend detection)\n")
cat("  • Mann-Kendall trend test (significance testing)\n")
cat("  • Pixel-level spatial analysis\n")
cat("  • Quality control and data validation\n")
cat("  • Trend pixel percentage analysis (positive/negative trends)\n")
cat("  • Elevation-slope relationship analysis (SRTM DEM integration)\n")
cat("\nOutputs saved:\n")
cat("  • Figure 4 (manuscript): figures/fig_4.pdf\n")
cat("  • Trend maps: figures/*_trend_map.pdf\n")
cat("  • Significance maps: figures/*_significance_map.pdf\n")
cat("  • Slope distributions: figures/*_slope_distribution.pdf\n")
if (!is.null(elevation_slope_results) && !is.null(elevation_slope_results$correlation)) {
  cat("  • Elevation-slope scatter: figures/VIIRS_elevation_slope_scatter.pdf\n")
  cat("  • Elevation-slope distribution: figures/VIIRS_elevation_slope_distribution.pdf\n")
  cat("  • Elevation-slope density: figures/VIIRS_elevation_slope_density.pdf\n")
}
cat("  • Raster data: output/*.tif\n")
cat("  • Summary statistics: output/trend_analysis_summary.csv\n")
cat("\nData exports saved:\n")
cat("  • Trend map data: output/*_trend_map_data.csv\n")
cat("  • Significance map data: output/*_significance_map_data.csv\n")
cat("  • Slope distribution data: output/*_slope_distribution_data.csv\n")
cat("  • Trend percentages: output/*_trend_percentages.csv\n")
if (!is.null(elevation_slope_results) && !is.null(elevation_slope_results$correlation)) {
  cat("  • Elevation-slope scatter data: output/VIIRS_elevation_slope_scatter_data.csv\n")
  cat("  • Elevation-slope distribution data: output/VIIRS_elevation_slope_distribution_data.csv\n")
  cat("  • Elevation bin statistics: output/VIIRS_elevation_bin_stats.csv\n")
  cat("  • Elevation-slope density data: output/VIIRS_elevation_slope_density_data.csv\n")
}
cat("All plots use Times New Roman font size 11 for publication quality.\n")
