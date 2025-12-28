# Script to Create Figure 4: Albedo Trend Map with SRTM DEM Contour Lines
# Heard Island Glacier Albedo Analysis
# This script generates Figure 4 showing pixel-level albedo trends with elevation contours
#
# Author: Automated script for HIG project
# Last updated: 2024

# =============================================================================
# CONFIGURATION
# =============================================================================

# Set working directory (adjust if needed)
setwd("")

# File paths
sensor_name <- "VIIRS"
slope_raster_path <- "output/VIIRS_slope.tif"
pvalue_raster_path <- "output/VIIRS_pvalue.tif"
srtm_dem_path <- "data/heard_island_srtm_dem.tif"
glacier_outline_path <- "shapefiles/hig/hig_dissolved_noholes.shp"
output_path <- "figures/fig_4.pdf"
output_dir <- "figures/fig_4"

# Glacier outline settings
glacier_outline_color <- "blue"  # Changed to blue as requested
glacier_outline_width <- 0.5
glacier_outline_alpha <- 0.8
show_glacier_outline <- TRUE  # Set to FALSE to hide outline

# Zero-trend line settings (separates positive from negative trends)
show_zero_trend_line <- TRUE  # Set to FALSE to hide zero-trend line
zero_trend_line_color <- "black"
zero_trend_line_width <- 0.5  # Slightly thicker for visibility
zero_trend_line_alpha <- 0.9
zero_trend_line_linetype <- "solid"
zero_trend_smooth <- TRUE  # Smooth the boundary for rougher separation
zero_trend_tolerance <- 0.001  # Tolerance for smoothing (affects roughness)

# Plot settings
fig_width <- 7.0
fig_height <- 6.0

# Contour settings
contour_interval <- 200  # Elevation interval in meters (adjustable)
contour_color <- "gray30"
contour_alpha <- 0.4
contour_linewidth <- 0.3

# Label settings (if metR is available)
label_size <- 3.5  # Increased for better visibility
label_color <- "black"  # Changed to black for better contrast
label_alpha <- 0.9  # Text opacity (0-1)
label_family <- "Times New Roman"  # Font family for labels
label_stroke <- 0  # Stroke width (0 = transparent/no background)
label_stroke_color <- "transparent"  # Transparent background
label_format <- function(x) paste0(x, " m")  # Format: "200 m" instead of "200"

# Export summary statistics
export_summary <- TRUE
summary_path <- "output/fig_4_summary_statistics.csv"

# =============================================================================
# LOAD LIBRARIES
# =============================================================================

# Required libraries
required_packages <- c("raster", "terra", "ggplot2", "extrafont", "sf")

cat("Checking required packages...\n")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), 
       "\nInstall with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n")
}

library(raster)
library(terra)
library(ggplot2)
library(extrafont)
library(sf)
library(metR)

# Optional package for contour labels
cat("Checking optional packages...\n")
if (!requireNamespace("metR", quietly = TRUE)) {
  cat("  Note: metR package not found. Contour labels will not be added.\n")
  cat("  Install with: install.packages('metR')\n")
  use_contour_labels <- FALSE
} else {
  library(metR)
  use_contour_labels <- TRUE
  cat("  ✓ metR package available - contour labels will be added\n")
}

# Load Times New Roman font with error handling
tryCatch({
  loadfonts(device = "pdf", quiet = TRUE)
  cat("✓ Times New Roman font loaded successfully\n")
}, error = function(e) {
  cat("Warning: Could not load Times New Roman font\n")
  cat("Using system default font instead\n")
})

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Custom theme function
theme_times <- function() {
  # Check if Times New Roman is available
  if ("Times New Roman" %in% fonts()) {
    font_family <- "Times New Roman"
  } else {
    font_family <- "serif"  # Fallback to serif font
  }
  
  theme_minimal() +
    theme(
      text = element_text(family = font_family, size = 11),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

# Function to check file existence with informative error
check_file <- function(file_path, file_type = "file") {
  if (!file.exists(file_path)) {
    stop("\nERROR: ", file_type, " not found at:\n  ", file_path, 
         "\n\nPlease ensure the file exists or update the path in the configuration section.\n",
         call. = FALSE)
  }
  return(TRUE)
}

# =============================================================================
# MAIN SCRIPT
# =============================================================================

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("CREATING FIGURE 4: ALBEDO TREND MAP WITH SRTM DEM CONTOURS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n\n")
}

# -----------------------------------------------------------------------------
# Load and prepare data
# -----------------------------------------------------------------------------

cat("Step 1: Loading raster data...\n")
cat("  Loading slope raster from:", slope_raster_path, "\n")
check_file(slope_raster_path, "Slope raster")
slope_raster <- raster(slope_raster_path)

cat("  Loading p-value raster from:", pvalue_raster_path, "\n")
check_file(pvalue_raster_path, "P-value raster")
pvalue_raster <- raster(pvalue_raster_path)

cat("  Converting rasters to data frames...\n")
slope_df <- as.data.frame(slope_raster, xy = TRUE)
pvalue_df <- as.data.frame(pvalue_raster, xy = TRUE)

# Merge data
plot_df <- merge(slope_df, pvalue_df, by = c("x", "y"))
names(plot_df) <- c("x", "y", "slope", "pvalue")

# Remove NA values
plot_df <- plot_df[!is.na(plot_df$slope) & !is.na(plot_df$pvalue), ]

if (nrow(plot_df) == 0) {
  stop("ERROR: No valid data found in rasters after removing NA values.\n",
       "Please check that the input rasters contain valid data.\n", call. = FALSE)
}

cat("  ✓ Loaded", format(nrow(plot_df), big.mark = ","), "valid pixels\n")
cat("  ✓ Slope range:", round(min(plot_df$slope), 6), "to", round(max(plot_df$slope), 6), "per year\n\n")

# -----------------------------------------------------------------------------
# Load and prepare SRTM DEM for contour lines
# -----------------------------------------------------------------------------

dem_df <- NULL
has_dem <- FALSE

if (file.exists(srtm_dem_path)) {
  cat("Step 2: Loading SRTM DEM for contour lines...\n")
  cat("  Loading DEM from:", srtm_dem_path, "\n")
  srtm <- rast(srtm_dem_path)
  
  # Convert slope raster to terra format for resampling
  slope_rast <- rast(slope_raster)
  
  # Check if resampling is needed
  ext_match <- identical(ext(srtm), ext(slope_rast))
  res_match <- identical(res(srtm), res(slope_rast))
  
  if (!ext_match || !res_match) {
    cat("  Resampling DEM to match slope raster resolution...\n")
    srtm_resampled <- resample(srtm, slope_rast, method = "bilinear")
    cat("  ✓ DEM resampled successfully\n")
  } else {
    srtm_resampled <- srtm
    cat("  ✓ DEM already matches slope raster resolution\n")
  }
  
  # Convert DEM to data frame for contours
  dem_df <- as.data.frame(srtm_resampled, xy = TRUE)
  names(dem_df) <- c("x", "y", "elevation")
  dem_df <- dem_df[!is.na(dem_df$elevation), ]
  
  if (nrow(dem_df) > 0) {
    has_dem <- TRUE
    elev_range <- range(dem_df$elevation, na.rm = TRUE)
    cat("  ✓ Loaded", format(nrow(dem_df), big.mark = ","), "DEM pixels\n")
    cat("  ✓ Elevation range:", round(elev_range[1]), "to", round(elev_range[2]), "m\n\n")
  } else {
    cat("  Warning: DEM data frame is empty after processing\n\n")
  }
} else {
  cat("Step 2: SRTM DEM not found\n")
  cat("  Warning: SRTM DEM file not found at:", srtm_dem_path, "\n")
  cat("  Figure will be created without contour lines.\n\n")
}

# -----------------------------------------------------------------------------
# Create the plot
# -----------------------------------------------------------------------------

cat("Step 3: Creating Figure 4...\n")

# Base plot with raster
p <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_raster(aes(fill = slope)) +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0,
    name = "Slope\n(per year)",
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Albedo Trend\n(per year)"
  ) +
  theme_times() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    legend.box = "vertical",  # Stack legends vertically
    legend.margin = margin(t = 10, r = 0, b = 0, l = 0)  # Add spacing between legends
  ) +
  coord_fixed()

# Initialize variables to track which elements are added for legend
has_zero_trend <- FALSE
has_elevation_contours <- FALSE
has_glacier_outline <- FALSE
contour_breaks <- NULL

# Add zero-trend contour line (separates positive from negative trends)
if (show_zero_trend_line) {
  cat("  Adding zero-trend line (separates positive/negative trends)...\n")
  
  # Create a simplified boundary by smoothing the zero contour
  # Use bin2D to create a coarser grid for smoother contours
  if (zero_trend_smooth && nrow(plot_df) > 1000) {
    # Create a coarser grid for smoother contours
    x_range <- range(plot_df$x, na.rm = TRUE)
    y_range <- range(plot_df$y, na.rm = TRUE)
    
    # Reduce resolution for smoother/rougher line
    nx <- 30  # Further reduced for rougher separation
    ny <- 30
    
    x_seq <- seq(x_range[1], x_range[2], length.out = nx)
    y_seq <- seq(y_range[1], y_range[2], length.out = ny)
    
    # Create grid and interpolate using spatial averaging for smoother boundaries
    grid_df <- expand.grid(x = x_seq, y = y_seq)
    
    # Use spatial averaging - take average slope within each grid cell
    grid_df$slope <- NA
    cell_size_x <- abs(diff(x_range)) / nx
    cell_size_y <- abs(diff(y_range)) / ny
    
    for (i in seq_len(nrow(grid_df))) {
      # Find points within this grid cell
      x_center <- grid_df$x[i]
      y_center <- grid_df$y[i]
      
      within_cell <- abs(plot_df$x - x_center) < cell_size_x * 0.6 &
                     abs(plot_df$y - y_center) < cell_size_y * 0.6
      
      if (sum(within_cell) > 0) {
        # Average slope within the cell for smoother boundary
        grid_df$slope[i] <- mean(plot_df$slope[within_cell], na.rm = TRUE)
      }
    }
    
    # Remove NA values
    grid_df <- grid_df[!is.na(grid_df$slope), ]
    
    if (nrow(grid_df) > 100) {
      # Use smoothed grid for contour
      p <- p + 
        geom_contour(
          data = grid_df,
          aes(x = x, y = y, z = slope, color = "Zero-trend line"),
          breaks = 0,  # Only draw contour at slope = 0
          linewidth = zero_trend_line_width,
          alpha = zero_trend_line_alpha,
          linetype = zero_trend_line_linetype,
          bins = 1,  # Single contour
          binwidth = NULL,
          inherit.aes = FALSE
        )
      cat("  ✓ Zero-trend line added (smoothed)\n\n")
    } else {
      # Fallback to regular contour
      p <- p + 
        geom_contour(
          data = plot_df,
          aes(x = x, y = y, z = slope, color = "Zero-trend line"),
          breaks = 0,
          linewidth = zero_trend_line_width,
          alpha = zero_trend_line_alpha,
          linetype = zero_trend_line_linetype,
          inherit.aes = FALSE
        )
      cat("  ✓ Zero-trend line added\n\n")
    }
  } else {
    # Regular contour without smoothing
    p <- p + 
      geom_contour(
        data = plot_df,
        aes(x = x, y = y, z = slope, color = "Zero-trend line"),
        breaks = 0,  # Only draw contour at slope = 0
        linewidth = zero_trend_line_width,
        alpha = zero_trend_line_alpha,
        linetype = zero_trend_line_linetype,
        inherit.aes = FALSE
      )
    
    cat("  ✓ Zero-trend line added\n\n")
  }
  
  has_zero_trend <- TRUE
}

# Add contour lines if DEM data is available
if (has_dem && !is.null(dem_df) && nrow(dem_df) > 0) {
  # Determine appropriate contour break intervals
  elev_range <- range(dem_df$elevation, na.rm = TRUE)
  elev_min <- floor(elev_range[1] / 100) * 100
  elev_max <- ceiling(elev_range[2] / 100) * 100
  
  # Create contour breaks
  contour_breaks <- seq(elev_min, elev_max, by = contour_interval)
  # Filter to only include breaks within actual elevation range
  contour_breaks <- contour_breaks[contour_breaks >= elev_range[1] & contour_breaks <= elev_range[2]]
  
  if (length(contour_breaks) > 0) {
    cat("  Adding contour lines at", contour_interval, "m intervals...\n")
    has_elevation_contours <- TRUE
    
    p <- p + 
      geom_contour(
        data = dem_df,
        aes(x = x, y = y, z = elevation, color = "Elevation contours"),
        breaks = contour_breaks,
        alpha = contour_alpha,
        linewidth = contour_linewidth,
        inherit.aes = FALSE
      )
    
    # Add elevation labels to contour lines if metR is available
    if (use_contour_labels) {
      cat("  Adding elevation labels to contour lines...\n")
      
      # Use geom_text_contour with improved spacing and placement
      # Try using label_placer_flattest for better placement, with fallback
      tryCatch({
        p <- p + 
          geom_text_contour(
            data = dem_df,
            aes(x = x, y = y, z = elevation),
            breaks = contour_breaks,
            size = label_size,
            color = label_color,
            alpha = label_alpha,
            family = label_family,  # Set font family
            stroke = label_stroke,  # 0 = transparent background
            stroke.color = label_stroke_color,  # Transparent background
            rotate = TRUE,
            skip = 2,  # Skip labels to reduce crowding (label every 3rd position)
            check_overlap = TRUE,  # Prevent overlapping labels
            min.size = 30,  # Minimum number of points required for labeling
            label.placer = label_placer_flattest(),  # Place labels on flattest parts
            inherit.aes = FALSE
          )
      }, error = function(e) {
        # Fallback if label_placer_flattest is not available
        cat("  Note: Using basic label placement\n")
        p <<- p + 
          geom_text_contour(
            data = dem_df,
            aes(x = x, y = y, z = elevation),
            breaks = contour_breaks,
            size = label_size,
            color = label_color,
            alpha = label_alpha,
            family = label_family,  # Set font family
            stroke = label_stroke,  # 0 = transparent background
            stroke.color = label_stroke_color,  # Transparent background
            rotate = TRUE,
            skip = 2,  # Skip labels to reduce crowding
            check_overlap = TRUE,
            min.size = 30,
            inherit.aes = FALSE
          )
      })
      
      cat("  ✓ Contour labels added (size:", label_size, ", color:", label_color, ")\n")
    } else {
      cat("  Note: Contour labels not added (metR package not available)\n")
      cat("  Install metR with: install.packages('metR')\n")
    }
    
    cat("  ✓ Added", length(contour_breaks), "contour lines (", 
        min(contour_breaks), "-", max(contour_breaks), "m)\n\n")
  }
} else {
  cat("  Note: Contour lines not added (DEM data not available)\n\n")
}

# -----------------------------------------------------------------------------
# Load and add glacier outline
# -----------------------------------------------------------------------------

glacier_outline <- NULL

if (show_glacier_outline && file.exists(glacier_outline_path)) {
  cat("Step 3b: Loading glacier outline...\n")
  cat("  Loading from:", glacier_outline_path, "\n")
  
  tryCatch({
    glacier_outline <- st_read(glacier_outline_path, quiet = TRUE)
    
    # Get CRS from slope raster and reproject if needed
    slope_crs <- crs(slope_raster)
    
    if (!is.na(slope_crs)) {
      # Convert raster CRS to sf CRS format
      slope_crs_string <- as.character(slope_crs@projargs)
      if (slope_crs_string != "" && !is.na(slope_crs_string)) {
        slope_crs_sf <- st_crs(slope_crs_string)
        glacier_crs <- st_crs(glacier_outline)
        
        # Reproject if CRS doesn't match
        if (!identical(glacier_crs, slope_crs_sf)) {
          cat("  Reprojecting glacier outline to match plot CRS...\n")
          glacier_outline <- st_transform(glacier_outline, crs = slope_crs_sf)
        }
      }
    }
    
    # Add glacier outline to plot
    has_glacier_outline <- TRUE
    p <- p + 
      geom_sf(
        data = glacier_outline,
        aes(color = "Glacier outline"),
        linewidth = glacier_outline_width,
        alpha = glacier_outline_alpha,
        fill = NA,  # No fill, just outline
        inherit.aes = FALSE
      )
    
    cat("  ✓ Glacier outline added to plot\n\n")
  }, error = function(e) {
    cat("  Warning: Could not load glacier outline:", e$message, "\n")
    cat("  Continuing without glacier outline...\n\n")
  })
} else if (show_glacier_outline) {
  cat("Step 3b: Glacier outline not found\n")
  cat("  Warning: Glacier outline file not found at:", glacier_outline_path, "\n")
  cat("  Continuing without glacier outline...\n\n")
}

# -----------------------------------------------------------------------------
# Create combined legend for all line elements
# -----------------------------------------------------------------------------

# Build combined legend values and labels
legend_values <- c()
legend_labels <- c()
legend_alphas <- c()
legend_widths <- c()

if (has_zero_trend && show_zero_trend_line) {
  legend_values <- c(legend_values, zero_trend_line_color)
  legend_labels <- c(legend_labels, "Zero-trend line")
  legend_alphas <- c(legend_alphas, zero_trend_line_alpha)
  legend_widths <- c(legend_widths, zero_trend_line_width)
}

if (has_elevation_contours && length(contour_breaks) > 0) {
  legend_values <- c(legend_values, contour_color)
  legend_labels <- c(legend_labels, "Elevation contours")
  legend_alphas <- c(legend_alphas, contour_alpha)
  legend_widths <- c(legend_widths, contour_linewidth)
}

if (has_glacier_outline && show_glacier_outline) {
  legend_values <- c(legend_values, glacier_outline_color)
  legend_labels <- c(legend_labels, "Glacier outline")
  legend_alphas <- c(legend_alphas, glacier_outline_alpha)
  legend_widths <- c(legend_widths, glacier_outline_width)
}

# Add combined legend if we have any line elements
if (length(legend_labels) > 0) {
  names(legend_values) <- legend_labels
  
  # Update the plot with combined legend
  # Use the last scale_color_manual, but combine all values
  p <- p + 
    scale_color_manual(
      name = "",
      values = legend_values,
      guide = guide_legend(
        override.aes = list(
          alpha = legend_alphas,
          linewidth = legend_widths
        ),
        order = 1
      )
    )
  
  cat("✓ Combined legend created for", length(legend_labels), "elements\n\n")
}

# -----------------------------------------------------------------------------
# Display and save the figure
# -----------------------------------------------------------------------------

cat("Step 4: Displaying and saving Figure 4...\n")
cat("  Output path:", output_path, "\n")
cat("  Dimensions:", fig_width, "x", fig_height, "inches\n")

# Display the plot in RStudio's Plots pane
cat("  Displaying plot in RStudio Plots pane...\n")
print(p)

cat("  ✓ Plot displayed in RStudio\n")

# Save the figure
cat("  Saving plot to file...\n")
tryCatch({
  ggsave(
    output_path,
    p,
    width = fig_width,
    height = fig_height,
    device = "pdf",
    useDingbats = FALSE
  )
  cat("  ✓ Figure 4 saved successfully!\n\n")
}, error = function(e) {
  stop("ERROR saving figure:\n  ", e$message, "\n\n", call. = FALSE)
})

# -----------------------------------------------------------------------------
# Calculate and print summary statistics
# -----------------------------------------------------------------------------

cat("Step 5: Summary statistics...\n")

# Calculate statistics
summary_stats <- data.frame(
  metric = c(
    "Total pixels",
    "Mean slope (per year)",
    "Median slope (per year)",
    "SD slope (per year)",
    "Min slope (per year)",
    "Max slope (per year)",
    "Pixels with negative trends",
    "Pixels with positive trends",
    "Pixels with zero trend",
    "Statistically significant (p < 0.05)",
    "Statistically significant negative (p < 0.05)",
    "Statistically significant positive (p < 0.05)"
  ),
  value = c(
    nrow(plot_df),
    round(mean(plot_df$slope), 6),
    round(median(plot_df$slope), 6),
    round(sd(plot_df$slope), 6),
    round(min(plot_df$slope), 6),
    round(max(plot_df$slope), 6),
    sum(plot_df$slope < 0),
    sum(plot_df$slope > 0),
    sum(plot_df$slope == 0),
    sum(plot_df$pvalue < 0.05),
    sum(plot_df$slope < 0 & plot_df$pvalue < 0.05),
    sum(plot_df$slope > 0 & plot_df$pvalue < 0.05)
  ),
  percentage = c(
    100,
    NA,
    NA,
    NA,
    NA,
    NA,
    round(sum(plot_df$slope < 0) / nrow(plot_df) * 100, 1),
    round(sum(plot_df$slope > 0) / nrow(plot_df) * 100, 1),
    round(sum(plot_df$slope == 0) / nrow(plot_df) * 100, 1),
    round(sum(plot_df$pvalue < 0.05) / nrow(plot_df) * 100, 1),
    round(sum(plot_df$slope < 0 & plot_df$pvalue < 0.05) / nrow(plot_df) * 100, 1),
    round(sum(plot_df$slope > 0 & plot_df$pvalue < 0.05) / nrow(plot_df) * 100, 1)
  )
)

# Print summary
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("FIGURE 4 SUMMARY STATISTICS\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat(sprintf("%-45s %12s %10s\n", "Metric", "Value", "Percent"))
cat(paste(rep("-", 70), collapse=""), "\n")

for (i in seq_len(nrow(summary_stats))) {
  if (is.na(summary_stats$percentage[i])) {
    cat(sprintf("%-45s %12s\n", 
                summary_stats$metric[i], 
                format(summary_stats$value[i], big.mark = ",")))
  } else {
    cat(sprintf("%-45s %12s %9.1f%%\n", 
                summary_stats$metric[i], 
                format(summary_stats$value[i], big.mark = ","),
                summary_stats$percentage[i]))
  }
}

cat(paste(rep("=", 70), collapse=""), "\n\n")

# Export summary if requested
if (export_summary) {
  # Create output directory if needed
  summary_dir <- dirname(summary_path)
  if (!dir.exists(summary_dir)) {
    dir.create(summary_dir, recursive = TRUE)
  }
  
  write.csv(summary_stats, summary_path, row.names = FALSE)
  cat("✓ Summary statistics exported to:", summary_path, "\n\n")
}

# -----------------------------------------------------------------------------
# Completion message
# -----------------------------------------------------------------------------

cat(paste(rep("=", 70), collapse=""), "\n")
cat("SCRIPT COMPLETED SUCCESSFULLY\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Figure 4 saved to:", output_path, "\n")
if (export_summary) {
  cat("Summary statistics saved to:", summary_path, "\n")
}
cat("\n")
