# Script to Create Figure 4: Albedo Trend Map (Multi-panel)
# Heard Island Glacier Albedo Analysis
# This script generates Figure 4 with three panels:
#   Panel a: Pixel-level albedo trend map (from figure_4_trend_map.R)
#   Panel b: Slope distribution histogram (from albedo_trend_map_linear.R)
#   Panel c: Elevation quartile analysis (from albedo_trend_map_linear.R)
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
elevation_bin_stats_path <- "output/VIIRS_elevation_bin_stats.csv"
elevation_slope_data_path <- "output/VIIRS_elevation_slope_distribution_data.csv"
output_path <- "figures/fig_4.pdf"
output_dir <- "figures/fig_4"

# Glacier outline settings
glacier_outline_color <- "blue"
glacier_outline_width <- 0.5
glacier_outline_alpha <- 0.8
show_glacier_outline <- TRUE

# Zero-trend line settings
show_zero_trend_line <- TRUE
zero_trend_line_color <- "black"
zero_trend_line_width <- 0.5
zero_trend_line_alpha <- 0.9
zero_trend_line_linetype <- "solid"
zero_trend_smooth <- TRUE

# Contour settings
contour_interval <- 200  # Elevation interval in meters
contour_color <- "gray30"
contour_alpha <- 0.4
contour_linewidth <- 0.3

# Label settings (if metR is available)
label_size <- 3.5
label_color <- "black"
label_alpha <- 0.9
label_family <- "Times New Roman"
label_stroke <- 0
label_stroke_color <- "transparent"

# Panel settings
fig_width <- 14.0  # Width for three panels side by side
fig_height <- 5.0  # Height for single row

# Also save as PNG for verification (optional)
save_png_version <- TRUE
output_path_png <- "figures/fig_4.png"

# Export summary statistics
export_summary <- TRUE
summary_path <- "output/fig_4_summary_statistics.csv"

# =============================================================================
# LOAD LIBRARIES
# =============================================================================

required_packages <- c("raster", "terra", "ggplot2", "extrafont", "sf", "patchwork", "dplyr", "viridis", "gridExtra")

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
library(patchwork)
library(dplyr)
library(viridis)
library(gridExtra)
library(grid)

# Optional package for contour labels
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
  if ("Times New Roman" %in% fonts()) {
    font_family <- "Times New Roman"
  } else {
    font_family <- "serif"
  }
  
  theme_minimal() +
    theme(
      text = element_text(family = font_family, size = 11),
      plot.margin = margin(5, 5, 5, 5, "mm"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid = element_blank()
    )
}

# Function to check file existence
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
cat("CREATING FIGURE 4: ALBEDO TREND MAP (MULTI-PANEL)\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n\n")
}

# -----------------------------------------------------------------------------
# Load and prepare data for Panel A and B
# -----------------------------------------------------------------------------

cat("Step 1: Loading raster data...\n")
check_file(slope_raster_path, "Slope raster")
check_file(pvalue_raster_path, "P-value raster")
slope_raster <- raster(slope_raster_path)
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
# Load SRTM DEM for contour lines
# -----------------------------------------------------------------------------

dem_df <- NULL
has_dem <- FALSE

if (file.exists(srtm_dem_path)) {
  cat("Step 2: Loading SRTM DEM for contour lines...\n")
  srtm <- rast(srtm_dem_path)
  slope_rast <- rast(slope_raster)
  
  if (!identical(ext(srtm), ext(slope_rast)) || !identical(res(srtm), res(slope_rast))) {
    cat("  Resampling DEM to match slope raster resolution...\n")
    srtm_resampled <- resample(srtm, slope_rast, method = "bilinear")
  } else {
    srtm_resampled <- srtm
  }
  
  dem_df <- as.data.frame(srtm_resampled, xy = TRUE)
  names(dem_df) <- c("x", "y", "elevation")
  dem_df <- dem_df[!is.na(dem_df$elevation), ]
  
  if (nrow(dem_df) > 0) {
    has_dem <- TRUE
    elev_range <- range(dem_df$elevation, na.rm = TRUE)
    cat("  ✓ Loaded", format(nrow(dem_df), big.mark = ","), "DEM pixels\n")
    cat("  ✓ Elevation range:", round(elev_range[1]), "to", round(elev_range[2]), "m\n\n")
  }
} else {
  cat("Step 2: SRTM DEM not found - contours will be omitted\n\n")
}

# -----------------------------------------------------------------------------
# Load glacier outline
# -----------------------------------------------------------------------------

glacier_outline <- NULL
has_glacier_outline <- FALSE

if (show_glacier_outline && file.exists(glacier_outline_path)) {
  cat("Step 3: Loading glacier outline...\n")
  tryCatch({
    glacier_outline <- st_read(glacier_outline_path, quiet = TRUE)
    slope_crs <- crs(slope_raster)
    
    if (!is.na(slope_crs)) {
      slope_crs_string <- as.character(slope_crs@projargs)
      if (slope_crs_string != "" && !is.na(slope_crs_string)) {
        slope_crs_sf <- st_crs(slope_crs_string)
        glacier_crs <- st_crs(glacier_outline)
        
        if (!identical(glacier_crs, slope_crs_sf)) {
          cat("  Reprojecting glacier outline to match plot CRS...\n")
          glacier_outline <- st_transform(glacier_outline, crs = slope_crs_sf)
        }
      }
    }
    
    has_glacier_outline <- TRUE
    cat("  ✓ Glacier outline loaded\n\n")
  }, error = function(e) {
    cat("  Warning: Could not load glacier outline:", e$message, "\n\n")
  })
}

# -----------------------------------------------------------------------------
# PANEL A: Trend Map (from figure_4_trend_map.R)
# -----------------------------------------------------------------------------

cat("Step 4: Creating Panel A (Trend Map)...\n")

# Base plot with raster
p_a <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_raster(aes(fill = slope)) +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0,
    name = "Slope\n(per year)",
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(5, "cm"),
      barheight = unit(0.4, "cm")
    )
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "(a) Albedo Trend Map"
  ) +
  theme_times() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.0, "cm"),
    legend.margin = margin(t = 5, r = 0, b = 0, l = 0)
  ) +
  coord_fixed()

# Add zero-trend contour line
if (show_zero_trend_line) {
  if (zero_trend_smooth && nrow(plot_df) > 1000) {
    x_range <- range(plot_df$x, na.rm = TRUE)
    y_range <- range(plot_df$y, na.rm = TRUE)
    nx <- 30
    ny <- 30
    
    x_seq <- seq(x_range[1], x_range[2], length.out = nx)
    y_seq <- seq(y_range[1], y_range[2], length.out = ny)
    grid_df <- expand.grid(x = x_seq, y = y_seq)
    
    grid_df$slope <- NA
    cell_size_x <- abs(diff(x_range)) / nx
    cell_size_y <- abs(diff(y_range)) / ny
    
    for (i in seq_len(nrow(grid_df))) {
      x_center <- grid_df$x[i]
      y_center <- grid_df$y[i]
      within_cell <- abs(plot_df$x - x_center) < cell_size_x * 0.6 &
                     abs(plot_df$y - y_center) < cell_size_y * 0.6
      
      if (sum(within_cell) > 0) {
        grid_df$slope[i] <- mean(plot_df$slope[within_cell], na.rm = TRUE)
      }
    }
    
    grid_df <- grid_df[!is.na(grid_df$slope), ]
    
    if (nrow(grid_df) > 100) {
      p_a <- p_a + 
        geom_contour(
          data = grid_df,
          aes(x = x, y = y, z = slope),
          breaks = 0,
          linewidth = zero_trend_line_width,
          alpha = zero_trend_line_alpha,
          linetype = zero_trend_line_linetype,
          color = zero_trend_line_color,
          inherit.aes = FALSE
        )
    }
  } else {
    p_a <- p_a + 
      geom_contour(
        data = plot_df,
        aes(x = x, y = y, z = slope),
        breaks = 0,
        linewidth = zero_trend_line_width,
        alpha = zero_trend_line_alpha,
        linetype = zero_trend_line_linetype,
        color = zero_trend_line_color,
        inherit.aes = FALSE
      )
  }
}

# Add elevation contours
if (has_dem && !is.null(dem_df) && nrow(dem_df) > 0) {
  elev_range <- range(dem_df$elevation, na.rm = TRUE)
  elev_min <- floor(elev_range[1] / 100) * 100
  elev_max <- ceiling(elev_range[2] / 100) * 100
  contour_breaks <- seq(elev_min, elev_max, by = contour_interval)
  contour_breaks <- contour_breaks[contour_breaks >= elev_range[1] & contour_breaks <= elev_range[2]]
  
  if (length(contour_breaks) > 0) {
    p_a <- p_a + 
      geom_contour(
        data = dem_df,
        aes(x = x, y = y, z = elevation),
        breaks = contour_breaks,
        alpha = contour_alpha,
        linewidth = contour_linewidth,
        color = contour_color,
        inherit.aes = FALSE
      )
    
    if (use_contour_labels) {
      tryCatch({
        p_a <- p_a + 
          geom_text_contour(
            data = dem_df,
            aes(x = x, y = y, z = elevation),
            breaks = contour_breaks,
            size = label_size,
            color = label_color,
            alpha = label_alpha,
            family = label_family,
            stroke = label_stroke,
            stroke.color = label_stroke_color,
            rotate = TRUE,
            skip = 2,
            check_overlap = TRUE,
            min.size = 30,
            inherit.aes = FALSE
          )
      }, error = function(e) {
        p_a <<- p_a + 
          geom_text_contour(
            data = dem_df,
            aes(x = x, y = y, z = elevation),
            breaks = contour_breaks,
            size = label_size,
            color = label_color,
            alpha = label_alpha,
            family = label_family,
            stroke = label_stroke,
            stroke.color = label_stroke_color,
            rotate = TRUE,
            skip = 2,
            check_overlap = TRUE,
            min.size = 30,
            inherit.aes = FALSE
          )
      })
    }
  }
}

# Add glacier outline
if (has_glacier_outline && show_glacier_outline) {
  p_a <- p_a + 
    geom_sf(
      data = glacier_outline,
      color = glacier_outline_color,
      linewidth = glacier_outline_width,
      alpha = glacier_outline_alpha,
      fill = NA,
      inherit.aes = FALSE
    )
}

cat("  ✓ Panel A created\n\n")

# -----------------------------------------------------------------------------
# PANEL B: Slope Distribution (from albedo_trend_map_linear.R)
# -----------------------------------------------------------------------------

cat("Step 5: Creating Panel B (Slope Distribution)...\n")

# Create plot data frame (same as in albedo_trend_map_linear.R)
slope_values <- plot_df$slope
pvalue_values <- plot_df$pvalue
plot_df_b <- data.frame(slope = slope_values, pvalue = pvalue_values)
plot_df_b$significant <- plot_df_b$pvalue < 0.05

# Calculate summary statistics
n_pixels <- length(slope_values)
significant_pixels <- sum(plot_df_b$significant)
significant_percent <- significant_pixels / n_pixels * 100

cat("  ✓ Total pixels:", n_pixels, "\n")
cat("  ✓ Significant (p < 0.05):", significant_pixels, "(", 
    round(significant_percent, 1), "%)\n\n")

# Create plot (matching albedo_trend_map_linear.R style)
p_b <- ggplot(plot_df_b, aes(x = slope)) +
  geom_histogram(aes(fill = significant), bins = 30, alpha = 0.7) +
  scale_fill_manual(
    values = c("FALSE" = "gray70", "TRUE" = "red"),
    name = "Significant\n(p < 0.05)",
    labels = c("FALSE" = "No", "TRUE" = "Yes")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "(b) Slope Distribution",
    subtitle = paste("n =", format(n_pixels, big.mark = ","), "pixels;", 
                     round(significant_percent, 1), "% significant"),
    x = "Linear Regression Slope (per year)",
    y = "Count"
  ) +
  theme_times() +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(size = 9)
  )

cat("  ✓ Panel B created\n\n")

# -----------------------------------------------------------------------------
# PANEL C: Elevation Quartile Analysis (from albedo_trend_map_linear.R)
# -----------------------------------------------------------------------------

cat("Step 6: Creating Panel C (Elevation Quartile Analysis)...\n")

# Load elevation bin statistics
check_file(elevation_bin_stats_path, "Elevation bin statistics")
bin_stats <- read.csv(elevation_bin_stats_path, stringsAsFactors = FALSE)

# Load elevation-slope distribution data
check_file(elevation_slope_data_path, "Elevation-slope distribution data")
elev_data <- read.csv(elevation_slope_data_path, stringsAsFactors = FALSE)

cat("  ✓ Loaded elevation statistics for", nrow(bin_stats), "quartiles\n")
cat("  ✓ Loaded elevation-slope data:", nrow(elev_data), "pixels\n\n")

# Create bin labels - shorter format to prevent overlap
bin_labels <- paste0(bin_stats$elevation_bin, "\n", 
                     round(bin_stats$mean_elevation), " m")

# Create plot (matching albedo_trend_map_linear.R style)
p_c <- ggplot(elev_data, aes(x = elevation_bin, y = slope, fill = elevation_bin)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.size = 0.3) +
  geom_point(data = bin_stats, aes(y = mean_slope, x = elevation_bin), 
            color = "red", size = 2) +
  scale_fill_viridis_d(option = "plasma") +
  scale_x_discrete(labels = bin_labels) +
  labs(
    x = "Elevation Quantile",
    y = "Albedo Trend (per year)",
    title = "(c) Elevation Quartile Analysis"
  ) +
  theme_times() +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 10, 5, 5, "mm"),
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5, vjust = 0.5)
  )

cat("  ✓ Panel C created\n\n")

# -----------------------------------------------------------------------------
# Combine panels
# -----------------------------------------------------------------------------

cat("Step 7: Combining panels...\n")

# Combine panels using patchwork - use wrap_plots for more control
# Give Panel C slightly more width to accommodate labels
combined_plot <- wrap_plots(p_a, p_b, p_c, ncol = 3, widths = c(1.1, 1, 1.1))

cat("  ✓ Panels combined\n\n")

# -----------------------------------------------------------------------------
# Display and save
# -----------------------------------------------------------------------------

cat("Step 8: Displaying and saving Figure 4...\n")
cat("  Output path:", output_path, "\n")
cat("  Dimensions:", fig_width, "x", fig_height, "inches\n")

# Display the plot
cat("  Displaying plot...\n")
print(combined_plot)

cat("  ✓ Plot displayed\n")

# Save the figure
# Use grid.arrange for PDF (more reliable than patchwork for PDF devices)
cat("  Saving plot to file...\n")

# Create arrangeGrob version for PDF saving (returns grob instead of printing)
cat("  Creating arrangeGrob layout for PDF...\n")
combined_plot_grid <- arrangeGrob(
  p_a, p_b, p_c,
  ncol = 3,
  widths = c(1.1, 1, 1.1)
)

# Method 1: Try cairo_pdf device with grid.arrange
if (capabilities("cairo")) {
  cat("  Attempting Method 1: cairo_pdf device with grid.arrange...\n")
  tryCatch({
    cairo_pdf(file = output_path, width = fig_width, height = fig_height, onefile = TRUE)
    grid.draw(combined_plot_grid)
    dev.off()
    cat("  ✓ Figure 4 saved successfully using cairo_pdf() with grid.arrange!\n\n")
  }, error = function(e) {
    cat("  Method 1 failed:", e$message, "\n")
    # Method 2: Try standard pdf device with grid.arrange
    cat("  Attempting Method 2: pdf() device with grid.arrange...\n")
    tryCatch({
      pdf(file = output_path, width = fig_width, height = fig_height, useDingbats = FALSE, onefile = TRUE)
      grid.draw(combined_plot_grid)
      dev.off()
      cat("  ✓ Figure 4 saved successfully using pdf() device with grid.arrange!\n\n")
    }, error = function(e2) {
      cat("  Method 2 failed:", e2$message, "\n")
      # Method 3: Try ggsave with patchwork (fallback)
      cat("  Attempting Method 3: ggsave() with patchwork...\n")
      tryCatch({
        ggsave(
          filename = output_path,
          plot = combined_plot,
          width = fig_width,
          height = fig_height,
          device = cairo_pdf,
          useDingbats = FALSE,
          limitsize = FALSE
        )
        cat("  ✓ Figure 4 saved successfully using ggsave() with cairo_pdf!\n\n")
      }, error = function(e3) {
        # Method 4: Final fallback - ggsave with standard pdf
        cat("  Attempting Method 4: ggsave() with standard pdf...\n")
        tryCatch({
          ggsave(
            filename = output_path,
            plot = combined_plot,
            width = fig_width,
            height = fig_height,
            device = "pdf",
            useDingbats = FALSE,
            limitsize = FALSE
          )
          cat("  ✓ Figure 4 saved successfully using ggsave() with pdf!\n\n")
        }, error = function(e4) {
          stop("ERROR: All save methods failed!\n",
               "  Method 1 (cairo_pdf grid.arrange): ", e$message, "\n",
               "  Method 2 (pdf grid.arrange): ", e2$message, "\n",
               "  Method 3 (ggsave cairo_pdf): ", e3$message, "\n",
               "  Method 4 (ggsave pdf): ", e4$message, "\n\n", call. = FALSE)
        })
      })
    })
  })
} else {
  # If cairo is not available, try standard methods with grid.arrange
  cat("  Cairo not available, using standard PDF methods with grid.arrange...\n")
  tryCatch({
    pdf(file = output_path, width = fig_width, height = fig_height, useDingbats = FALSE, onefile = TRUE)
    grid.draw(combined_plot_grid)
    dev.off()
    cat("  ✓ Figure 4 saved successfully using pdf() device with grid.arrange!\n\n")
  }, error = function(e) {
    cat("  pdf() device failed, trying ggsave()...\n")
    ggsave(
      filename = output_path,
      plot = combined_plot,
      width = fig_width,
      height = fig_height,
      device = "pdf",
      useDingbats = FALSE,
      limitsize = FALSE
    )
    cat("  ✓ Figure 4 saved successfully using ggsave()!\n\n")
  })
}

# Save PNG version for verification (helps diagnose PDF issues)
# Use patchwork for PNG (works fine) and grid.arrange for PDF
if (save_png_version) {
  cat("  Saving PNG version for verification...\n")
  tryCatch({
    png(file = output_path_png, width = fig_width * 300, height = fig_height * 300, res = 300)
    print(combined_plot)  # Use patchwork for PNG
    dev.off()
    cat("  ✓ PNG version saved to:", output_path_png, "\n")
    cat("  (PNG uses patchwork, PDF uses grid.arrange for compatibility)\n\n")
  }, error = function(e) {
    cat("  Warning: Could not save PNG version:", e$message, "\n\n")
  })
}

# -----------------------------------------------------------------------------
# Summary statistics
# -----------------------------------------------------------------------------

cat("Step 9: Summary statistics...\n")

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
cat("Figure 4 (multi-panel) saved to:", output_path, "\n")
if (export_summary) {
  cat("Summary statistics saved to:", summary_path, "\n")
}
cat("\n")

