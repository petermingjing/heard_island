# Create Graphic Abstract for Heard Island Glacier Albedo Study
# Journal of Advanced Research format: 1062 x 2656 pixels (height x width) - doubled size
# Created: 2024

# Load required libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)
library(png)  # For reading PNG images

# Set working directory
setwd("/Users/jingming/Desktop/HIG")

# Load Times New Roman font with error handling
tryCatch({
  loadfonts(device = "pdf", quiet = TRUE)
  cat("✓ Times New Roman font loaded successfully\n")
}, error = function(e) {
  cat("Warning: Could not load Times New Roman font:", e$message, "\n")
  cat("Using system default font instead\n")
})

# Define font family
get_font_family <- function() {
  if ("Times New Roman" %in% fonts()) {
    return("Times New Roman")
  } else {
    return("serif")
  }
}

font_family <- get_font_family()

# Color scheme
colors <- list(
  glacier = "#E8F4F8",
  declining = "#D73027",
  increasing = "#4575B4",
  radiation = "#FFA500",
  background = "#FFFFFF",
  text = "#2C3E50",
  accent = "#1A73E8",
  precipitation = "#4A90E2",
  temperature = "#E74C3C"
)

# ============================================================================
# LEFT PANEL: Study Area & Data Source
# ============================================================================

create_left_panel <- function() {
  # Load Heard Island image
  island_img_path <- "figures/fig_1/b.png"
  if (file.exists(island_img_path)) {
    island_img <- readPNG(island_img_path)
    use_image <- TRUE
  } else {
    use_image <- FALSE
    cat("Warning: Island image not found at", island_img_path, "\n")
  }
  
  p <- ggplot() +
    # Title
    annotate("text", x = 0.5, y = 0.98, label = "STUDY AREA & DATA",
             hjust = 0.5, vjust = 1, size = 3.2, 
             color = colors$text, fontface = "bold") +
    # Satellite representation (moved up)
    annotate("rect", xmin = 0.42, xmax = 0.58, ymin = 0.88, ymax = 0.92,
             fill = "#666666", color = "#333333", linewidth = 0.5) +
    annotate("rect", xmin = 0.40, xmax = 0.60, ymin = 0.86, ymax = 0.88,
             fill = colors$accent, color = "#0D47A1", linewidth = 0.5) +
    # Signal lines (shortened)
    annotate("segment", x = 0.47, xend = 0.47, y = 0.86, yend = 0.82,
             linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.5) +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.86, yend = 0.82,
             linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.5) +
    annotate("segment", x = 0.53, xend = 0.53, y = 0.86, yend = 0.82,
             linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.5) +
    annotate("text", x = 0.5, y = 0.80, label = "VIIRS Satellite",
             hjust = 0.5, vjust = 0, size = 2.0, color = colors$text, fontface = "italic")
  
  # Add Heard Island image or fallback polygon
  if(use_image) {
    p <- p + annotation_raster(island_img, xmin = 0.15, xmax = 0.85, ymin = 0.20, ymax = 0.65)
  } else {
    # Fallback to polygon if image not found
    island_x <- c(0.2, 0.5, 0.8, 0.7, 0.5, 0.3, 0.2)
    island_y <- c(0.25, 0.15, 0.30, 0.60, 0.70, 0.55, 0.25)
    island_df <- data.frame(x = island_x, y = island_y)
    p <- p + geom_polygon(data = island_df, aes(x = x, y = y), 
                          fill = colors$glacier, color = colors$text, linewidth = 0.5)
  }
  
  p <- p +
    # Location label with timeline in bracket (moved to bottom)
    annotate("text", x = 0.5, y = 0.12, 
             label = "Heard Island\n53°06'S, 73°31'E",
             hjust = 0.5, vjust = 0, size = 2.2, color = colors$text) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(0.05, 0.05, 0.05, 0.05, "inches"))
  
  return(p)
}

# ============================================================================
# CENTER PANEL: Methodology & Key Findings
# ============================================================================

create_center_panel <- function() {
  p <- ggplot() +
    # Title
    annotate("text", x = 0.5, y = 0.98, label = "KEY FINDINGS",
             hjust = 0.5, vjust = 1, size = 3.2, 
             color = colors$text, fontface = "bold") +
    # Input box (top, symmetric positioning)
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.80, ymax = 0.90,
             fill = "#F0F0F0", color = colors$text, linewidth = 0.5,
             linetype = "solid") +
    annotate("text", x = 0.5, y = 0.87, label = "VIIRS Albedo Observations",
             hjust = 0.5, vjust = 0.5, size = 2.6, 
             color = colors$text, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.83, 
             label = "Validated VNP43 (2012-2024)",
             hjust = 0.5, vjust = 0.5, size = 2.0, color = colors$text) +
    # Arrow 1 (equal spacing)
    annotate("segment", x = 0.5, xend = 0.5, y = 0.80, yend = 0.72,
             arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
             color = colors$text, linewidth = 0.5) +
    # Analysis box (center, symmetric size)
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.60, ymax = 0.70,
             fill = colors$glacier, color = colors$accent, linewidth = 0.5,
             linetype = "solid") +
    annotate("text", x = 0.5, y = 0.66, label = "Multi-Method Analysis",
             hjust = 0.5, vjust = 0.5, size = 2.6, 
             color = colors$text, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.62, 
             label = "Linear Regression + RF + GAM",
             hjust = 0.5, vjust = 0.5, size = 2.0, color = colors$text) +
    # Arrow 2 (equal spacing)
    annotate("segment", x = 0.5, xend = 0.5, y = 0.60, yend = 0.52,
             arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
             color = colors$text, linewidth = 0.5) +
    # Results box (bottom, symmetric positioning)
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.25, ymax = 0.50,
             fill = "#FFF4E6", color = colors$declining, linewidth = 0.5,
             linetype = "solid") +
    annotate("text", x = 0.5, y = 0.42, 
             label = "Albedo Decline: -0.015 per decade",
             hjust = 0.5, vjust = 0.5, size = 2.8, 
             color = colors$declining, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.36, 
             label = "77.3% of pixels show darkening",
             hjust = 0.5, vjust = 0.5, size = 2.4, color = colors$text) +
    annotate("text", x = 0.5, y = 0.30, 
             label = "Solar radiation: 53.5-60.3% variance",
             hjust = 0.5, vjust = 0.5, size = 2.0, color = colors$text) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(0.05, 0.05, 0.05, 0.05, "inches"))
  
  return(p)
}

# ============================================================================
# RIGHT PANEL: Climate Attribution
# ============================================================================

create_right_panel <- function() {
  # Variable importance data
  vars_df <- data.frame(
    variable = c("Solar\nRadiation", "Precipitation", "Surface\nTemperature"),
    importance = c(58, 25, 17),
    y_pos = c(0.70, 0.50, 0.30),
    color = c(colors$radiation, colors$precipitation, colors$temperature)
  )
  
  p <- ggplot(vars_df) +
    # Title
    annotate("text", x = 0.5, y = 0.98, label = "CLIMATE DRIVERS",
             hjust = 0.5, vjust = 1, size = 3.2, 
             color = colors$text, fontface = "bold") +
    # Sun icon for radiation (moved up, smaller)
    annotate("point", x = 0.5, y = 0.88, size = 6, 
             color = colors$radiation, fill = colors$radiation, shape = 21) +
    # Sun rays (shorter)
    annotate("segment", x = 0.5, xend = 0.47, y = 0.88, yend = 0.84,
             color = colors$radiation, linewidth = 0.5) +
    annotate("segment", x = 0.5, xend = 0.53, y = 0.88, yend = 0.84,
             color = colors$radiation, linewidth = 0.5) +
    annotate("segment", x = 0.5, xend = 0.5, y = 0.88, yend = 0.82,
             color = colors$radiation, linewidth = 0.5) +
    annotate("segment", x = 0.5, xend = 0.5, y = 0.88, yend = 0.94,
             color = colors$radiation, linewidth = 0.5) +
    # Bars (moved right to avoid text overlap)
    geom_rect(aes(xmin = 0.30, xmax = 0.30 + importance/100 * 0.65,
                  ymin = y_pos - 0.05, ymax = y_pos + 0.02),
              fill = vars_df$color, color = colors$text, linewidth = 0.5) +
    # Labels (moved right to prevent cropping, using left alignment)
    annotate("text", x = 0.02, y = vars_df$y_pos, 
             label = vars_df$variable,
             hjust = 0, vjust = 0.5, size = 2.0, 
             color = colors$text, fontface = "bold") +
    # Percentages (moved right to align with new bar position)
    annotate("text", x = 0.93, y = vars_df$y_pos,
             label = paste0(vars_df$importance, "%"),
             hjust = 0, vjust = 0.5, size = 2.0, 
             color = colors$text, fontface = "bold") +
    # Summary box (moved down)
    annotate("rect", xmin = 0.1, xmax = 0.9, ymin = 0.08, ymax = 0.20,
             fill = "#FFF9E6", color = colors$radiation, linewidth = 0.5,
             linetype = "solid") +
    annotate("text", x = 0.5, y = 0.14, 
             label = "Primary Driver:\nSolar Radiation",
             hjust = 0.5, vjust = 0.5, size = 2.4, 
             color = colors$text, fontface = "bold") +
    xlim(-0.05, 1) + ylim(0, 1) +  # Extended xlim to show left labels
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(0.05, 0.1, 0.05, 0.1, "inches"))  # Increased left/right margins
  
  return(p)
}

# ============================================================================
# MAIN TITLE
# ============================================================================

create_title <- function() {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Glacier Albedo Decline on Heard Island (2012-2024):\nVIIRS Observations and Climate Attribution",
             hjust = 0.5, vjust = 0.5, size = 3.5,
             color = colors$text, fontface = "bold") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(0.02, 0.05, 0.02, 0.05, "inches"))
  
  return(p)
}

# ============================================================================
# ASSEMBLE GRAPHIC ABSTRACT
# ============================================================================

cat("Creating graphic abstract...\n")

# Create panels
left_panel <- create_left_panel()
center_panel <- create_center_panel()
right_panel <- create_right_panel()
title_panel <- create_title()

# Dimensions: 1062 x 2656 pixels (height x width) - doubled size
# At 300 DPI: 1062/300 = 3.54 inches height, 2656/300 = 8.86 inches width
DPI <- 300
width_inches <- 2656 / DPI  # 8.86 inches
height_inches <- 1062 / DPI   # 3.54 inches

# Create main plot with title and three panels
# Use grid.arrange for layout with better spacing
final_plot <- arrangeGrob(
  title_panel,
  arrangeGrob(
    left_panel,
    center_panel,
    right_panel,
    ncol = 3,
    widths = c(1.2, 1.8, 1.2),
    padding = unit(0.05, "inches")
  ),
  nrow = 2,
  heights = c(0.12, 0.88),
  padding = unit(0.03, "inches")
)

# Save as PDF (vector format) only
output_path_pdf <- "figures/graphic_abstract.pdf"
ggsave(output_path_pdf, final_plot, 
       width = width_inches, height = height_inches, 
       units = "in", dpi = DPI, device = "pdf")
cat("✓ Graphic abstract (PDF) saved to:", output_path_pdf, "\n")

cat("\nGraphic abstract created successfully!\n")
cat("Dimensions:", 1062, "x", 2656, "pixels\n")
cat("Format: PDF (vector)\n")

