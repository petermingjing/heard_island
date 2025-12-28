# Heard Island Glacier Albedo Figure 3: Time Series
# Script for generating Figure 3 showing temporal trends
# Created: 2024
# Purpose: Generate time series plot with 31-point running average and linear regression trend

# Load required libraries
library(tidyverse)
library(lubridate)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(scales)
library(extrafont)
library(zoo)  # For rolling average calculations
library(Kendall)  # For Mann-Kendall test

# Set working directory
setwd("")

# Create output directories
if (!dir.exists("output")) {
  dir.create("output")
}
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Load Times New Roman font with error handling
tryCatch({
  loadfonts(device = "pdf", quiet = TRUE)
  cat("✓ Times New Roman font loaded successfully\n")
}, error = function(e) {
  cat("Warning: Could not load Times New Roman font:", e$message, "\n")
  cat("Using system default font instead\n")
})

# Define Times New Roman theme with fallback (same settings as Figure 2)
theme_times <- function() {
  if ("Times New Roman" %in% fonts()) {
    font_family <- "Times New Roman"
  } else {
    font_family <- "serif"
    cat("Warning: Times New Roman not available, using serif font\n")
  }
  
  theme_minimal() +
    theme(
      text = element_text(family = font_family, size = 11),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0)
    )
}

# Load time series data
cat("Loading time series data...\n")
time_series <- read_csv("output/time_series_data.csv", show_col_types = FALSE)
time_series$date <- as.Date(time_series$date)
time_series$year <- as.integer(time_series$year)
time_series$month <- as.integer(time_series$month)

# Filter for VIIRS only
time_series <- time_series %>% filter(sensor == "VIIRS")

# Sort by date for running average calculation
time_series <- time_series %>% arrange(date)

cat("✓ Loaded", nrow(time_series), "VIIRS observations\n")
cat("Date range:", min(time_series$date), "to", max(time_series$date), "\n")

# Calculate 31-point running average
cat("\nCalculating 31-point running average...\n")
time_series$running_avg <- zoo::rollmean(time_series$albedo, k = 31, fill = NA, align = "center")
cat("✓ Running average calculated\n")

# Calculate linear regression trend for annotation
cat("\nCalculating linear regression trend...\n")
annual_means <- time_series %>%
  group_by(year) %>%
  summarise(mean_albedo = mean(albedo, na.rm = TRUE),
            n_obs = n(),
            .groups = "drop") %>%
  filter(n_obs >= 3)  # At least 3 observations per year

if (nrow(annual_means) >= 4) {
  # Linear regression on annual means
  ols_model <- lm(mean_albedo ~ year, data = annual_means)
  ols_slope <- coef(ols_model)[2] * 10  # per decade
  ols_intercept <- coef(ols_model)[1]
  ols_r_squared <- summary(ols_model)$r.squared
  ols_p_value <- summary(ols_model)$coefficients[2, 4]
  
  # Get confidence intervals
  ols_ci <- confint(ols_model, level = 0.95)
  ols_ci_lower <- ols_ci[2, 1] * 10  # per decade
  ols_ci_upper <- ols_ci[2, 2] * 10  # per decade
  
  # Mann-Kendall test (non-parametric)
  mk_result <- MannKendall(annual_means$mean_albedo)
  mk_tau <- mk_result$tau
  mk_p_value <- mk_result$sl
  
  cat(sprintf("Linear trend: %.3f per decade [95%% CI: %.3f, %.3f]\n", 
              ols_slope, ols_ci_lower, ols_ci_upper))
  cat(sprintf("R² = %.3f, p = %.3f\n", ols_r_squared, ols_p_value))
  cat(sprintf("Mann-Kendall: τ = %.3f, p = %.3f\n", mk_tau, mk_p_value))
} else {
  cat("Warning: Insufficient data for trend calculation\n")
  ols_slope <- NA
  ols_r_squared <- NA
  ols_p_value <- NA
  ols_ci_lower <- NA
  ols_ci_upper <- NA
  mk_tau <- NA
  mk_p_value <- NA
}

# Create Figure 3: Time series plot
cat("\n=== Creating Figure 3 ===\n")

# Create annotation text with key statistics
if (!is.na(ols_slope)) {
  annotation_text <- paste0(
    "Trend: ", sprintf("%.3f", ols_slope), " per decade\n",
    "95% CI: [", sprintf("%.3f", ols_ci_lower), ", ", sprintf("%.3f", ols_ci_upper), "]\n",
    "Mann-Kendall: τ = ", sprintf("%.3f", mk_tau), ", p = ", sprintf("%.3f", mk_p_value)
  )
} else {
  annotation_text <- ""
}

# Get date and y-axis ranges for annotation positioning
max_date <- max(time_series$date)
max_y <- 0.8  # Using the y-axis limit

# Create the plot with individual observations, 31-point running average, and linear trend
figure_3 <- ggplot(time_series, aes(x = date, y = albedo)) +
  # Individual observations as points (enlarged)
  geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
  # 31-point running average line
  geom_line(aes(y = running_avg), color = "darkblue", linewidth = 1.2, na.rm = TRUE) +
  # Linear regression 95% confidence interval ribbon (shaded area)
  geom_ribbon(stat = "smooth", method = "lm", se = TRUE,
              fill = "coral", alpha = 0.3, color = NA) +
  # Linear regression trend line
  geom_smooth(method = "lm", se = FALSE,
              color = "red", linewidth = 1.0, linetype = "dashed") +
  # Add statistical annotation (using actual date values)
  annotate("text", x = max_date, y = max_y, 
           label = annotation_text,
           hjust = 1.05, vjust = 1.2,
           size = 3.5, family = ifelse("Times New Roman" %in% fonts(), "Times New Roman", "serif"),
           color = "black") +
  labs(
    x = "Date",
    y = "Albedo"
  ) +
  theme_times() +
  theme(
    plot.margin = ggplot2::margin(15, 15, 15, 15, unit = "pt")
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
  scale_y_continuous(limits = c(0, 0.8))

# Save Figure 3 as PDF only (harmonized size)
ggsave("figures/figure_3.pdf", 
       figure_3, width = 10, height = 6, device = "pdf")

cat("✓ Saved: figure_3.pdf\n")

cat("\n=== Script completed successfully ===\n")

