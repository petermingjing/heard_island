# Heard Island Glacier Albedo Descriptive Statistics
# Script for generating descriptive statistics and box plots for Section 3.1
# Created: 2024
# Purpose: Generate comprehensive descriptive statistics and visualizations
#          for albedo distributions on annual, seasonal, and monthly basis

# Load required libraries
library(tidyverse)
library(lubridate)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(scales)
library(extrafont)

# Set working directory
setwd("/Users/jingming/Desktop/HIG")

# Create output directories
if (!dir.exists("output")) {
  dir.create("output")
}
if (!dir.exists("figures")) {
  dir.create("figures")
}
if (!dir.exists("figures/descriptive_statistics")) {
  dir.create("figures/descriptive_statistics")
}

# Load Times New Roman font with error handling
tryCatch({
  loadfonts(device = "pdf", quiet = TRUE)
  cat("✓ Times New Roman font loaded successfully\n")
}, error = function(e) {
  cat("Warning: Could not load Times New Roman font:", e$message, "\n")
  cat("Using system default font instead\n")
})

# Define Times New Roman theme with fallback
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

cat("✓ Loaded", nrow(time_series), "VIIRS observations\n")
cat("Date range:", min(time_series$date), "to", max(time_series$date), "\n")

# Calculate descriptive statistics
cat("\n=== Calculating Descriptive Statistics ===\n")

# 1. Overall statistics
  overall_stats <- time_series %>%
  summarise(
    n_obs = n(),
    mean_albedo = mean(albedo, na.rm = TRUE),
    median_albedo = median(albedo, na.rm = TRUE),
    sd_albedo = sd(albedo, na.rm = TRUE),
    min_albedo = min(albedo, na.rm = TRUE),
    max_albedo = max(albedo, na.rm = TRUE),
    q25_albedo = quantile(albedo, 0.25, na.rm = TRUE),
    q75_albedo = quantile(albedo, 0.75, na.rm = TRUE),
    iqr_albedo = IQR(albedo, na.rm = TRUE),
    cv_albedo = sd(albedo, na.rm = TRUE) / mean(albedo, na.rm = TRUE) * 100
  ) %>%
  mutate(stat_type = "Overall")

cat("Overall statistics:\n")
print(overall_stats)

# 2. Annual statistics
annual_stats <- time_series %>%
  group_by(year) %>%
  summarise(
    n_obs = n(),
    mean_albedo = mean(albedo, na.rm = TRUE),
    median_albedo = median(albedo, na.rm = TRUE),
    sd_albedo = sd(albedo, na.rm = TRUE),
    min_albedo = min(albedo, na.rm = TRUE),
    max_albedo = max(albedo, na.rm = TRUE),
    q25_albedo = quantile(albedo, 0.25, na.rm = TRUE),
    q75_albedo = quantile(albedo, 0.75, na.rm = TRUE),
    iqr_albedo = IQR(albedo, na.rm = TRUE),
    cv_albedo = sd(albedo, na.rm = TRUE) / mean(albedo, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\nAnnual statistics:\n")
print(annual_stats)

# 3. Seasonal statistics
seasonal_stats <- time_series %>%
  group_by(season) %>%
  summarise(
    n_obs = n(),
    mean_albedo = mean(albedo, na.rm = TRUE),
    median_albedo = median(albedo, na.rm = TRUE),
    sd_albedo = sd(albedo, na.rm = TRUE),
    min_albedo = min(albedo, na.rm = TRUE),
    max_albedo = max(albedo, na.rm = TRUE),
    q25_albedo = quantile(albedo, 0.25, na.rm = TRUE),
    q75_albedo = quantile(albedo, 0.75, na.rm = TRUE),
    iqr_albedo = IQR(albedo, na.rm = TRUE),
    cv_albedo = sd(albedo, na.rm = TRUE) / mean(albedo, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(match(season, c("Summer", "Autumn", "Winter", "Spring")))

cat("\nSeasonal statistics:\n")
print(seasonal_stats)

# 4. Monthly statistics
monthly_stats <- time_series %>%
  group_by(month) %>%
  summarise(
    n_obs = n(),
    mean_albedo = mean(albedo, na.rm = TRUE),
    median_albedo = median(albedo, na.rm = TRUE),
    sd_albedo = sd(albedo, na.rm = TRUE),
    min_albedo = min(albedo, na.rm = TRUE),
    max_albedo = max(albedo, na.rm = TRUE),
    q25_albedo = quantile(albedo, 0.25, na.rm = TRUE),
    q75_albedo = quantile(albedo, 0.75, na.rm = TRUE),
    iqr_albedo = IQR(albedo, na.rm = TRUE),
    cv_albedo = sd(albedo, na.rm = TRUE) / mean(albedo, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(month_name = month.abb[month])

cat("\nMonthly statistics:\n")
print(monthly_stats)

# Save statistics to CSV files
write_csv(overall_stats, "output/descriptive_statistics_overall.csv")
write_csv(annual_stats, "output/descriptive_statistics_annual.csv")
write_csv(seasonal_stats, "output/descriptive_statistics_seasonal.csv")
write_csv(monthly_stats, "output/descriptive_statistics_monthly.csv")

cat("\n✓ Statistics saved to output files\n")

# Create box plots
cat("\n=== Creating Box Plots ===\n")

# Prepare data for plotting
time_series <- time_series %>%
  mutate(
    year_factor = factor(year),
    season_factor = factor(season, levels = c("Summer", "Autumn", "Winter", "Spring")),
    month_name = factor(month.abb[month], levels = month.abb)
  )

# 1. Overall distribution box plot
p1 <- ggplot(time_series, aes(x = "", y = albedo)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
  labs(
    title = "Overall Albedo Distribution",
    subtitle = paste("VIIRS (2012-2024), n =", nrow(time_series), "observations"),
    x = "",
    y = "Albedo"
  ) +
  theme_times() +
  theme(axis.text.x = element_blank())

ggsave("figures/descriptive_statistics/boxplot_overall.pdf", 
       p1, width = 6, height = 6, device = "pdf")

cat("✓ Saved: boxplot_overall.pdf\n")

# 2. Annual box plot
p2 <- ggplot(time_series, aes(x = year_factor, y = albedo)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "red") +
  labs(
    title = "Annual Albedo Distribution",
    subtitle = "VIIRS observations (2012-2024)",
    x = "Year",
    y = "Albedo"
  ) +
  theme_times() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/descriptive_statistics/boxplot_annual.pdf", 
       p2, width = 10, height = 6, device = "pdf")

cat("✓ Saved: boxplot_annual.pdf\n")

# 3. Seasonal box plot
p3 <- ggplot(time_series, aes(x = season_factor, y = albedo)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
  labs(
    title = "Seasonal Albedo Distribution",
    subtitle = paste("VIIRS observations (2012-2024), n =", nrow(time_series)),
    x = "Season",
    y = "Albedo"
  ) +
  theme_times()

ggsave("figures/descriptive_statistics/boxplot_seasonal.pdf", 
       p3, width = 8, height = 6, device = "pdf")

cat("✓ Saved: boxplot_seasonal.pdf\n")

# 4. Monthly box plot
p4 <- ggplot(time_series, aes(x = month_name, y = albedo)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "red") +
  labs(
    title = "Monthly Albedo Distribution",
    subtitle = paste("VIIRS observations (2012-2024), n =", nrow(time_series)),
    x = "Month",
    y = "Albedo"
  ) +
  theme_times() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/descriptive_statistics/boxplot_monthly.pdf", 
       p4, width = 10, height = 6, device = "pdf")

cat("✓ Saved: boxplot_monthly.pdf\n")

# 5. Combined annual and seasonal plot
p5 <- grid.arrange(
  p2 + theme(plot.title = element_text(size = 11)),
  p3 + theme(plot.title = element_text(size = 11)),
  ncol = 2
)

ggsave("figures/descriptive_statistics/boxplot_combined.pdf", 
       p5, width = 14, height = 6, device = "pdf")

cat("✓ Saved: boxplot_combined.pdf\n")

# 6. Create Figure 2 with panels a, b, c (Overall, Seasonal, Annual)
# Create panel versions without titles, with panel labels
# Note: gridExtra is already loaded at the top

# Panel a: Overall distribution
p1_panel <- ggplot(time_series, aes(x = "", y = albedo)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.alpha = 0.3, width = 0.5, linewidth = 1) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
  labs(x = "", y = "Albedo") +
  theme_times() +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt")
  ) +
  annotate("text", x = 0.5, y = Inf, label = "(a)", vjust = 1.5, hjust = 0.5, 
           size = 4, fontface = "bold")

# Panel b: Seasonal distribution
p3_panel <- ggplot(time_series, aes(x = season_factor, y = albedo)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.alpha = 0.3, linewidth = 1) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
  labs(x = "Season", y = "Albedo") +
  theme_times() +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt")
  ) +
  annotate("text", x = 0.5, y = Inf, label = "(b)", vjust = 1.5, hjust = 0, 
           size = 4, fontface = "bold")

# Panel c: Annual distribution
p2_panel <- ggplot(time_series, aes(x = year_factor, y = albedo)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.alpha = 0.3, linewidth = 1) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "red") +
  labs(x = "Year", y = "Albedo") +
  theme_times() +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt")
  ) +
  annotate("text", x = 0.5, y = Inf, label = "(c)", vjust = 1.5, hjust = 0, 
           size = 4, fontface = "bold")

# Create Figure 2: Arrange panels in a single row (1x3 layout)
# Adjust panel widths: overall (narrower), seasonal (medium), annual (wider)
figure_2 <- grid.arrange(
  p1_panel, p3_panel, p2_panel,
  ncol = 3,
  widths = c(1, 1.2, 2.2)  # Adjust relative widths
)

ggsave("figures/figure_2.pdf", 
       figure_2, width = 16, height = 5.5, device = "pdf")

cat("✓ Saved: figure_2.pdf\n")

# Create comprehensive summary table
summary_table <- bind_rows(
  overall_stats %>% mutate(group = "Overall", group_value = "All"),
  annual_stats %>% mutate(group = "Annual", group_value = as.character(year)),
  seasonal_stats %>% mutate(group = "Seasonal", group_value = season),
  monthly_stats %>% mutate(group = "Monthly", group_value = month_name)
) %>%
  dplyr::select(group, group_value, n_obs, mean_albedo, median_albedo, sd_albedo, 
         min_albedo, max_albedo, q25_albedo, q75_albedo)

write_csv(summary_table, "output/descriptive_statistics_summary.csv")

cat("\n✓ Summary table saved\n")

# Print summary statistics for manuscript
cat("\n=== Summary Statistics for Manuscript ===\n")
cat("\nOverall:\n")
cat(sprintf("  Mean: %.3f ± %.3f (SD)\n", overall_stats$mean_albedo, overall_stats$sd_albedo))
cat(sprintf("  Median: %.3f\n", overall_stats$median_albedo))
cat(sprintf("  Range: %.3f - %.3f\n", overall_stats$min_albedo, overall_stats$max_albedo))
cat(sprintf("  IQR: %.3f\n", overall_stats$iqr_albedo))

cat("\nBy Season:\n")
for (i in 1:nrow(seasonal_stats)) {
  cat(sprintf("  %s: %.3f ± %.3f (mean ± SD), n = %d\n", 
              seasonal_stats$season[i], 
              seasonal_stats$mean_albedo[i], 
              seasonal_stats$sd_albedo[i],
              seasonal_stats$n_obs[i]))
}

cat("\nBy Year:\n")
for (i in 1:nrow(annual_stats)) {
  cat(sprintf("  %d: %.3f ± %.3f (mean ± SD), n = %d\n", 
              annual_stats$year[i], 
              annual_stats$mean_albedo[i], 
              annual_stats$sd_albedo[i],
              annual_stats$n_obs[i]))
}

cat("\n=== Script completed successfully ===\n")

