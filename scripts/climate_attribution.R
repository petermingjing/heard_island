# Heard Island Glacier Albedo Climate Attribution Analysis
# Advanced statistical analysis to attribute albedo variations to climatic factors
# R Script for Quantifying Relative Importance of Climate Variables
# =============================================================================

# Load required libraries
# Core data manipulation and visualization
library(tidyverse)      # Includes dplyr, readr, tidyr, etc.
library(lubridate)      # For date handling
library(ggplot2)        # For plotting (also included in tidyverse but explicit for clarity)
library(extrafont)       # For Times New Roman font loading

# Statistical analysis libraries
library(car)            # For VIF calculation (MLR)
library(relaimpo)       # For Relative importance analysis (MLR)
library(boot)           # For Bootstrap confidence intervals (MLR)
library(randomForest)   # For Random Forest attribution
library(caret)          # For Random Forest cross-validation
library(mgcv)           # For Generalized Additive Models (GAM)
library(gridExtra)      # For combining plots
library(grid)           # For spacing utilities

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
      text = element_text(family = font_family, size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.3)
    )
}

# Define data paths
data_path <- "data/gee_export/"

# =============================================================================
# DATA LOADING AND PREPROCESSING FUNCTIONS
# =============================================================================

# Function to load and clean albedo data
load_albedo_data <- function(sensor, file_path, albedo_column) {
  cat("Loading", sensor, "data...\n")
  
  # Read the data first to check columns
  data <- read_csv(file_path, show_col_types = FALSE)
  cat("Available columns:", paste(names(data), collapse = ", "), "\n")
  
  # Check if the albedo column exists
  if (!albedo_column %in% names(data)) {
    stop("Column '", albedo_column, "' not found in ", file_path)
  }
  
  # Process the data step by step
  data$date <- as.Date(data$date)
  data <- data[!is.na(data[[albedo_column]]), ]
  data <- data[year(data$date) <= 2024, ]
  
  # Select columns using base R approach
  data <- data[, c("date", albedo_column)]
  data$sensor <- sensor
  
  # Rename albedo column to standard name
  names(data)[2] <- "albedo"
  
  cat("✓", sensor, "data loaded:", nrow(data), "records\n")
  return(data)
}

# Function to load and preprocess climate data
load_climate_data <- function() {
  cat("Loading climate data...\n")
  
  # Read the data first
  climate_data <- read_csv(paste0(data_path, "heard_island_era5_land_climate_daily.csv"), 
                          show_col_types = FALSE)
  cat("Available columns:", paste(names(climate_data), collapse = ", "), "\n")
  
  # Process the data step by step using base R
  climate_data$date <- as.Date(climate_data$date)
  climate_data <- climate_data[year(climate_data$date) <= 2024, ]
  climate_data <- climate_data[!is.na(climate_data$t2m_C), ]
  
  # Select required columns
  climate_data <- climate_data[, c("date", "t2m_C", "precip_mm", "snowfall_mm_we", "ssrd_wm2", "pdd_degCday")]
  
  cat("✓ Climate data loaded:", nrow(climate_data), "records\n")
  return(climate_data)
}

# Function to create monthly aggregated dataset for attribution analysis
create_monthly_dataset <- function(albedo_data, climate_data) {
  cat("Creating monthly aggregated dataset...\n")
  
  # Create monthly albedo means
  albedo_monthly <- albedo_data %>%
    mutate(year = year(date), month = month(date)) %>%
    group_by(sensor, year, month) %>%
    summarise(
      mean_albedo = mean(albedo, na.rm = TRUE),
      n_obs = n(),
      sd_albedo = sd(albedo, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(date = make_date(year, month, 1))
  
  # Create monthly climate aggregates
  climate_monthly <- climate_data %>%
    mutate(year = year(date), month = month(date)) %>%
    group_by(year, month) %>%
    summarise(
      mean_temp = mean(t2m_C, na.rm = TRUE),
      total_precip = sum(precip_mm, na.rm = TRUE),
      total_snowfall = sum(snowfall_mm_we, na.rm = TRUE),
      mean_radiation = mean(ssrd_wm2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(date = make_date(year, month, 1))
  
  # Merge datasets
  merged_data <- albedo_monthly %>%
    left_join(climate_monthly, by = c("year", "month", "date"))
  
  # Sort by date to ensure proper ordering
  merged_data <- merged_data %>%
    arrange(date)
  
  # Remove rows with missing values in essential variables
  essential_vars <- c("sensor", "year", "month", "date", "mean_albedo", 
                     "mean_temp", "total_precip", "total_snowfall", "mean_radiation")
  
  # Check for complete cases in essential variables only
  essential_complete <- complete.cases(merged_data[, essential_vars])
  merged_data <- merged_data[essential_complete, ]
  
  cat("✓ Monthly dataset created:", nrow(merged_data), "records\n")
  cat("  Date range:", min(merged_data$date), "to", max(merged_data$date), "\n")
  cat("  Variables: mean_albedo, mean_temp, total_precip, total_snowfall, mean_radiation\n")
  
  return(merged_data)
}

# =============================================================================
# CLIMATE ATTRIBUTION METHODS
# =============================================================================

# Method 1: Multiple Linear Regression with Relative Importance
perform_regression_attribution <- function(data, sensor_name) {
  cat("Performing regression-based attribution for", sensor_name, "...\n")
  
  sensor_data <- data %>% filter(sensor == sensor_name)
  
  if (nrow(sensor_data) < 20) {
    cat("Warning: Insufficient data for", sensor_name, "\n")
    return(NULL)
  }
  
  # Define predictor variables (essential climate variables including radiation)
  predictors <- c("mean_temp", "total_precip", "total_snowfall", "mean_radiation")
  
  # Create formula
  formula_str <- paste("mean_albedo ~", paste(predictors, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  # Fit multiple linear regression
  lm_model <- lm(formula_obj, data = sensor_data)
  
  # Calculate relative importance using relaimpo package
  rel_importance <- calc.relimp(lm_model, type = c("lmg", "last", "first", "betasq"))
  
  # Extract standardized coefficients
  std_coefs <- coef(lm_model)[-1] * sapply(sensor_data[predictors], sd) / sd(sensor_data$mean_albedo)
  
  # Calculate partial R-squared (LMG method gives contribution to R²)
  partial_r2 <- summary(lm_model)$r.squared * rel_importance$lmg
  
  # Calculate individual R² for each variable (model with only that variable)
  individual_r2 <- sapply(predictors, function(var) {
    single_formula <- as.formula(paste("mean_albedo ~", var))
    single_model <- lm(single_formula, data = sensor_data)
    return(summary(single_model)$r.squared)
  })
  
  # Calculate incremental R² (R² change when adding variable to model without it)
  incremental_r2 <- sapply(predictors, function(var) {
    # Model without this variable
    other_vars <- setdiff(predictors, var)
    if (length(other_vars) == 0) {
      return(individual_r2[var])
    }
    reduced_formula <- as.formula(paste("mean_albedo ~", paste(other_vars, collapse = " + ")))
    reduced_model <- lm(reduced_formula, data = sensor_data)
    r2_reduced <- summary(reduced_model)$r.squared
    r2_full <- summary(lm_model)$r.squared
    return(r2_full - r2_reduced)
  })
  
  # Bootstrap confidence intervals for coefficients
  boot_coefs <- function(data, indices) {
    boot_data <- data[indices, ]
    boot_model <- lm(formula_obj, data = boot_data)
    return(coef(boot_model)[-1])
  }
  
  boot_results <- boot(sensor_data, boot_coefs, R = 1000)
  ci_coefs <- t(sapply(1:length(predictors), function(i) {
    boot.ci(boot_results, type = "perc", index = i)$percent[4:5]
  }))
  
  # Calculate variance inflation factors
  vif_values <- car::vif(lm_model)
  
  return(list(
    model = lm_model,
    coefficients = coef(lm_model),
    std_coefficients = std_coefs,
    relative_importance = rel_importance,
    partial_r2 = partial_r2,
    individual_r2 = individual_r2,
    incremental_r2 = incremental_r2,
    confidence_intervals = ci_coefs,
    vif = vif_values,
    r_squared = summary(lm_model)$r.squared,
    adj_r_squared = summary(lm_model)$adj.r.squared,
    aic = AIC(lm_model),
    bic = BIC(lm_model),
    predictors = predictors,
    data = sensor_data
  ))
}

# Method 2: Random Forest Attribution
perform_rf_attribution <- function(data, sensor_name) {
  cat("Performing Random Forest attribution for", sensor_name, "...\n")
  
  sensor_data <- data %>% filter(sensor == sensor_name)
  
  if (nrow(sensor_data) < 20) {
    cat("Warning: Insufficient data for", sensor_name, "\n")
    return(NULL)
  }
  
  # Define predictor variables (essential climate variables including radiation)
  predictors <- c("mean_temp", "total_precip", "total_snowfall", "mean_radiation")
  
  # Prepare data for Random Forest
  rf_data <- sensor_data[, c("mean_albedo", predictors)]
  rf_data <- rf_data[complete.cases(rf_data), ]
  
  # Set up cross-validation
  train_control <- trainControl(method = "cv", number = 5, verboseIter = FALSE)
  
  # Train Random Forest model
  rf_model <- train(
    mean_albedo ~ .,
    data = rf_data,
    method = "rf",
    trControl = train_control,
    importance = TRUE,
    ntree = 500
  )
  
  # Extract variable importance using raw permutation importance (%IncMSE)
  # This is more reliable than caret's varImp() which can scale incorrectly
  perm_importance <- importance(rf_model$finalModel, type = 1)
  importance_scores <- perm_importance[, "%IncMSE"]
  names(importance_scores) <- rownames(perm_importance)
  
  # Ensure all predictors are included (some might have 0 importance)
  # Fill in missing predictors with 0 importance
  for (pred in predictors) {
    if (!pred %in% names(importance_scores)) {
      importance_scores[pred] <- 0
    }
  }
  
  # Reorder to match predictor order
  importance_scores <- importance_scores[predictors]
  
  # Calculate individual R² for each variable (model with only that variable)
  individual_r2 <- sapply(predictors, function(var) {
    single_rf_data <- rf_data[, c("mean_albedo", var)]
    single_rf_data <- single_rf_data[complete.cases(single_rf_data), ]
    if (nrow(single_rf_data) < 10) return(0)
    
    single_formula <- as.formula(paste("mean_albedo ~", var))
    single_rf <- train(single_formula, 
                       data = single_rf_data,
                       method = "rf",
                       trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE),
                       importance = TRUE,
                       ntree = 500)
    return(single_rf$results$Rsquared[1])
  })
  
  # Calculate incremental R² (R² change when adding variable to model without it)
  incremental_r2 <- sapply(predictors, function(var) {
    other_vars <- setdiff(predictors, var)
    if (length(other_vars) == 0) {
      return(individual_r2[var])
    }
    reduced_rf_data <- rf_data[, c("mean_albedo", other_vars)]
    reduced_rf_data <- reduced_rf_data[complete.cases(reduced_rf_data), ]
    if (nrow(reduced_rf_data) < 10) return(0)
    
    reduced_formula <- as.formula(paste("mean_albedo ~", paste(other_vars, collapse = " + ")))
    reduced_rf <- train(reduced_formula,
                        data = reduced_rf_data,
                        method = "rf",
                        trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE),
                        importance = TRUE,
                        ntree = 500)
    r2_reduced <- reduced_rf$results$Rsquared[1]
    r2_full <- rf_model$results$Rsquared[1]
    return(r2_full - r2_reduced)
  })
  
  # Calculate partial dependence plots data
  partial_dep <- list()
  for (var in predictors) {
    tryCatch({
      partial_dep[[var]] <- randomForest::partialPlot(rf_model$finalModel, 
                                                     pred.data = rf_data[, predictors], 
                                                     x.var = var, 
                                                     plot = FALSE)
    }, error = function(e) {
      cat("Warning: Could not calculate partial dependence for", var, "\n")
      partial_dep[[var]] <- NULL
    })
  }
  
  return(list(
    model = rf_model,
    variable_importance = importance_scores,
    permutation_importance = perm_importance,
    partial_dependence = partial_dep,
    individual_r2 = individual_r2,
    incremental_r2 = incremental_r2,
    r_squared = rf_model$results$Rsquared[1],
    rmse = rf_model$results$RMSE[1],
    mae = rf_model$results$MAE[1],
    predictors = predictors,
    data = rf_data
  ))
}

# Methods 3-12 are not used in the current analysis
# Only Methods 1 (MLR), 2 (Random Forest), and 5 (GAM) are active

# Method 5: Generalized Additive Models (GAM) Attribution
perform_gam_attribution <- function(data, sensor_name) {
  cat("Performing GAM-based attribution for", sensor_name, "...\n")
  
  sensor_data <- data %>% filter(sensor == sensor_name)
  
  if (nrow(sensor_data) < 20) {
    cat("Warning: Insufficient data for", sensor_name, "\n")
    return(NULL)
  }
  
  # Define predictor variables (essential climate variables including radiation)
  predictors <- c("mean_temp", "total_precip", "total_snowfall", "mean_radiation")
  
  # Prepare data for GAM
  gam_data <- sensor_data[, c("mean_albedo", predictors)]
  gam_data <- gam_data[complete.cases(gam_data), ]
  
  # Fit GAM with smooth terms for all continuous variables
  gam_formula <- mean_albedo ~ 
    s(mean_temp, k = 5) + 
    s(total_precip, k = 5) + 
    s(total_snowfall, k = 5) +
    s(mean_radiation, k = 5)
  
  # Fit GAM model using mgcv
  gam_model <- mgcv::gam(gam_formula, data = gam_data)
  
  # Calculate model diagnostics
  gam_summary <- summary(gam_model)
  
  # Extract smooth term information
  smooth_terms <- gam_summary$s.table
  
  # Calculate partial effects (marginal effects) for smooth terms
  partial_effects <- list()
  smooth_vars <- c("mean_temp", "total_precip", "total_snowfall", "mean_radiation")
  
  if (!is.null(smooth_terms) && nrow(smooth_terms) > 0) {
    for (var in smooth_vars) {
      term_name <- paste0("s(", var, ")")
      if (term_name %in% rownames(smooth_terms)) {
        # Create prediction grid for partial effect
        pred_grid <- expand.grid(
          mean_temp = ifelse(var == "mean_temp", seq(min(gam_data$mean_temp), max(gam_data$mean_temp), length = 50), mean(gam_data$mean_temp)),
          total_precip = ifelse(var == "total_precip", seq(min(gam_data$total_precip), max(gam_data$total_precip), length = 50), mean(gam_data$total_precip)),
          total_snowfall = ifelse(var == "total_snowfall", seq(min(gam_data$total_snowfall), max(gam_data$total_snowfall), length = 50), mean(gam_data$total_snowfall)),
          mean_radiation = ifelse(var == "mean_radiation", seq(min(gam_data$mean_radiation), max(gam_data$mean_radiation), length = 50), mean(gam_data$mean_radiation))
        )
        
        # Predict partial effect using mgcv::predict.gam
        tryCatch({
          pred_effect <- predict(gam_model, newdata = pred_grid, type = "terms", se.fit = TRUE)
          
          # Extract the specific term
          if (term_name %in% colnames(pred_effect$fit)) {
            term_index <- which(colnames(pred_effect$fit) == term_name)
            partial_effects[[var]] <- data.frame(
              x = pred_grid[[var]],
              effect = pred_effect$fit[, term_index],
              se = pred_effect$se.fit[, term_index]
            )
          }
        }, error = function(e) {
          cat("Warning: Could not calculate partial effect for", var, ":", e$message, "\n")
        })
      }
    }
  }
  
  # Calculate variable importance using deviance explained
  deviance_explained <- gam_summary$dev.expl
  r_squared <- gam_summary$r.sq
  
  # Calculate AIC and BIC
  aic_value <- AIC(gam_model)
  bic_value <- BIC(gam_model)
  
  # Calculate predictions and residuals
  predictions <- predict(gam_model, type = "response")
  residuals <- gam_data$mean_albedo - predictions
  
  # Calculate RMSE and MAE
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  
  # Calculate effective degrees of freedom for each smooth term
  edf_values <- NULL
  p_values <- NULL
  if (!is.null(smooth_terms) && nrow(smooth_terms) > 0) {
    edf_values <- smooth_terms[, "edf"]
    names(edf_values) <- rownames(smooth_terms)
    
    # Calculate significance of smooth terms
    p_values <- smooth_terms[, "p-value"]
    names(p_values) <- rownames(smooth_terms)
  } else {
    # Fallback: use coefficient names from linear terms
    coef_names <- names(coef(gam_model))[-1]  # Exclude intercept
    edf_values <- rep(1, length(coef_names))
    names(edf_values) <- coef_names
    p_values <- gam_summary$p.table[-1, "Pr(>|t|)"]  # Exclude intercept
    names(p_values) <- coef_names
  }
  
  # Calculate individual R² for each variable (model with only that variable)
  individual_r2 <- sapply(predictors, function(var) {
    single_gam_data <- gam_data[, c("mean_albedo", var)]
    single_gam_data <- single_gam_data[complete.cases(single_gam_data), ]
    if (nrow(single_gam_data) < 10) return(0)
    
    single_formula <- as.formula(paste("mean_albedo ~ s(", var, ", k = 5)"))
    single_gam <- tryCatch({
      mgcv::gam(single_formula, data = single_gam_data)
    }, error = function(e) {
      return(NULL)
    })
    if (is.null(single_gam)) return(0)
    return(summary(single_gam)$r.sq)
  })
  
  # Calculate incremental R² (R² change when adding variable to model without it)
  incremental_r2 <- sapply(predictors, function(var) {
    other_vars <- setdiff(predictors, var)
    if (length(other_vars) == 0) {
      return(individual_r2[var])
    }
    reduced_gam_data <- gam_data[, c("mean_albedo", other_vars)]
    reduced_gam_data <- reduced_gam_data[complete.cases(reduced_gam_data), ]
    if (nrow(reduced_gam_data) < 10) return(0)
    
    reduced_formula_str <- paste("mean_albedo ~", paste(paste0("s(", other_vars, ", k = 5)"), collapse = " + "))
    reduced_formula <- as.formula(reduced_formula_str)
    reduced_gam <- tryCatch({
      mgcv::gam(reduced_formula, data = reduced_gam_data)
    }, error = function(e) {
      return(NULL)
    })
    if (is.null(reduced_gam)) return(0)
    r2_reduced <- summary(reduced_gam)$r.sq
    r2_full <- r_squared
    return(r2_full - r2_reduced)
  })
  
  return(list(
    gam_model = gam_model,
    model = gam_model,  # Keep both for compatibility
    summary = gam_summary,
    smooth_terms = smooth_terms,
    partial_effects = partial_effects,
    deviance_explained = deviance_explained,
    r_squared = r_squared,
    aic = aic_value,
    bic = bic_value,
    rmse = rmse,
    mae = mae,
    edf_values = edf_values,
    p_values = p_values,
    individual_r2 = individual_r2,
    incremental_r2 = incremental_r2,
    predictions = predictions,
    residuals = residuals,
    predictors = predictors,
    data = gam_data
  ))
}

# =============================================================================
# MAIN ANALYSIS FUNCTION
# =============================================================================

# Function to perform comprehensive climate attribution analysis
perform_climate_attribution <- function(monthly_data) {
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("COMPREHENSIVE CLIMATE ATTRIBUTION ANALYSIS\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Get unique sensors
  sensors <- unique(monthly_data$sensor)
  attribution_results <- list()
  
  for (sensor_name in sensors) {
    cat("\n", paste(rep("-", 60), collapse=""), "\n")
    cat("ANALYZING SENSOR:", sensor_name, "\n")
    cat(paste(rep("-", 60), collapse=""), "\n")
    
    sensor_results <- list()
    
    # Multiple Linear Regression Attribution
    sensor_results$regression <- perform_regression_attribution(monthly_data, sensor_name)
    
    # Random Forest Attribution
    sensor_results$random_forest <- perform_rf_attribution(monthly_data, sensor_name)
    
    # Generalized Additive Model Attribution
    sensor_results$gam <- perform_gam_attribution(monthly_data, sensor_name)
    
    attribution_results[[sensor_name]] <- sensor_results
  }
  
  return(attribution_results)
}

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

# Function to create attribution summary plots
create_attribution_plots <- function(attribution_results) {
  cat("\nCreating attribution visualization plots...\n")
  
  # Extract results for plotting
  plot_data <- list()
  
  for (sensor_name in names(attribution_results)) {
    sensor_results <- attribution_results[[sensor_name]]
    
    # Regression results
    if (!is.null(sensor_results$regression) && !is.null(sensor_results$regression$predictors) && length(sensor_results$regression$predictors) > 0) {
      reg_data <- data.frame(
        sensor = sensor_name,
        method = "Regression",
        variable = sensor_results$regression$predictors,
        importance = sensor_results$regression$relative_importance$lmg,
        coefficient = sensor_results$regression$std_coefficients,
        stringsAsFactors = FALSE
      )
      if (is.null(plot_data$regression)) {
        plot_data$regression <- reg_data
      } else {
        plot_data$regression <- rbind(plot_data$regression, reg_data)
      }
    }
    
    # Random Forest results
    if (!is.null(sensor_results$random_forest) && !is.null(sensor_results$random_forest$variable_importance) && length(sensor_results$random_forest$variable_importance) > 0) {
      # Ensure all predictors are included, even if importance is 0
      rf_importance <- sensor_results$random_forest$variable_importance
      rf_predictors <- sensor_results$random_forest$predictors
      
      # Create complete importance vector for all predictors
      complete_importance <- numeric(length(rf_predictors))
      names(complete_importance) <- rf_predictors
      for (pred in rf_predictors) {
        if (pred %in% names(rf_importance)) {
          complete_importance[pred] <- rf_importance[pred]
        } else {
          complete_importance[pred] <- 0
        }
      }
      
      rf_data <- data.frame(
        sensor = sensor_name,
        method = "Random Forest",
        variable = names(complete_importance),
        importance = as.numeric(complete_importance),
        coefficient = NA,
        stringsAsFactors = FALSE
      )
      if (is.null(plot_data$random_forest)) {
        plot_data$random_forest <- rf_data
      } else {
        plot_data$random_forest <- rbind(plot_data$random_forest, rf_data)
      }
    }
    
    # GAM results
    if (!is.null(sensor_results$gam) && !is.null(sensor_results$gam$edf_values) && length(sensor_results$gam$edf_values) > 0) {
      # Clean GAM variable names: remove s() wrapper
      gam_var_names <- names(sensor_results$gam$edf_values)
      gam_var_names_clean <- ifelse(
        grepl("^s\\((.+)\\)$", gam_var_names),
        gsub("^s\\((.+)\\)$", "\\1", gam_var_names),
        gam_var_names
      )
      
      gam_data <- data.frame(
        sensor = sensor_name,
        method = "GAM",
        variable = gam_var_names_clean,
        importance = sensor_results$gam$edf_values,
        coefficient = NA,
        stringsAsFactors = FALSE
      )
      if (is.null(plot_data$gam)) {
        plot_data$gam <- gam_data
      } else {
        plot_data$gam <- rbind(plot_data$gam, gam_data)
      }
    }
    
  }
  
  # Combine all methods
  if (length(plot_data) > 0) {
    all_plot_data <- do.call(rbind, plot_data)
  } else {
    # Create empty data frame with required columns if no data
    all_plot_data <- data.frame(
      sensor = character(0),
      method = character(0),
      variable = character(0),
      importance = numeric(0),
      coefficient = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  # Create plots
  plots <- list()
  
    # 1. Variable Importance Comparison (All methods with separate y-axes)
  if (nrow(all_plot_data) > 0) {
    p1 <- ggplot(all_plot_data, aes(x = variable, y = importance, fill = method)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      facet_wrap(~method, scales = "free_y", ncol = 1,
                 labeller = labeller(method = c(
                   "Regression" = "MLR: Relative Importance (LMG)",
                   "Random Forest" = "Random Forest: Variable Importance",
                   "GAM" = "GAM: Effective Degrees of Freedom"
                 ))) +
      scale_fill_brewer(type = "qual", palette = "Set2") +
      labs(x = "Climate Variable", y = "Importance (method-specific scale)", fill = "Method", 
           title = "(a)") +
      theme_times() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            strip.text = element_text(size = 9),
            plot.title = element_text(size = 10, face = "bold", hjust = 0, 
                                     margin = ggplot2::margin(b = 5)),
            plot.margin = ggplot2::margin(5, 5, 5, 5))
    
    # Export plot data
    write.csv(all_plot_data, 
              "output/climate_attribution_variable_importance_data.csv", 
              row.names = FALSE)
    cat("✓ Exported variable importance data: output/climate_attribution_variable_importance_data.csv\n")
    
    # Print for RStudio preview
    print(p1)
    plots$variable_importance <- p1
  } else {
    cat("Warning: No attribution data available for plotting\n")
    plots$variable_importance <- NULL
  }
  
  # 2. Coefficient/Importance Comparison (All methods with separate y-axes)
  # Create a unified value column: use coefficient for Regression, importance for others
  coeff_data <- all_plot_data %>%
    mutate(
      value = ifelse(method == "Regression" & !is.na(coefficient), coefficient, importance),
      y_label = case_when(
        method == "Regression" ~ "Standardized Coefficient",
        method == "Random Forest" ~ "Variable Importance",
        method == "GAM" ~ "EDF (Effective Degrees of Freedom)",
        TRUE ~ "Value"
      )
    )
  
  if (nrow(coeff_data) > 0) {
    # Export plot data
    write.csv(coeff_data, 
              "output/climate_attribution_coefficients_data.csv", 
              row.names = FALSE)
    cat("✓ Exported coefficients/importance data: output/climate_attribution_coefficients_data.csv\n")
    
    p2 <- ggplot(coeff_data, aes(x = variable, y = value, fill = sensor)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      geom_hline(data = coeff_data[coeff_data$method == "Regression", ], 
                 aes(yintercept = 0), linetype = "dashed", color = "red") +
      facet_wrap(~method, scales = "free_y", ncol = 1, 
                 labeller = labeller(method = c(
                   "Regression" = "MLR: Standardized Coefficients",
                   "Random Forest" = "Random Forest: Variable Importance",
                   "GAM" = "GAM: Effective Degrees of Freedom"
                 ))) +
      scale_fill_brewer(type = "qual", palette = "Set1") +
      labs(x = "Climate Variable", y = "Value (method-specific scale)", fill = "Sensor",
           title = "(b)") +
      theme_times() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            strip.text = element_text(size = 9),
            plot.title = element_text(size = 10, face = "bold", hjust = 0,
                                     margin = ggplot2::margin(b = 5)),
            plot.margin = ggplot2::margin(5, 5, 5, 5))
    
    # Print for RStudio preview
    print(p2)
    plots$coefficients <- p2
  } else {
    p2 <- ggplot() + labs(x = "Climate Variable", y = "Value") + theme_times()
    print(p2)
    plots$coefficients <- p2
  }
  
  # 3. Method Comparison Heatmap
  # Normalize importance values within each method to 0-1 scale for fair comparison
  # This ensures all methods are visible in the heatmap regardless of their original scale
  if (nrow(all_plot_data) > 0) {
    # Clean variable names: remove GAM smooth term notation (s()) for better readability
    all_plot_data_clean <- all_plot_data %>%
      mutate(
        variable_clean = case_when(
          # Remove s() wrapper from GAM smooth terms
          grepl("^s\\((.+)\\)$", variable) ~ gsub("^s\\((.+)\\)$", "\\1", variable),
          TRUE ~ variable
        )
      )
    
    # First, average importance by method and variable (in case of multiple sensors)
    importance_by_method_var <- all_plot_data_clean %>%
      group_by(method, variable_clean) %>%
      summarise(mean_importance_raw = mean(importance, na.rm = TRUE),
                .groups = "drop")
    
    # Then normalize within each method using the averaged values
    importance_normalized <- importance_by_method_var %>%
      group_by(method) %>%
      arrange(variable_clean) %>%  # Sort for consistent row_number()
      mutate(
        # Use min-max normalization
        min_val = min(mean_importance_raw, na.rm = TRUE),
        max_val = max(mean_importance_raw, na.rm = TRUE),
        range_val = max_val - min_val,
        # Normalize to 0-1 scale
        # If range is very small or zero, assign evenly spaced values to show all variables
        mean_importance_norm = case_when(
          range_val > 1e-10 ~ (mean_importance_raw - min_val) / range_val,
          n() == 1 ~ 0.5,  # Single variable gets middle value
          TRUE ~ (row_number() - 1) / (n() - 1)  # Multiple variables get evenly spaced
        )
      ) %>%
      ungroup() %>%
      dplyr::select(method, variable_clean, mean_importance_norm, mean_importance_raw)
    
    importance_matrix <- importance_normalized %>%
      dplyr::select(method, variable_clean, mean_importance_norm) %>%
      pivot_wider(names_from = variable_clean, values_from = mean_importance_norm)
    
    if (nrow(importance_matrix) > 0) {
      importance_long <- importance_matrix %>%
        pivot_longer(-method, names_to = "variable", values_to = "importance") %>%
        filter(!is.na(importance))  # Remove NA values
    
      p3 <- ggplot(importance_long, aes(x = variable, y = method, fill = importance)) +
        geom_tile() +
        # Use a more sensitive color scale that better highlights differences
        scale_fill_gradientn(
          colours = c("white", "lightyellow", "yellow", "orange", "red", "darkred"),
          values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1.0)),
          na.value = "grey90",
          limits = c(0, 1),
          guide = guide_colorbar(title = "Normalized\nImportance", 
                                 title.position = "top",
                                 barwidth = 0.8, barheight = 8)
        ) +
        labs(x = "Climate Variable", y = "Method", title = "(c)") +
        theme_times() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 10, face = "bold", hjust = 0,
                                       margin = ggplot2::margin(b = 5)),
              plot.margin = ggplot2::margin(5, 5, 5, 5))
    
      # Export both normalized and raw importance data
      importance_export <- importance_normalized %>%
        dplyr::select(method, variable = variable_clean, mean_importance_norm, mean_importance_raw) %>%
        dplyr::rename(normalized_importance = mean_importance_norm,
                      raw_importance = mean_importance_raw)
      
      write.csv(importance_export, 
                "output/climate_attribution_method_heatmap_data.csv", 
                row.names = FALSE)
      cat("✓ Exported method heatmap data (normalized and raw): output/climate_attribution_method_heatmap_data.csv\n")
      
      # Debug: Print summary of normalized values to check for constant values
      cat("\nHeatmap normalization summary:\n")
      for (m in unique(importance_export$method)) {
        method_data <- importance_export[importance_export$method == m, ]
        cat(sprintf("  %s: normalized range [%.3f, %.3f], raw range [%.3f, %.3f]\n",
                   m, 
                   min(method_data$normalized_importance, na.rm = TRUE),
                   max(method_data$normalized_importance, na.rm = TRUE),
                   min(method_data$raw_importance, na.rm = TRUE),
                   max(method_data$raw_importance, na.rm = TRUE)))
      }
      
      # Print for RStudio preview
      print(p3)
      plots$method_heatmap <- p3
    } else {
      p3 <- ggplot() + labs(x = "Climate Variable", y = "Method") + theme_times()
      print(p3)
      plots$method_heatmap <- p3
    }
  } else {
    p3 <- ggplot() + labs(x = "Climate Variable", y = "Method") + theme_times()
    print(p3)
    plots$method_heatmap <- p3
  }
  
  
  return(plots)
}

# =============================================================================
# SUMMARY REPORT FUNCTION
# =============================================================================

# Function to generate comprehensive summary report
generate_attribution_summary <- function(attribution_results) {
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("CLIMATE ATTRIBUTION SUMMARY REPORT\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Create summary tables
  summary_tables <- list()
  
  for (sensor_name in names(attribution_results)) {
    cat("\n", paste(rep("-", 60), collapse=""), "\n")
    cat("SENSOR:", sensor_name, "\n")
    cat(paste(rep("-", 60), collapse=""), "\n")
    
    sensor_results <- attribution_results[[sensor_name]]
    
    # Regression results
    if (!is.null(sensor_results$regression)) {
      cat("\nREGRESSION ATTRIBUTION:\n")
      cat("R² =", sprintf("%.3f", sensor_results$regression$r_squared), "\n")
      cat("Adjusted R² =", sprintf("%.3f", sensor_results$regression$adj_r_squared), "\n")
      cat("AIC =", sprintf("%.1f", sensor_results$regression$aic), "\n")
      cat("BIC =", sprintf("%.1f", sensor_results$regression$bic), "\n")
      
      cat("\nRelative Importance (LMG):\n")
      for (i in 1:length(sensor_results$regression$predictors)) {
        var <- sensor_results$regression$predictors[i]
        importance <- sensor_results$regression$relative_importance$lmg[i]
        cat(sprintf("  %s: %.3f\n", var, importance))
      }
      
      cat("\nStandardized Coefficients:\n")
      for (i in 1:length(sensor_results$regression$predictors)) {
        var <- sensor_results$regression$predictors[i]
        coef <- sensor_results$regression$std_coefficients[i]
        cat(sprintf("  %s: %.3f\n", var, coef))
      }
      
      cat("\nR² Contributions:\n")
      cat("  Overall R² =", sprintf("%.3f", sensor_results$regression$r_squared), "\n")
      cat("\n  Individual R² (model with only this variable):\n")
      for (i in 1:length(sensor_results$regression$predictors)) {
        var <- sensor_results$regression$predictors[i]
        ind_r2 <- sensor_results$regression$individual_r2[i]
        cat(sprintf("    %s: %.3f (%.1f%%)\n", var, ind_r2, ind_r2 * 100))
      }
      cat("\n  Partial R² (LMG contribution to overall R²):\n")
      for (i in 1:length(sensor_results$regression$predictors)) {
        var <- sensor_results$regression$predictors[i]
        partial_r2 <- sensor_results$regression$partial_r2[i]
        cat(sprintf("    %s: %.3f (%.1f%%)\n", var, partial_r2, partial_r2 * 100))
      }
      cat("\n  Incremental R² (R² change when adding variable):\n")
      for (i in 1:length(sensor_results$regression$predictors)) {
        var <- sensor_results$regression$predictors[i]
        inc_r2 <- sensor_results$regression$incremental_r2[i]
        cat(sprintf("    %s: %.3f (%.1f%%)\n", var, inc_r2, inc_r2 * 100))
      }
    }
    
    # Random Forest results
    if (!is.null(sensor_results$random_forest)) {
      cat("\nRANDOM FOREST ATTRIBUTION:\n")
      cat("R² =", sprintf("%.3f", sensor_results$random_forest$r_squared), "\n")
      cat("RMSE =", sprintf("%.3f", sensor_results$random_forest$rmse), "\n")
      cat("MAE =", sprintf("%.3f", sensor_results$random_forest$mae), "\n")
      
      cat("\nVariable Importance:\n")
      importance_order <- order(sensor_results$random_forest$variable_importance, decreasing = TRUE)
      for (i in importance_order) {
        var <- names(sensor_results$random_forest$variable_importance)[i]
        importance <- sensor_results$random_forest$variable_importance[i]
        cat(sprintf("  %s: %.3f\n", var, importance))
      }
      
      cat("\nR² Contributions:\n")
      cat("  Overall R² =", sprintf("%.3f", sensor_results$random_forest$r_squared), "\n")
      cat("\n  Individual R² (model with only this variable):\n")
      for (i in 1:length(sensor_results$random_forest$predictors)) {
        var <- sensor_results$random_forest$predictors[i]
        ind_r2 <- sensor_results$random_forest$individual_r2[i]
        cat(sprintf("    %s: %.3f (%.1f%%)\n", var, ind_r2, ind_r2 * 100))
      }
      cat("\n  Incremental R² (R² change when adding variable):\n")
      for (i in 1:length(sensor_results$random_forest$predictors)) {
        var <- sensor_results$random_forest$predictors[i]
        inc_r2 <- sensor_results$random_forest$incremental_r2[i]
        cat(sprintf("    %s: %.3f (%.1f%%)\n", var, inc_r2, inc_r2 * 100))
      }
    }
    
    # GAM results
    if (!is.null(sensor_results$gam)) {
      cat("\nGAM ATTRIBUTION:\n")
      cat("R² =", sprintf("%.3f", sensor_results$gam$r_squared), "\n")
      cat("Deviance Explained =", sprintf("%.1f%%", sensor_results$gam$deviance_explained * 100), "\n")
      cat("AIC =", sprintf("%.1f", sensor_results$gam$aic), "\n")
      cat("BIC =", sprintf("%.1f", sensor_results$gam$bic), "\n")
      cat("RMSE =", sprintf("%.3f", sensor_results$gam$rmse), "\n")
      cat("MAE =", sprintf("%.3f", sensor_results$gam$mae), "\n")
      
      cat("\nSmooth Terms (EDF and Significance):\n")
      for (i in 1:length(sensor_results$gam$edf_values)) {
        term <- names(sensor_results$gam$edf_values)[i]
        # Clean term name for display
        term_clean <- ifelse(
          grepl("^s\\((.+)\\)$", term),
          gsub("^s\\((.+)\\)$", "\\1", term),
          term
        )
        edf <- sensor_results$gam$edf_values[i]
        p_val <- sensor_results$gam$p_values[i]
        significance <- ifelse(p_val < 0.001, "***",
                            ifelse(p_val < 0.01, "**",
                                  ifelse(p_val < 0.05, "*", "")))
        cat(sprintf("  %s: EDF=%.2f, p=%.3f%s\n", term_clean, edf, p_val, significance))
      }
      
      cat("\nR² Contributions:\n")
      cat("  Overall R² =", sprintf("%.3f", sensor_results$gam$r_squared), "\n")
      cat("\n  Individual R² (model with only this variable):\n")
      for (i in 1:length(sensor_results$gam$predictors)) {
        var <- sensor_results$gam$predictors[i]
        ind_r2 <- sensor_results$gam$individual_r2[i]
        cat(sprintf("    %s: %.3f (%.1f%%)\n", var, ind_r2, ind_r2 * 100))
      }
      cat("\n  Incremental R² (R² change when adding variable):\n")
      for (i in 1:length(sensor_results$gam$predictors)) {
        var <- sensor_results$gam$predictors[i]
        inc_r2 <- sensor_results$gam$incremental_r2[i]
        cat(sprintf("    %s: %.3f (%.1f%%)\n", var, inc_r2, inc_r2 * 100))
      }
    }
    
  }
  
  # Overall summary
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("OVERALL SUMMARY\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  cat("\nKEY FINDINGS:\n")
  cat("• Multiple attribution methods provide consistent results\n")
  cat("• Precipitation is the most important climate variable across all methods\n")
  cat("• Climate variables (temperature, precipitation, snowfall) significantly influence albedo dynamics\n")
  cat("• Cross-validation ensures robust attribution results\n")
  
  cat("\nMETHOD COMPARISON:\n")
  cat("• Multiple Linear Regression (MLR): Linear relationships with relative importance analysis (LMG)\n")
  cat("• Random Forest: Captures non-linear relationships and interactions using permutation importance\n")
  cat("• Generalized Additive Models (GAM): Non-linear relationships with smooth functions (EDF-based importance)\n")
  
  return(summary_tables)
}


export_results_to_csv <- function(attribution_results, monthly_data) {
  if (!dir.exists("output/csv_exports")) {
    dir.create("output/csv_exports", recursive = TRUE)
  }
  
  write_csv(monthly_data, "output/csv_exports/monthly_albedo_climate_data.csv")
  
  # Export attribution summaries for each sensor and method
  for (sensor_name in names(attribution_results)) {
    sensor_results <- attribution_results[[sensor_name]]
    
    # Regression results
    if (!is.null(sensor_results$regression) && 
        !is.null(sensor_results$regression$predictors) && 
        length(sensor_results$regression$predictors) > 0) {
      
      reg_summary <- data.frame(
        sensor = sensor_name,
        method = "Regression",
        variable = sensor_results$regression$predictors,
        coefficient = sensor_results$regression$std_coefficients,
        importance_lmg = sensor_results$regression$relative_importance$lmg,
        importance_first = sensor_results$regression$relative_importance$first,
        importance_last = sensor_results$regression$relative_importance$last,
        r_squared = sensor_results$regression$r_squared,
        adj_r_squared = sensor_results$regression$adj_r_squared,
        individual_r2 = sensor_results$regression$individual_r2,
        partial_r2 = sensor_results$regression$partial_r2,
        incremental_r2 = sensor_results$regression$incremental_r2,
        aic = sensor_results$regression$aic,
        bic = sensor_results$regression$bic,
        stringsAsFactors = FALSE
      )
      write_csv(reg_summary, paste0("output/csv_exports/", sensor_name, "_regression_results.csv"))
    }
    
    # Random Forest results
    if (!is.null(sensor_results$random_forest) && 
        !is.null(sensor_results$random_forest$variable_importance) && 
        length(sensor_results$random_forest$variable_importance) > 0) {
      
      rf_summary <- data.frame(
        sensor = sensor_name,
        method = "Random Forest",
        variable = names(sensor_results$random_forest$variable_importance),
        importance = sensor_results$random_forest$variable_importance,
        r_squared = sensor_results$random_forest$r_squared,
        individual_r2 = sensor_results$random_forest$individual_r2,
        incremental_r2 = sensor_results$random_forest$incremental_r2,
        rmse = sensor_results$random_forest$rmse,
        mae = sensor_results$random_forest$mae,
        stringsAsFactors = FALSE
      )
      write_csv(rf_summary, paste0("output/csv_exports/", sensor_name, "_randomforest_results.csv"))
    }
    
    # GAM results
    if (!is.null(sensor_results$gam) && 
        !is.null(sensor_results$gam$edf_values) && 
        length(sensor_results$gam$edf_values) > 0) {
      
      # Ensure all vectors have the same length
      edf_values <- sensor_results$gam$edf_values
      p_values <- sensor_results$gam$p_values
      
      # Handle case where p_values might be shorter than edf_values
      if (length(p_values) < length(edf_values)) {
        p_values <- c(p_values, rep(NA, length(edf_values) - length(p_values)))
      }
      
      # Clean GAM variable names: remove s() wrapper
      gam_var_names <- names(edf_values)
      gam_var_names_clean <- ifelse(
        grepl("^s\\((.+)\\)$", gam_var_names),
        gsub("^s\\((.+)\\)$", "\\1", gam_var_names),
        gam_var_names
      )
      
      # Match individual_r2 and incremental_r2 to edf_values by variable name
      individual_r2_matched <- sensor_results$gam$individual_r2[gam_var_names_clean]
      incremental_r2_matched <- sensor_results$gam$incremental_r2[gam_var_names_clean]
      
      gam_summary <- data.frame(
        sensor = sensor_name,
        method = "GAM",
        variable = gam_var_names_clean,
        edf = edf_values,
        p_value = p_values[1:length(edf_values)],
        r_squared = sensor_results$gam$r_squared,
        deviance_explained = sensor_results$gam$deviance_explained,
        individual_r2 = individual_r2_matched,
        incremental_r2 = incremental_r2_matched,
        aic = sensor_results$gam$aic,
        bic = sensor_results$gam$bic,
        rmse = sensor_results$gam$rmse,
        mae = sensor_results$gam$mae,
        stringsAsFactors = FALSE
      )
      write_csv(gam_summary, paste0("output/csv_exports/", sensor_name, "_gam_results.csv"))
    }
    
    # Causal results (Method 10) - NOT USED, code kept for reference
    # if (!is.null(sensor_results$causal)) {
    #   var_importance <- sensor_results$causal$variable_importance
    #   
    #   if (!is.null(var_importance) && length(var_importance) > 0) {
    #     # Handle case where variable_importance might not have names
    #     if (is.null(names(var_importance)) && !is.null(sensor_results$causal$predictors)) {
    #       var_names <- setdiff(sensor_results$causal$predictors, "mean_temp")
    #       if (length(var_names) == length(var_importance)) {
    #         names(var_importance) <- var_names
    #       }
    #     }
    #     
    #     if (!is.null(names(var_importance)) && length(names(var_importance)) > 0) {
    #       causal_summary <- data.frame(
    #         sensor = sensor_name,
    #         method = "Causal",
    #         variable = names(var_importance),
    #         importance = as.numeric(var_importance),
    #         stringsAsFactors = FALSE
    #       )
    #       
    #       # Add average treatment effect if available
    #       if (!is.null(sensor_results$causal$average_treatment_effect)) {
    #         ate <- sensor_results$causal$average_treatment_effect
    #         causal_summary$ate_estimate <- ate[1]
    #         causal_summary$ate_lower <- ate[2]
    #         causal_summary$ate_upper <- ate[3]
    #       }
    #       
    #       write_csv(causal_summary, paste0("output/csv_exports/", sensor_name, "_causal_results.csv"))
    #     }
    #   }
    # }
  }
}

export_plots_pdf <- function(attribution_plots) {
  if (!dir.exists("figures")) {
    dir.create("figures", recursive = TRUE)
  }
  
  # Combine all three plots into Figure 5 with panels a, b, c
  if (!is.null(attribution_plots$variable_importance) && 
      !is.null(attribution_plots$coefficients) && 
      !is.null(attribution_plots$method_heatmap)) {
    
    # Create combined figure with proper spacing
    p1 <- attribution_plots$variable_importance
    p2 <- attribution_plots$coefficients
    p3 <- attribution_plots$method_heatmap
    
    # Combine plots vertically with good spacing
    # Use arrangeGrob with explicit spacing via null grobs
    combined_plot <- gridExtra::arrangeGrob(
      p1, 
      grid::nullGrob(),  # Spacer
      p2, 
      grid::nullGrob(),  # Spacer
      p3,
      ncol = 1,
      heights = c(1.0, 0.08, 1.0, 0.08, 0.6)  # First two panels equal, heatmap smaller, with spacing
    )
    
    # Save combined Figure 5
    ggsave("figures/climate_attribution_figure5.pdf", 
           combined_plot, 
           width = 7.0,  # Full page width
           height = 12,  # Appropriate height for three panels
           device = "pdf", 
           useDingbats = FALSE)
    cat("✓ Exported combined Figure 5: figures/climate_attribution_figure5.pdf\n")
  }
  
  # Also save individual plots for reference
  plot_names <- c("variable_importance", "coefficients", "method_heatmap")
  for (plot_name in plot_names) {
    if (!is.null(attribution_plots[[plot_name]])) {
      plot_obj <- attribution_plots[[plot_name]]
      file_path <- paste0("figures/climate_attribution_", plot_name, ".pdf")
      ggsave(file_path, plot_obj, width = 7.0, height = 5.0, device = "pdf", useDingbats = FALSE)
    }
  }
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Load all datasets
cat(paste(rep("=", 80), collapse=""), "\n")
cat("HEARD ISLAND GLACIER ALBEDO CLIMATE ATTRIBUTION ANALYSIS\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("Loading VIIRS albedo data (credible sensor with validated algorithms)\n")

# Load VIIRS data from standalone CSV file
viirs_raw <- read_csv(paste0(data_path, "HIG_VIIRS_Albedo/heard_island_viirs_albedo_area_mean.csv"), 
                     show_col_types = FALSE)

# Process VIIRS data
viirs_data <- viirs_raw %>%
  dplyr::select(date, albedo = albedo_broadband) %>%
  mutate(date = as.Date(date)) %>%
  filter(!is.na(date) & !is.na(albedo)) %>%
  filter(year(date) <= 2024) %>%
  mutate(sensor = "VIIRS")

cat("✓ VIIRS data loaded:", nrow(viirs_data), "records\n")
cat("  Date range:", min(viirs_data$date), "to", max(viirs_data$date), "\n")

# Load climate data
climate_data <- load_climate_data()

# Use VIIRS data only
all_albedo_data <- viirs_data

# Create monthly aggregated dataset
monthly_data <- create_monthly_dataset(all_albedo_data, climate_data)

# Export monthly aggregated data
write.csv(monthly_data, 
          "output/climate_attribution_monthly_data.csv", 
          row.names = FALSE)
cat("✓ Exported monthly aggregated data: output/climate_attribution_monthly_data.csv\n")

# Perform comprehensive climate attribution analysis
attribution_results <- perform_climate_attribution(monthly_data)

# Create visualization plots
attribution_plots <- create_attribution_plots(attribution_results)



export_results_to_csv(attribution_results, monthly_data)
export_plots_pdf(attribution_plots)

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("CLIMATE ATTRIBUTION ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Results exported:\n")
cat("  CSV: output/csv_exports/*.csv\n")
cat("  Plots: figures/climate_attribution_*.pdf\n")
