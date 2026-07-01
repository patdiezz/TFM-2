################################################################################
##################### SENSITIVITY ANALYSIS BY SEX ##############################
##################### MODEL-BASED SENSITIVITY ANALYSIS #########################
################################################################################

library(dplyr)
library(ggplot2)
library(dbrobust)
library(dbstats)
library(haven)
library(stringr)

# ============================================================
# 1. LOAD VARIABLE SETS
# ============================================================

var_sets <- readRDS("var_sets.rds")

to_chr_vec <- function(x) {
  unique(as.character(unlist(x, use.names = FALSE)))
}

cont_vars <- to_chr_vec(var_sets$cont_vars)
bin_vars  <- to_chr_vec(var_sets$bin_vars)
cat_vars  <- to_chr_vec(var_sets$cat_vars)

vars <- to_chr_vec(c(cont_vars, bin_vars, cat_vars))


# ============================================================
# 2. AUXILIARY FUNCTIONS
# ============================================================

to_plain <- function(x) {
  haven::zap_labels(x)
}

to_plain_numeric <- function(x) {
  as.numeric(haven::zap_labels(x))
}

weighted_median <- function(x, w) {
  
  x <- to_plain_numeric(x)
  w <- to_plain_numeric(w)
  
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  
  if (length(x) == 0 || length(w) == 0) {
    return(NA_real_)
  }
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  cw <- cumsum(w) / sum(w)
  
  x[which(cw >= 0.5)[1]]
}

weighted_mode <- function(x, w) {
  
  x <- to_plain(x)
  w <- to_plain_numeric(w)
  
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  
  if (length(x) == 0 || length(w) == 0) {
    return(NA)
  }
  
  tab <- tapply(w, x, sum)
  
  names(tab)[which.max(tab)]
}

safe_filename <- function(x) {
  stringr::str_replace_all(x, "[^A-Za-z0-9_]", "_")
}


# ============================================================
# 3. STRATIFIED TRAIN/TEST SPLIT
# ============================================================

train_test_split_stratified <- function(data, strata_var = "sphus_bin", train_prop = 0.7) {
  
  data <- data %>% mutate(.row_id = row_number())
  
  split_list <- data %>%
    group_by(.data[[strata_var]]) %>%
    group_split()
  
  result <- lapply(split_list, function(df) {
    
    n_train   <- round(nrow(df) * train_prop)
    train_idx <- sample(seq_len(nrow(df)), size = n_train, replace = FALSE)
    
    list(
      train = df[train_idx, ],
      test  = df[-train_idx, ]
    )
  })
  
  list(
    train = bind_rows(lapply(result, `[[`, "train")),
    test  = bind_rows(lapply(result, `[[`, "test"))
  )
}


# ============================================================
# 4. GRID FUNCTION
# ============================================================

make_sensitivity_grid <- function(x, var, cont_vars, grid_n = 30, max_unique_observed = 100) {
  
  x <- to_plain(x)
  
  if (is.numeric(x)) {
    x <- as.numeric(x)
  }
  
  vals_unicos <- sort(unique(x))
  vals_unicos <- vals_unicos[!is.na(vals_unicos)]
  
  if (length(vals_unicos) == 0) {
    return(NULL)
  }
  
  # Binary and categorical variables: always use real observed values
  if (!(var %in% cont_vars)) {
    return(vals_unicos)
  }
  
  # Continuous variables that are actually discrete scales:
  # use real observed values when the number of unique values is not too high
  if (length(vals_unicos) <= max_unique_observed) {
    return(vals_unicos)
  }
  
  # Truly continuous variables with many unique values:
  # use quantile-based grid
  x_num <- as.numeric(x)
  
  grid <- quantile(
    x_num,
    probs = seq(0.05, 0.95, length.out = grid_n),
    na.rm = TRUE,
    names = FALSE
  )
  
  grid <- sort(unique(as.numeric(grid)))
  
  return(grid)
}


# ============================================================
# 5. CREATE BASE INDIVIDUAL BY SEX
# ============================================================

create_base_individual <- function(data_train, vars, cont_vars, bin_vars, cat_vars, gender_value) {
  
  x_base <- data_train[1, vars, drop = FALSE]
  
  # Continuous variables: weighted median
  for (v in cont_vars) {
    x_base[[v]] <- weighted_median(
      data_train[[v]],
      data_train$cciw_w9
    )
  }
  
  # Binary variables: weighted mode
  for (v in bin_vars) {
    x_base[[v]] <- weighted_mode(
      data_train[[v]],
      data_train$cciw_w9
    )
  }
  
  # Categorical variables: weighted mode
  for (v in cat_vars) {
    x_base[[v]] <- weighted_mode(
      data_train[[v]],
      data_train$cciw_w9
    )
  }
  
  # Fix sex profile
  # gender: 1 = Male, 2 = Female
  x_base$gender <- gender_value
  
  # Force same types as training set
  for (v in vars) {
    
    if (is.numeric(data_train[[v]]) || inherits(data_train[[v]], "haven_labelled")) {
      
      x_base[[v]] <- as.numeric(x_base[[v]])
      
    } else if (is.factor(data_train[[v]])) {
      
      x_base[[v]] <- factor(
        x_base[[v]],
        levels = levels(data_train[[v]])
      )
      
    } else {
      
      x_base[[v]] <- as.character(x_base[[v]])
    }
  }
  
  return(x_base)
}


# ============================================================
# 6. PREDICTION FUNCTION:
#    MALE + FEMALE TOGETHER, ONE GENERALIZED GOWER MATRIX PER VARIABLE
# ============================================================

# The baseline profile is defined separately for men and women.
# Continuous predictors are fixed at their weighted median,
# while binary and categorical predictors are fixed at their weighted mode.

predecir_grid_variable_ggower_two_profiles <- function(var, grid, modelo, data_train,
                                                       cont_vars, bin_vars, cat_vars,
                                                       x_base_male, x_base_female,
                                                       peso_nuevo = 1,
                                                       threshold_main = 0.60) {
  
  vars <- to_chr_vec(c(cont_vars, bin_vars, cat_vars))
  B <- length(modelo$models_list)
  
  # ----------------------------------------------------------
  # Create male grid
  # ----------------------------------------------------------
  
  new_grid_male <- x_base_male[rep(1, length(grid)), ]
  new_grid_male[[var]] <- grid
  new_grid_male$sex_profile <- "Male"
  
  # ----------------------------------------------------------
  # Create female grid
  # ----------------------------------------------------------
  
  new_grid_female <- x_base_female[rep(1, length(grid)), ]
  new_grid_female[[var]] <- grid
  new_grid_female$sex_profile <- "Female"
  
  # ----------------------------------------------------------
  # Join both profiles
  # ----------------------------------------------------------
  
  new_grid <- bind_rows(new_grid_male, new_grid_female)
  
  sex_profile <- new_grid$sex_profile
  new_grid$sex_profile <- NULL
  
  # Force same types as training set
  for (v in vars) {
    
    if (is.numeric(data_train[[v]]) || inherits(data_train[[v]], "haven_labelled")) {
      
      new_grid[[v]] <- as.numeric(new_grid[[v]])
      
    } else if (is.factor(data_train[[v]])) {
      
      new_grid[[v]] <- factor(
        new_grid[[v]],
        levels = levels(data_train[[v]])
      )
      
    } else {
      
      new_grid[[v]] <- as.character(new_grid[[v]])
    }
  }
  
  n_new <- nrow(new_grid)
  
  x_train_full <- data_train[, vars, drop = FALSE]
  w_train_full <- as.numeric(haven::zap_labels(data_train$cciw_w9))
  
  # ----------------------------------------------------------
  # Synthetic individuals + training set
  # Synthetic individuals have weight 1
  # ----------------------------------------------------------
  
  x_combined <- dplyr::bind_rows(new_grid, x_train_full)
  w_combined <- c(rep(peso_nuevo, n_new), w_train_full)
  
  # ----------------------------------------------------------
  # GENERALIZED GOWER DISTANCE
  # ----------------------------------------------------------
  
  # A new G-Gower distance matrix is computed for the synthetic grid
  # together with the training observations, so that predictions are
  # obtained in the same distance-based space used by the fitted models.
  
  D_all_ggower <- dbrobust::robust_distances(
    data      = x_combined,
    cont_vars = cont_vars,
    bin_vars  = bin_vars,
    cat_vars  = cat_vars,
    w         = w_combined,
    alpha     = 0.10,
    method    = "ggower"
  )
  
  D_all_ggower_euc <- dbrobust::make_euclidean(
    D_all_ggower,
    w = w_combined
  )$D_euc
  
  D_all_ggower_D2 <- D_all_ggower_euc^2
  
  pred_mat <- matrix(NA, nrow = n_new, ncol = B)
  
  #Predictions are obtained from all models in the sub-bagging ensemble
  #and then averaged to obtain the final predicted probability.
  
  for (b in 1:B) {
    
    idx_sub <- modelo$idx_sub_list[[b]]
    
    cols_b <- n_new + idx_sub
    
    D_cross_b <- D_all_ggower_D2[1:n_new, cols_b, drop = FALSE]
    class(D_cross_b) <- "D2"
    
    pred_mat[, b] <- predict(
      modelo$models_list[[b]],
      newdata   = D_cross_b,
      type.pred = "response",
      type.var  = "D2"
    )
  }
  
  prob <- rowMeans(pred_mat, na.rm = TRUE)
  
  data.frame(
    variable    = var,
    valor       = new_grid[[var]],
    prob        = prob,
    clase_06    = as.integer(prob >= threshold_main),
    clase_05    = as.integer(prob >= 0.50),
    sex_profile = sex_profile
  )
}


# ============================================================
# 7. SUMMARY FUNCTION
# ============================================================

summarise_sensitivity <- function(sens_var, threshold_main = 0.60) {
  
  sens_var %>%
    group_by(variable, sex_profile) %>%
    arrange(valor, .by_group = TRUE) %>%
    mutate(
      clase_06  = as.integer(prob >= threshold_main),
      clase_05  = as.integer(prob >= 0.50),
      cambio_06 = clase_06 != lag(clase_06),
      cambio_05 = clase_05 != lag(clase_05)
    ) %>%
    summarise(
      min_prob   = min(prob, na.rm = TRUE),
      max_prob   = max(prob, na.rm = TRUE),
      delta_prob = max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE),
      
      cruza_06 = any(cambio_06 == TRUE, na.rm = TRUE),
      umbral_06_aprox = ifelse(
        any(cambio_06 == TRUE, na.rm = TRUE),
        as.character(valor[which(cambio_06 == TRUE)[1]]),
        NA
      ),
      
      cruza_05 = any(cambio_05 == TRUE, na.rm = TRUE),
      umbral_05_aprox = ifelse(
        any(cambio_05 == TRUE, na.rm = TRUE),
        as.character(valor[which(cambio_05 == TRUE)[1]]),
        NA
      ),
      
      .groups = "drop"
    )
}


# ============================================================
# 8. PLOT FUNCTION: MALE + FEMALE IN SAME GRAPH
# ============================================================

plot_sensitivity_pretty <- function(sens_var, var_name, country_name, output_path,
                                    threshold_main = 0.60) {
  
  p <- ggplot(
    sens_var,
    aes(
      x = valor,
      y = prob,
      color = sex_profile,
      group = sex_profile
    )
  ) +
    geom_hline(
      yintercept = 0.50,
      linetype = "dotted",
      linewidth = 0.45,
      color = "grey45"
    ) +
    geom_hline(
      yintercept = threshold_main,
      linetype = "dashed",
      linewidth = 0.50,
      color = "#D62728"
    ) +
    geom_line(linewidth = 0.65, alpha = 0.95) +
    geom_point(size = 1.4, alpha = 0.95) +
    scale_color_manual(
      values = c(
        "Male"   = "#1F78B4",
        "Female" = "#E31A1C"
      )
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
    labs(
      title    = paste0("Sensitivity analysis: ", var_name),
      subtitle = paste0(country_name, " | Sensitivity analysis"),
      x        = var_name,
      y        = "Predicted probability of poor self-perceived health",
      color    = "Baseline profile",
      caption  = "Dashed red line: selected threshold = 0.60 | Dotted grey line: reference threshold = 0.50"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey25"),
      plot.caption = element_text(size = 9, color = "grey35", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  ggsave(
    filename = output_path,
    plot     = p,
    width    = 7.5,
    height   = 4.8,
    dpi      = 300
  )
  
  return(p)
}


# ============================================================
# 9. RUN SENSITIVITY ANALYSIS FOR ONE COUNTRY
# ============================================================

run_sensitivity_country <- function(country_name,
                                    data_file,
                                    model_file,
                                    output_dir,
                                    threshold_main = 0.60,
                                    grid_n = 30,
                                    vars_to_run = NULL) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("Running sensitivity analysis for:", country_name, "\n")
  cat(strrep("=", 70), "\n")
  
  # ----------------------------------------------------------
  # Load data and model
  # ----------------------------------------------------------
  
  data_model <- readRDS(data_file)
  modelo <- readRDS(model_file)
  
  # ----------------------------------------------------------
  # Define variables robustly
  # ----------------------------------------------------------
  
  if (!is.null(modelo$vars_model)) {
    vars <- to_chr_vec(modelo$vars_model)
  } else {
    vars <- to_chr_vec(c(cont_vars, bin_vars, cat_vars))
  }
  
  vars <- vars[vars %in% names(data_model)]
  
  cont_vars_run <- intersect(cont_vars, vars)
  bin_vars_run  <- intersect(bin_vars, vars)
  cat_vars_run  <- intersect(cat_vars, vars)
  
  cat("Number of variables used:", length(vars), "\n")
  
  # ----------------------------------------------------------
  # Remove SPSS labels safely
  # ----------------------------------------------------------
  
  cols_to_zap <- unique(c(vars, "sphus_bin", "cciw_w9"))
  cols_to_zap <- cols_to_zap[cols_to_zap %in% names(data_model)]
  
  data_model <- data_model %>%
    mutate(across(
      all_of(cols_to_zap),
      ~ haven::zap_labels(.x)
    ))
  
  # ----------------------------------------------------------
  # Use original train/test split if available
  # Otherwise, recreate same split
  # ----------------------------------------------------------
  
  country_key <- tolower(country_name)
  
  train_object_name <- paste0("data_", country_key, "_train")
  test_object_name  <- paste0("data_", country_key, "_test")
  
  if (exists(train_object_name, envir = .GlobalEnv) &&
      exists(test_object_name,  envir = .GlobalEnv)) {
    
    cat("Using existing train/test objects from environment.\n")
    
    data_train <- get(train_object_name, envir = .GlobalEnv)
    data_test  <- get(test_object_name,  envir = .GlobalEnv)
    
    data_train <- data_train %>%
      mutate(across(
        all_of(cols_to_zap[cols_to_zap %in% names(data_train)]),
        ~ haven::zap_labels(.x)
      ))
    
    data_test <- data_test %>%
      mutate(across(
        all_of(cols_to_zap[cols_to_zap %in% names(data_test)]),
        ~ haven::zap_labels(.x)
      ))
    
  } else {
    
    cat("Train/test objects not found. Recreating split with set.seed(1234).\n")
    
    set.seed(1234)
    split_country <- train_test_split_stratified(data_model)
    
    data_train <- split_country$train
    data_test  <- split_country$test
  }
  
  cat("Train size:", nrow(data_train), "\n")
  cat("Test size:", nrow(data_test), "\n")
  
  # ----------------------------------------------------------
  # Output folder
  # ----------------------------------------------------------
  
  dir.create(output_dir, showWarnings = FALSE)
  
  dir_all <- file.path(output_dir, "all_plots")
  dir.create(dir_all, showWarnings = FALSE)
  
  # ----------------------------------------------------------
  # Create male and female baseline individuals
  # ----------------------------------------------------------
  
  x_base_male <- create_base_individual(
    data_train   = data_train,
    vars         = vars,
    cont_vars    = cont_vars_run,
    bin_vars     = bin_vars_run,
    cat_vars     = cat_vars_run,
    gender_value = 1
  )
  
  x_base_female <- create_base_individual(
    data_train   = data_train,
    vars         = vars,
    cont_vars    = cont_vars_run,
    bin_vars     = bin_vars_run,
    cat_vars     = cat_vars_run,
    gender_value = 2
  )
  
  # Synthetic individuals get the median sampling weight from the training set
  peso_sintetico <- median(
    as.numeric(haven::zap_labels(data_train$cciw_w9)),
    na.rm = TRUE
  )
  
  cat("Synthetic weight:", peso_sintetico, "\n")
  
  # Variables to analyse
  if (is.null(vars_to_run)) {
    vars_sens <- setdiff(vars, "gender")
  } else {
    vars_sens <- intersect(vars_to_run, setdiff(vars, "gender"))
  }
  
  cat("Variables to run:", length(vars_sens), "\n")
  
  resultados <- list()
  tabla_resumen <- data.frame()
  
  # ----------------------------------------------------------
  # Loop over variables
  # ----------------------------------------------------------
  
  for (var in vars_sens) {
    
    cat("\n============================\n")
    cat("Variable:", var, "\n")
    cat("============================\n")
    
    grid <- make_sensitivity_grid(
      x                   = data_train[[var]],
      var                 = var,
      cont_vars           = cont_vars_run,
      grid_n              = grid_n,
      max_unique_observed = 100
    )
    
    if (is.null(grid) || length(grid) == 0) {
      cat("Skipping", var, "- no valid values\n")
      next
    }
    
    cat("Grid points:", length(grid), "\n")
    print(grid)
    
    # --------------------------------------------------------
    # Prediction
    # --------------------------------------------------------
    
    sens_var <- predecir_grid_variable_ggower_two_profiles(
      var            = var,
      grid           = grid,
      modelo         = modelo,
      data_train     = data_train,
      cont_vars      = cont_vars_run,
      bin_vars       = bin_vars_run,
      cat_vars       = cat_vars_run,
      x_base_male    = x_base_male,
      x_base_female  = x_base_female,
      peso_nuevo     = peso_sintetico,
      threshold_main = threshold_main
    )
    
    sens_var$country <- country_name
    
    resultados[[var]] <- sens_var
    
    # --------------------------------------------------------
    # Summary table
    # --------------------------------------------------------
    
    resumen_var <- summarise_sensitivity(
      sens_var       = sens_var,
      threshold_main = threshold_main
    ) %>%
      mutate(country = country_name) %>%
      select(country, everything())
    
    tabla_resumen <- bind_rows(tabla_resumen, resumen_var)
    
    # --------------------------------------------------------
    # Save plot
    # --------------------------------------------------------
    
    safe_var <- safe_filename(var)
    
    output_path <- file.path(
      dir_all,
      paste0(country_name, "_", safe_var, ".png")
    )
    
    plot_sensitivity_pretty(
      sens_var       = sens_var,
      var_name       = var,
      country_name   = country_name,
      output_path    = output_path,
      threshold_main = threshold_main
    )
    
    cat("Done:", var, "\n")
  }
  
  # ============================================================
  # 10. SUMMARY OF VARIABLES CROSSING THRESHOLDS
  # ============================================================
  
  variables_cruzan <- tabla_resumen %>%
    group_by(country, variable) %>%
    summarise(
      cruza_06_any = any(cruza_06, na.rm = TRUE),
      cruza_05_any = any(cruza_05, na.rm = TRUE),
      max_delta    = max(delta_prob, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(cruza_06_any), desc(max_delta))
  
  variables_main <- variables_cruzan %>%
    filter(cruza_06_any == TRUE)
  
  # ============================================================
  # 11. SAVE OUTPUTS
  # ============================================================
  
  tabla_resumen <- tabla_resumen %>%
    arrange(variable, sex_profile)
  
  write.csv(
    tabla_resumen,
    file.path(output_dir, paste0("tabla_resumen_sensibilidad_", country_name, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    variables_cruzan,
    file.path(output_dir, paste0("variables_cruzan_sensibilidad_", country_name, ".csv")),
    row.names = FALSE
  )
  
  saveRDS(
    resultados,
    file.path(output_dir, paste0("resultados_sensibilidad_", country_name, ".rds"))
  )
  
  cat("\n", strrep("=", 70), "\n")
  cat("Finished sensitivity analysis for:", country_name, "\n")
  cat(strrep("=", 70), "\n")
  
  cat("\nVariables crossing threshold 0.60:\n")
  print(variables_main)
  
  invisible(list(
    resultados       = resultados,
    tabla_resumen    = tabla_resumen,
    variables_cruzan = variables_cruzan,
    variables_main   = variables_main
  ))
}


################################################################################
############################ RUN EXAMPLES ######################################
################################################################################

# ============================================================
# LATVIA
# ============================================================

sens_latvia <- run_sensitivity_country(
  country_name    = "Latvia",
  data_file       = "data_latvia_model.rds",
  model_file      = "subbagging_latvia_ggower_70.rds",
  output_dir      = "sensitivity_latvia_ggower3",
  threshold_main  = 0.60,
  grid_n          = 30
)

# ============================================================
# SPAIN
# ============================================================

sens_spain <- run_sensitivity_country(
  country_name    = "Spain",
  data_file       = "data_spain_model.rds",
  model_file      = "subbagging_spain_ggower_70.rds",
  output_dir      = "sensitivity_spain_ggower4",
  threshold_main  = 0.60,
  grid_n          = 30
)

# ============================================================
# FINLAND
# ============================================================

sens_finland <- run_sensitivity_country(
  country_name    = "Finland",
  data_file       = "data_finland_model.rds",
  model_file      = "subbagging_finland_ggower_70.rds",
  output_dir      = "sensitivity_finland_ggower4",
  threshold_main  = 0.60,
  grid_n          = 30
)


# ============================================================
# GERMANY
# ============================================================

sens_germany <- run_sensitivity_country(
  country_name    = "Germany",
  data_file       = "data_germany_model.rds",
  model_file      = "subbagging_germany_ggower_70.rds",
  output_dir      = "sensitivity_germany_ggower4",
  threshold_main  = 0.60,
  grid_n          = 30
)
