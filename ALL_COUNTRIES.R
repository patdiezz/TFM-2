################################################################################
          ############### SUB-BAGGING MODELS ###########
################################################################################

#Upload the librarys 
library(dplyr)
library(tidyr)
library(dbrobust)
library(dbstats)
library(pROC)

#Set.seed for reproducibility
set.seed(1234)

# ============================================================
# 0. LOAD VARIABLE SETS
# ============================================================

var_sets <- readRDS("var_sets.rds") 

cont_vars <- var_sets$cont_vars
bin_vars  <- var_sets$bin_vars
cat_vars  <- var_sets$cat_vars

vars <- c(cont_vars, bin_vars, cat_vars)


# ============================================================
# 1. GENERAL PARAMETERS
# ============================================================

B        <- 30  #number of bags
bag_size <- 200 #size
gvar     <- 0.70 #geometric variability


# ============================================================
# 2. STRATIFIED TRAIN/TEST SPLIT FUNCTION
# ============================================================

#So train and test have the same proportion of classes
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
# 3. EVALUATION FUNCTION
# ============================================================

#Evaluate 4 thresholds (0.5, 0.6, 0.7, 0.75)
evaluate_threshold_unweighted <- function(probs, y_true, threshold) {
  
  pred_class <- as.integer(probs >= threshold)
  
  conf_mat <- table(
    Predicted = factor(pred_class, levels = c(0, 1)),
    Real      = factor(y_true,     levels = c(0, 1))
  )
  
  TP <- sum(pred_class == 1 & y_true == 1, na.rm = TRUE)
  TN <- sum(pred_class == 0 & y_true == 0, na.rm = TRUE)
  FP <- sum(pred_class == 1 & y_true == 0, na.rm = TRUE)
  FN <- sum(pred_class == 0 & y_true == 1, na.rm = TRUE)
  
  accuracy <- (TP + TN) / length(y_true)
  
  sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  
  ppv <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  npv <- ifelse((TN + FN) == 0, NA, TN / (TN + FN))
  
  baseline_1 <- mean(y_true == 1, na.rm = TRUE)
  baseline_majority <- max(prop.table(table(y_true)))
  
  invisible(list(
    threshold         = threshold,
    conf_mat          = conf_mat,
    accuracy          = accuracy,
    sensitivity       = sensitivity,
    specificity       = specificity,
    ppv               = ppv,
    npv               = npv,
    baseline_1        = baseline_1,
    baseline_majority = baseline_majority
  ))
}


# ============================================================
# 4. FUNCTION TO RUN ONE COUNTRY AND ONE DISTANCE
# ============================================================

run_country_subbagging <- function(country,
                                   distance_type = c("ggower", "classic_gower"),
                                   B = 30,
                                   bag_size = 200,
                                   gvar = 0.70) {
  
  distance_type <- match.arg(distance_type)
  
  cat("\n", strrep("#", 90), "\n")
  cat("COUNTRY:", toupper(country), "| MODEL:", distance_type, "\n")
  cat(strrep("#", 90), "\n")
  
  # ------------------------------------------------------------
  # 4.1 Load country data
  # ------------------------------------------------------------
  
  data_file <- paste0("data_", country, "_model.rds")
  data_country_model <- readRDS(data_file)
  
  set.seed(1234)
  split_country <- train_test_split_stratified(data_country_model)
  
  data_train <- split_country$train
  data_test  <- split_country$test
  
  cat("\nTrain proportions:\n")
  print(prop.table(table(data_train$sphus_bin)))
  
  cat("\nTest proportions:\n")
  print(prop.table(table(data_test$sphus_bin)))
  
  # ------------------------------------------------------------
  # 4.2 Prepare data
  # ------------------------------------------------------------
  
  x_test  <- data_test  %>% select(all_of(vars))
  x_train <- data_train %>% select(all_of(vars))
  
  y_test  <- as.integer(as.character(data_test$sphus_bin))
  y_train <- as.integer(as.character(data_train$sphus_bin))
  
  w_test  <- data_test$cciw_w9
  w_train <- data_train$cciw_w9
  
  n_test  <- nrow(data_test)
  n_train <- nrow(data_train)
  
  X_all <- bind_rows(x_test, x_train)
  w_all <- c(w_test, w_train)
  
  idx_test_all  <- 1:n_test
  idx_train_all <- (n_test + 1):(n_test + n_train)
  
  prop_strata <- prop.table(table(data_train$sphus_bin))
  
  # ------------------------------------------------------------
  # 4.3 Start total execution time
  # ------------------------------------------------------------
  
  total_time_start <- proc.time()
  
  # ------------------------------------------------------------
  # 4.4 Compute distance matrix
  # ------------------------------------------------------------
  
  if (distance_type == "ggower") {
    
    cat("\nComputing generalized Gower distance matrix...\n")
    
    D_all_raw <- robust_distances(
      data      = X_all,
      cont_vars = cont_vars,
      bin_vars  = bin_vars,
      cat_vars  = cat_vars,
      w         = w_all,
      alpha     = 0.10,
      method    = "ggower"
    )
    
    cat("Applying Euclidean correction to generalized Gower...\n")
    
    D_all_euc <- make_euclidean(
      D_all_raw,
      w = w_all
    )$D_euc
    
    D_all_D2 <- D_all_euc^2
    
    model_print_name <- "Generalized Gower"
    distance_label <- "generalized_gower_make_euclidean"
    euclidean_correction <- TRUE
    
  } else if (distance_type == "classic_gower") {
    
    cat("\nComputing classic Gower distance matrix...\n")
    
    D_all_raw <- calculate_distances(
      x                = X_all,
      method           = "gower",
      output_format    = "matrix",
      squared          = FALSE,
      continuous_cols  = cont_vars,
      binary_cols      = bin_vars,
      categorical_cols = cat_vars
    )
    
    cat("Classic Gower: no Euclidean correction applied.\n")
    
    D_all_D2 <- D_all_raw^2
    
    model_print_name <- "Classic Gower"
    distance_label <- "classic_gower_no_euclidean_correction"
    euclidean_correction <- FALSE
  }
  
  cat("\nDistance matrix dimensions:", dim(D_all_D2), "\n")
  
  # ------------------------------------------------------------
  # 4.5 Prepare sub-bagging objects
  # ------------------------------------------------------------
  
  pred_matrix <- matrix(NA, nrow = n_test, ncol = B)
  
  models_list  <- vector("list", B)
  idx_sub_list <- vector("list", B)
  
  dims_used       <- numeric(B)
  relgvar_reached <- numeric(B)
  
  # ------------------------------------------------------------
  # 4.6 Sub-bagging loop
  # ------------------------------------------------------------
  
  set.seed(1234)
  
  for (b in 1:B) {
    
    cat("\n", country, "|", model_print_name, "| Bag", b, "of", B, "\n")
    
    idx_sub <- data_train %>%
      mutate(.row = row_number()) %>%
      group_by(sphus_bin) %>%
      group_modify(~ {
        
        stratum_value <- .y$sphus_bin
        
        n_b <- round(
          as.numeric(prop_strata[as.character(stratum_value)]) * bag_size
        )
        
        slice_sample(
          .x,
          n = min(n_b, nrow(.x)),
          replace = FALSE
        )
        
      }) %>%
      ungroup() %>%
      pull(.row)
    
    y_b <- as.integer(as.character(data_train$sphus_bin[idx_sub]))
    w_b <- data_train$cciw_w9[idx_sub]
    
    idx_train_b_all <- idx_train_all[idx_sub]
    
    D_train_b_D2 <- D_all_D2[idx_train_b_all, idx_train_b_all]
    class(D_train_b_D2) <- "D2"
    
    model_b <- dbglm(
      D2       = D_train_b_D2,
      y        = y_b,
      family   = binomial(link = "logit"),
      method   = "rel.gvar",
      rel.gvar = gvar,
      weights  = w_b
    )
    
    dims_used[b]       <- model_b$eff.rank
    relgvar_reached[b] <- model_b$rel.gvar
    
    D_cross_D2 <- D_all_D2[idx_test_all, idx_train_b_all]
    class(D_cross_D2) <- "D2"
    
    pred_matrix[, b] <- predict(
      model_b,
      newdata   = D_cross_D2,
      type.pred = "response",
      type.var  = "D2"
    )
    
    models_list[[b]]  <- model_b
    idx_sub_list[[b]] <- idx_sub
    
    cat("  Bag size:", length(idx_sub), "\n")
    cat("  Dimensions used:", dims_used[b], "\n")
    cat("  Actual rel.gvar:", round(relgvar_reached[b], 4), "\n")
  }
  
  # ------------------------------------------------------------
  # 4.7 Total execution time
  # ------------------------------------------------------------
  
  total_time <- (proc.time() - total_time_start)["elapsed"]
  
  cat("\nTotal execution time:", round(total_time, 2), "seconds\n")
  
  # ------------------------------------------------------------
  # 4.8 Aggregate predictions
  # ------------------------------------------------------------
  
  prob_bagging <- rowMeans(pred_matrix, na.rm = TRUE)
  
  cat("\nNumber of NA predictions:\n")
  print(sum(is.na(pred_matrix)))
  
  cat("\nNA predictions per bag:\n")
  print(colSums(is.na(pred_matrix)))
  
  # ------------------------------------------------------------
  # 4.9 Unweighted evaluation at several thresholds
  # ------------------------------------------------------------
  
  results_50 <- evaluate_threshold_unweighted(
    probs     = prob_bagging,
    y_true    = y_test,
    threshold = 0.50
  )
  
  results_60 <- evaluate_threshold_unweighted(
    probs     = prob_bagging,
    y_true    = y_test,
    threshold = 0.60
  )
  
  results_70 <- evaluate_threshold_unweighted(
    probs     = prob_bagging,
    y_true    = y_test,
    threshold = 0.70
  )
  
  results_75 <- evaluate_threshold_unweighted(
    probs     = prob_bagging,
    y_true    = y_test,
    threshold = 0.75
  )
  
  summary_thresholds <- data.frame(
    Country     = country,
    Model       = model_print_name,
    Threshold   = c(0.50, 0.60, 0.70, 0.75),
    Accuracy    = c(results_50$accuracy,    results_60$accuracy,
                    results_70$accuracy,    results_75$accuracy),
    Sensitivity = c(results_50$sensitivity, results_60$sensitivity,
                    results_70$sensitivity, results_75$sensitivity),
    Specificity = c(results_50$specificity, results_60$specificity,
                    results_70$specificity, results_75$specificity),
    PPV         = c(results_50$ppv,         results_60$ppv,
                    results_70$ppv,         results_75$ppv),
    NPV         = c(results_50$npv,         results_60$npv,
                    results_70$npv,         results_75$npv),
    Baseline_1  = c(results_50$baseline_1, results_60$baseline_1,
                    results_70$baseline_1, results_75$baseline_1),
    Baseline_Majority = c(results_50$baseline_majority,
                          results_60$baseline_majority,
                          results_70$baseline_majority,
                          results_75$baseline_majority)
  )
  
  cat("\nSummary thresholds - unweighted evaluation:\n")
  print(
    summary_thresholds %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  )
  
  # ------------------------------------------------------------
  # 4.10 ROC and AUC
  # ------------------------------------------------------------
  
  roc_unweighted <- pROC::roc(
    response  = y_test,
    predictor = prob_bagging,
    levels    = c(0, 1),
    direction = "<",
    quiet     = TRUE
  )
  
  auc_unweighted <- as.numeric(pROC::auc(roc_unweighted))
  
  cat("\nUnweighted AUC:", round(auc_unweighted, 3), "\n")
  
  # ------------------------------------------------------------
  # 4.11 Dimensions by bag  (in case we need them)
  # ------------------------------------------------------------
  
  dimensions_by_bag <- data.frame(
    Country = country,
    Model = model_print_name,
    Bag = 1:B,
    Dimensions = dims_used,
    Actual_rel_gvar = relgvar_reached
  )
  
  # ------------------------------------------------------------
  # 4.12 Final computational summary 
  # ------------------------------------------------------------
  
  computational_summary <- data.frame(
    Country = country,
    Model = model_print_name,
    Distance = distance_label,
    Euclidean_correction = euclidean_correction,
    Requested_rel_gvar = gvar,
    Mean_actual_rel_gvar = mean(relgvar_reached),
    Mean_dimensions = mean(dims_used),
    Min_dimensions = min(dims_used),
    Max_dimensions = max(dims_used),
    Total_time_seconds = as.numeric(total_time), #computational cost
    Total_time_minutes = as.numeric(total_time) / 60,
    AUC_unweighted = auc_unweighted
  )
  
  cat("\nComputational summary:\n")
  print(
    computational_summary %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  )
  
  # ------------------------------------------------------------
  # 4.13 Final object
  # ------------------------------------------------------------
  
  result <- list(
    country = country,
    model = model_print_name,
    
    models_list = models_list,
    idx_sub_list = idx_sub_list,
    
    pred_matrix = pred_matrix,
    prob_bagging = prob_bagging,
    
    results_50 = results_50,
    results_60 = results_60,
    results_70 = results_70,
    results_75 = results_75,
    
    summary_thresholds = summary_thresholds,
    
    roc_unweighted = roc_unweighted,
    auc_unweighted = auc_unweighted,
    
    dimensions_by_bag = dimensions_by_bag,
    computational_summary = computational_summary,
    
    gvar = gvar,
    B = B,
    bag_size = bag_size,
    vars_model = vars,
    distance = distance_label,
    euclidean_correction = euclidean_correction,
    evaluation = "unweighted"
  )
  
  # ------------------------------------------------------------
  # 4.14 Save individual object
  # ------------------------------------------------------------
  
  output_file <- paste0(
    "subbagging_",
    country,
    "_",
    ifelse(distance_type == "ggower", "ggower", "classic_gower"),
    "_70.rds"
  )
  
  saveRDS(result, output_file)
  
  cat("\nSaved object:", output_file, "\n")
  
  return(result)
}


################################################################################
############################ 5. RUN MODELS ONE BY ONE ###########################
################################################################################

#Now we apply the function to each country with G-Gower and Classic Gower.
# ============================================================
# LATVIA
# ============================================================

# Latvia - Generalized Gower
result_latvia_ggower <- run_country_subbagging(
  country       = "latvia",
  distance_type = "ggower",
  B             = B,
  bag_size      = bag_size,
  gvar          = gvar
)

# Latvia - Classic Gower
result_latvia_classic_gower <- run_country_subbagging(
  country       = "latvia",
  distance_type = "classic_gower",
  B             = B,
  bag_size      = bag_size,
  gvar          = gvar
)


# ============================================================
# SPAIN
# ============================================================

# Spain - Generalized Gower
result_spain_ggower <- run_country_subbagging(
  country       = "spain",
  distance_type = "ggower",
  B             = B,
  bag_size      = bag_size,
  gvar          = gvar
)

# Spain - Classic Gower
result_spain_classic_gower <- run_country_subbagging(
  country       = "spain",
  distance_type = "classic_gower",
  B             = B,
  bag_size      = bag_size,
  gvar          = gvar
)


# ============================================================
# GERMANY
# ============================================================

# Germany - Generalized Gower
result_germany_ggower <- run_country_subbagging(
  country       = "germany",
  distance_type = "ggower",
  B             = B,
  bag_size      = bag_size,
  gvar          = gvar
)

# Germany - Classic Gower
result_germany_classic_gower <- run_country_subbagging(
  country       = "germany",
  distance_type = "classic_gower",
  B             = B,
  bag_size      = bag_size,
  gvar          = gvar
)


# ============================================================
# FINLAND
# ============================================================

# Finland - Generalized Gower
result_finland_ggower <- run_country_subbagging(
  country       = "finland",
  distance_type = "ggower",
  B             = B,
  bag_size      = bag_size,
  gvar          = gvar
)

# Finland - Classic Gower
result_finland_classic_gower <- run_country_subbagging(
  country       = "finland",
  distance_type = "classic_gower",
  B             = B,
  bag_size      = bag_size,
  gvar          = gvar
)





################################################################################
############ BALANCED SUB-BAGGING MODELS - GENERALIZED GOWER ONLY ##############
########################## UNWEIGHTED EVALUATION ###############################
################################################################################

#A NOTE ON THE BALANCED CASE- What would happen if the classes weren´t imbalanced?
#We apply the same procedure , but we take random samples of the imbalanced class ,
#with the balanced class size. 

set.seed(1234)

# ============================================================
# 0. LOAD VARIABLE SETS
# ============================================================

var_sets <- readRDS("var_sets.rds")

cont_vars <- var_sets$cont_vars
bin_vars  <- var_sets$bin_vars
cat_vars  <- var_sets$cat_vars

vars <- c(cont_vars, bin_vars, cat_vars)

# ============================================================
# 1. GENERAL PARAMETERS
# ============================================================

B    <- 30
gvar <- 0.70

thresholds_eval <- c(0.50, 0.60, 0.70, 0.75) #Only going to look at 0.5

# ============================================================
# 2. STRATIFIED TRAIN/TEST SPLIT FUNCTION
# ============================================================

train_test_split_stratified <- function(data, strata_var = "sphus_bin", train_prop = 0.7) {
  
  data <- data %>% 
    mutate(.row_id = row_number())
  
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
# 3. UNWEIGHTED EVALUATION FUNCTION
# ============================================================

evaluate_threshold_unweighted <- function(probs, y_true, threshold) {
  
  pred_class <- as.integer(probs >= threshold)
  
  conf_mat <- table(
    Predicted = factor(pred_class, levels = c(0, 1)),
    Real      = factor(y_true,     levels = c(0, 1))
  )
  
  TP <- sum(pred_class == 1 & y_true == 1, na.rm = TRUE)
  TN <- sum(pred_class == 0 & y_true == 0, na.rm = TRUE)
  FP <- sum(pred_class == 1 & y_true == 0, na.rm = TRUE)
  FN <- sum(pred_class == 0 & y_true == 1, na.rm = TRUE)
  
  accuracy <- (TP + TN) / length(y_true)
  
  sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  
  ppv <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  npv <- ifelse((TN + FN) == 0, NA, TN / (TN + FN))
  
  baseline_1 <- mean(y_true == 1, na.rm = TRUE)
  baseline_majority <- max(prop.table(table(y_true)))
  
  invisible(list(
    threshold         = threshold,
    conf_mat          = conf_mat,
    accuracy          = accuracy,
    sensitivity       = sensitivity,
    specificity       = specificity,
    ppv               = ppv,
    npv               = npv,
    baseline_1        = baseline_1,
    baseline_majority = baseline_majority
  ))
}

# ============================================================
# 4. FUNCTION TO RUN BALANCED GENERALIZED GOWER MODEL
# ============================================================

run_country_balanced_ggower <- function(country,
                                        B = 2,
                                        gvar = 0.70,
                                        thresholds = c(0.50, 0.60, 0.70, 0.75)) {
  
  cat("\n", strrep("#", 90), "\n")
  cat("COUNTRY:", toupper(country), "| MODEL: BALANCED GENERALIZED GOWER\n")
  cat(strrep("#", 90), "\n")
  
  # ------------------------------------------------------------
  # 4.1 Load country data
  # ------------------------------------------------------------
  
  data_file <- paste0("data_", country, "_model.rds")
  data_country_model <- readRDS(data_file)
  
  set.seed(1234)
  split_country <- train_test_split_stratified(data_country_model)
  
  data_train <- split_country$train %>%
    mutate(.train_pos = row_number())
  
  data_test <- split_country$test
  
  cat("\nTrain proportions:\n")
  print(prop.table(table(data_train$sphus_bin)))
  
  cat("\nTest proportions:\n")
  print(prop.table(table(data_test$sphus_bin)))
  
  # ------------------------------------------------------------
  # 4.2 Prepare data
  # ------------------------------------------------------------
  
  x_test  <- data_test  %>% select(all_of(vars))
  x_train <- data_train %>% select(all_of(vars))
  
  y_test  <- as.integer(as.character(data_test$sphus_bin))
  y_train <- as.integer(as.character(data_train$sphus_bin))
  
  w_test  <- data_test$cciw_w9
  w_train <- data_train$cciw_w9
  
  n_test  <- nrow(data_test)
  n_train <- nrow(data_train)
  
  X_all <- bind_rows(x_test, x_train)
  w_all <- c(w_test, w_train)
  
  idx_test_all <- 1:n_test
  
  # ------------------------------------------------------------
  # 4.3 Identify minority and majority classes in train
  # ------------------------------------------------------------
  
  tab_train <- table(y_train)
  
  minority_class <- as.integer(names(which.min(tab_train)))
  majority_class <- as.integer(names(which.max(tab_train)))
  
  idx_minority_train <- data_train %>%
    filter(as.integer(as.character(sphus_bin)) == minority_class) %>%
    pull(.train_pos)
  
  idx_majority_train <- data_train %>%
    filter(as.integer(as.character(sphus_bin)) == majority_class) %>%
    pull(.train_pos)
  
  n_minority <- length(idx_minority_train)
  
  cat("\nMinority class:", minority_class, "| n =", n_minority, "\n")
  cat("Majority class:", majority_class, "| n =", length(idx_majority_train), "\n")
  cat("Balanced sample size per bag:", 2 * n_minority, "\n")
  
  # ------------------------------------------------------------
  # 4.4 Start total execution time
  # ------------------------------------------------------------
  
  total_time_start <- proc.time()
  
  # ------------------------------------------------------------
  # 4.5 Compute full generalized Gower distance matrix once
  # ------------------------------------------------------------
  
  cat("\nComputing full generalized Gower distance matrix...\n")
  
  D_all_raw <- robust_distances(
    data      = X_all,
    cont_vars = cont_vars,
    bin_vars  = bin_vars,
    cat_vars  = cat_vars,
    w         = w_all,
    alpha     = 0.10,
    method    = "ggower"
  )
  
  cat("Applying Euclidean correction...\n")
  
  D_all_euc <- make_euclidean(
    D_all_raw,
    w = w_all
  )$D_euc
  
  D_all_D2 <- D_all_euc^2
  
  cat("\nDistance matrix dimensions:", dim(D_all_D2), "\n")
  
  stopifnot(nrow(D_all_D2) == n_test + n_train)
  
  # ------------------------------------------------------------
  # 4.6 Balanced sub-bagging loop
  # ------------------------------------------------------------
  
  pred_matrix  <- matrix(NA, nrow = n_test, ncol = B)
  class_matrix <- matrix(NA, nrow = n_test, ncol = B)
  
  models_list  <- vector("list", B)
  idx_sub_list <- vector("list", B)
  
  dims_used       <- numeric(B)
  relgvar_reached <- numeric(B)
  
  set.seed(1234)
  
  for (b in 1:B) {
    
    cat("\n", country, "| Balanced Generalized Gower | Bag", b, "of", B, "\n")
    
    # ----------------------------------------------------------
    # Balanced sample:
    # Minority class fixed + random sample of majority class
    # ----------------------------------------------------------
    
    idx_majority_sampled <- sample(
      idx_majority_train,
      size    = n_minority,
      replace = FALSE
    )
    
    idx_balanced_train <- sort(c(idx_minority_train, idx_majority_sampled))
    
    y_b <- y_train[idx_balanced_train]
    w_b <- w_train[idx_balanced_train]
    
    cat("  Balanced sample size:", length(idx_balanced_train),
        "| Class 0:", sum(y_b == 0),
        "| Class 1:", sum(y_b == 1), "\n")
    
    # ----------------------------------------------------------
    # Translate train indices to global indices in D_all_D2
    # test = 1:n_test
    # train = n_test + train_position
    # ----------------------------------------------------------
    
    idx_balanced_global <- n_test + idx_balanced_train
    
    stopifnot(max(idx_balanced_global) <= nrow(D_all_D2))
    
    # ----------------------------------------------------------
    # Extract distance matrices
    # ----------------------------------------------------------
    
    D_train_b_D2 <- D_all_D2[idx_balanced_global, idx_balanced_global]
    class(D_train_b_D2) <- "D2"
    
    D_cross_D2 <- D_all_D2[idx_test_all, idx_balanced_global]
    class(D_cross_D2) <- "D2"
    
    # ----------------------------------------------------------
    # Fit dbglm
    # ----------------------------------------------------------
    
    model_b <- dbglm(
      D2       = D_train_b_D2,
      y        = y_b,
      family   = binomial(link = "logit"),
      method   = "rel.gvar",
      rel.gvar = gvar,
      weights  = w_b
    )
    
    dims_used[b]       <- model_b$eff.rank
    relgvar_reached[b] <- model_b$rel.gvar
    
    # ----------------------------------------------------------
    # Predict probabilities on test set
    # ----------------------------------------------------------
    
    prob_b <- predict(
      model_b,
      newdata   = D_cross_D2,
      type.pred = "response",
      type.var  = "D2"
    )
    
    pred_matrix[, b]  <- prob_b
    class_matrix[, b] <- as.integer(prob_b >= 0.50)
    
    models_list[[b]]  <- model_b
    idx_sub_list[[b]] <- idx_balanced_train
    
    cat("  Dimensions used:", dims_used[b], "\n")
    cat("  Actual rel.gvar:", round(relgvar_reached[b], 4), "\n")
  }
  
  # ------------------------------------------------------------
  # 4.7 Total execution time
  # ------------------------------------------------------------
  
  total_time <- (proc.time() - total_time_start)["elapsed"]
  
  cat("\nTotal execution time:", round(total_time, 2), "seconds\n")
  
  # ------------------------------------------------------------
  # 4.8 Aggregate predictions
  # ------------------------------------------------------------
  
  prob_balanced <- rowMeans(pred_matrix, na.rm = TRUE)
  
  cat("\nNumber of NA predictions:\n")
  print(sum(is.na(pred_matrix)))
  
  cat("\nNA predictions per bag:\n")
  print(colSums(is.na(pred_matrix)))
  
  # ------------------------------------------------------------
  # 4.9 Evaluate thresholds - UNWEIGHTED
  # ------------------------------------------------------------
  
  results_list <- lapply(thresholds, function(th) {
    evaluate_threshold_unweighted(
      probs     = prob_balanced,
      y_true    = y_test,
      threshold = th
    )
  })
  
  names(results_list) <- paste0("results_", thresholds * 100)
  
  summary_thresholds <- data.frame(
    Country     = country,
    Model       = "Balanced Generalized Gower",
    Threshold   = thresholds,
    Accuracy    = sapply(results_list, function(x) x$accuracy),
    Sensitivity = sapply(results_list, function(x) x$sensitivity),
    Specificity = sapply(results_list, function(x) x$specificity),
    PPV         = sapply(results_list, function(x) x$ppv),
    NPV         = sapply(results_list, function(x) x$npv),
    Baseline_1  = sapply(results_list, function(x) x$baseline_1),
    Baseline_Majority = sapply(results_list, function(x) x$baseline_majority)
  )
  
  cat("\nSummary thresholds - unweighted evaluation:\n")
  print(
    summary_thresholds %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  )
  
  # ------------------------------------------------------------
  # 4.10 Unweighted ROC and AUC
  # ------------------------------------------------------------
  
  roc_unweighted <- pROC::roc(
    response  = y_test,
    predictor = prob_balanced,
    levels    = c(0, 1),
    direction = "<",
    quiet     = TRUE
  )
  
  auc_unweighted <- as.numeric(pROC::auc(roc_unweighted))
  
  cat("\nUnweighted AUC:", round(auc_unweighted, 3), "\n")
  
  # ------------------------------------------------------------
  # 4.11 Dimensions by bag
  # ------------------------------------------------------------
  
  dimensions_by_bag <- data.frame(
    Country = country,
    Model = "Balanced Generalized Gower",
    Bag = 1:B,
    Dimensions = dims_used,
    Actual_rel_gvar = relgvar_reached
  )
  
  # ------------------------------------------------------------
  # 4.12 Computational summary
  # ------------------------------------------------------------
  
  computational_summary <- data.frame(
    Country = country,
    Model = "Balanced Generalized Gower",
    Distance = "balanced_generalized_gower_make_euclidean",
    Euclidean_correction = TRUE,
    Requested_rel_gvar = gvar,
    Mean_actual_rel_gvar = mean(relgvar_reached),
    Mean_dimensions = mean(dims_used),
    Min_dimensions = min(dims_used),
    Max_dimensions = max(dims_used),
    Total_time_seconds = as.numeric(total_time),
    Total_time_minutes = as.numeric(total_time) / 60,
    AUC_unweighted = auc_unweighted
  )
  
  cat("\nComputational summary:\n")
  print(
    computational_summary %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  )
  
  # ------------------------------------------------------------
  # 4.13 Final object
  # ------------------------------------------------------------
  
  result <- list(
    country = country,
    model = "Balanced Generalized Gower",
    
    models_list = models_list,
    idx_sub_list = idx_sub_list,
    
    pred_matrix = pred_matrix,
    class_matrix = class_matrix,
    prob_balanced = prob_balanced,
    
    results_list = results_list,
    summary_thresholds = summary_thresholds,
    
    roc_unweighted = roc_unweighted,
    auc_unweighted = auc_unweighted,
    
    dimensions_by_bag = dimensions_by_bag,
    computational_summary = computational_summary,
    
    gvar = gvar,
    B = B,
    vars_model = vars,
    distance = "balanced_generalized_gower_make_euclidean",
    euclidean_correction = TRUE,
    evaluation = "unweighted",
    balancing_strategy = "minority_fixed_majority_undersampled"
  )
  
  # ------------------------------------------------------------
  # 4.14 Save object
  # ------------------------------------------------------------
  
  output_file <- paste0(
    "subbagging_",
    country,
    "_balanced_ggower_70.rds"
  )
  
  saveRDS(result, output_file)
  
  cat("\nSaved object:", output_file, "\n")
  
  return(result)
}

################################################################################
############################ 5. RUN BALANCED MODELS ############################
################################################################################

# ============================================================
# LATVIA
# ============================================================

balanced_latvia_ggower <- run_country_balanced_ggower(
  country    = "latvia",
  B          = B,
  gvar       = gvar,
  thresholds = thresholds_eval
)

# ============================================================
# SPAIN
# ============================================================

balanced_spain_ggower <- run_country_balanced_ggower(
  country    = "spain",
  B          = B,
  gvar       = gvar,
  thresholds = thresholds_eval
)

 # ============================================================
# GERMANY
# ============================================================

balanced_germany_ggower <- run_country_balanced_ggower(
  country    = "germany",
  B          = B,
  gvar       = gvar,
  thresholds = thresholds_eval
)

# ============================================================
# FINLAND
# ============================================================

balanced_finland_ggower <- run_country_balanced_ggower(
  country    = "finland",
  B          = B,
  gvar       = gvar,
  thresholds = thresholds_eval
)

################################################################################
############################ 6. CHECK RESULTS ##################################
################################################################################

balanced_latvia_ggower$summary_thresholds
balanced_spain_ggower$summary_thresholds
balanced_germany_ggower$summary_thresholds
balanced_finland_ggower$summary_thresholds

balanced_latvia_ggower$auc_unweighted
balanced_spain_ggower$auc_unweighted
balanced_germany_ggower$auc_unweighted
balanced_finland_ggower$auc_unweighted

