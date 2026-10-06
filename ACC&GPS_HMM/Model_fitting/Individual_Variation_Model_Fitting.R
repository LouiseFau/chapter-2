#' -----------------------------------------------------------------------------
# Title: Individual variation model fitting ----
#' Authors : Louise Faure
#' Date : 24.07.26
#' 
#' Info : this script follow the Transition_Matrix_Preparation.R script where 
#' covariates are standardized and individuals weighted.  
#' 
#' **Questions:**
#' (i) Is an increase in HFI associated with an increase in the probability of remaining aerial, assuming a common response among individuals? 
#' remain_aerial ~ backbone + **hfi_raw_z** + elevation_100m_z + ruggedness_100m_z + **(1 | individual_id)**
#' 
#' (ii) Do individuals varies in aerial persistence and slope to hfi, given that the two are correlated ?
#' remain_aerial ~ backbone + elevation_100m_z + **hfi_within_z** + ruggedness_100m_z + **(1 + hfi_within_z | individual_id)**
#' 
#' (iii) Are the individual using high HFI areas, more likely to remain aerial ? 
#' remain_aerial ~ backbone + **hfi_between_z** + elevation_100m_z + ruggedness_100m_z + **(1 | individual_id)**
#' 
#' (iv) When an individual encounters a higher HFI than its usual environment, does it modify its probability of remaining aerial? 
#' remain_aerial ~ backbone + **hfi_within_z + hfi_between_z** + elevation_100m_z + ruggedness_100m_z + **(1 | individual_id)**
#' 
#' (v) Does this within-individual HFI response vary among individuals? (assume an association between the response and the individual)
#' remain_aerial ~ backbone + **hfi_within_z + hfi_between_z** + elevation_100m_z + ruggedness_100m_z + **(1 + hfi_within_z | individual_id)**
#' 
#' (vi) Does this within-individual HFI response vary among individuals? (assume abs of association btw the response and the individual)
#' remain_aerial ~ backbone + **hfi_within_z + hfi_between_z** + elevation_100m_z + ruggedness_100m_z + **(1 + hfi_within_z || individual_id)**
#' 
#' **Main steps:**
#' (1) fit several models for the gps dataset (only) and resume the results with a 
#' table for model comparison (stability of HFI, AIC and EDF) and select the model 
#' with the highest AIC.
#' 
#' (2) fit individual variation model for the acc dataset using the retained form
#' in 1 and compare individuals name to the one obtained in 1. 
#' 
#' (3) test whether natal experience explains individual variation in the HFI slope
#' (two-step approach on natal_q90_NT, as for the iSSF) and compare the three familiarity
#' covariates in the retained model 2 (GPS only).
#' -----------------------------------------------------------------------------

# library
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(ggplot2)
library(lme4)

# data
GE_gps_weighted <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/gps_20_weighted.rds")
GE_acc_weighted <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/acc_20_weighted.rds")

# output directory
results_directory_60 <- paste0( "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/","THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results")
dir.create(results_directory_60,recursive = TRUE,showWarnings = FALSE)


#------------------------------------------------------------------------------- STEP 1: global, within-individual and individual-varying HFI responses ----
#' **Steps**:
#' (i) prepare dataset 
#' (ii) define three models 
#' (iii) compare models

# Parameters
hfi_variable_60 <- "settlement_density_z"
backbone_terms_60 <- c("cos_diel_c","sin_diel_c","duration_z", "duration_z2")
environmental_terms_60 <- c("elevation_100m_z","ruggedness_100m_z")
confidence_level_60 <- 0.95
confidence_critical_value_60 <- stats::qnorm(1 - (1 - confidence_level_60) / 2)
contrast_probabilities_60 <- c(0.05,0.95)
estimate_intercept_slope_correlation_60 <- FALSE
control_glmer_60 <- lme4::glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))
# Define and check variables required for STEP 1 ----
control_required_variables_60 <- unique(c(
  "individual_id","remain_aerial",
  backbone_terms_60,
  environmental_terms_60,
  hfi_variable_60))


# 1.1 Prepare one common model dataset ----
gps_hfi_model_data_60 <- GE_gps_weighted %>%
  dplyr::select(dplyr::all_of(control_required_variables_60)) %>%
  tidyr::drop_na() %>%
  dplyr::filter(dplyr::if_all(-individual_id,is.finite)) %>%
  dplyr::mutate(
    individual_id = droplevels(factor(individual_id)),
    hfi_raw_z = .data[[hfi_variable_60]]) %>%
  dplyr::group_by(individual_id) %>%
  dplyr::mutate(
    hfi_between_z = mean(hfi_raw_z),
    hfi_within_z = hfi_raw_z - hfi_between_z) %>%
  dplyr::ungroup()

# 1.2 Define the six human-pressure models ----
fixed_adjustment_terms_60 <- paste(
  c(backbone_terms_60,environmental_terms_60),
  collapse = " + ")

hfi_formulas_60 <- list(
  `1. Common global HFI slope` = stats::as.formula(
    paste("remain_aerial ~",fixed_adjustment_terms_60,
          "+ hfi_raw_z + (1 | individual_id)")),
  `2. Correlated individual global HFI slopes` = stats::as.formula(
    paste("remain_aerial ~",fixed_adjustment_terms_60,
          "+ hfi_raw_z + (1 + hfi_raw_z | individual_id)")),
  `3. Between-individual HFI effect` = stats::as.formula(
    paste("remain_aerial ~",fixed_adjustment_terms_60,
          "+ hfi_between_z + (1 | individual_id)")),
  `4. Common within-individual HFI slope` = stats::as.formula(
    paste( "remain_aerial ~",fixed_adjustment_terms_60,
           "+ hfi_within_z + hfi_between_z + (1 | individual_id)")),
  `5. Correlated individual within-HFI slopes` = stats::as.formula(
    paste("remain_aerial ~",fixed_adjustment_terms_60,
          "+ hfi_within_z + hfi_between_z +",
          "(1 + hfi_within_z | individual_id)")),
  `6. Uncorrelated individual within-HFI slopes` = stats::as.formula(
    paste("remain_aerial ~",fixed_adjustment_terms_60,
          "+ hfi_within_z + hfi_between_z +",
          "(1 + hfi_within_z || individual_id)")))

# 1.3 Fit the six unweighted binomial GLMMs ----
fit_hfi_glmer_60 <- function(model_formula) {
  lme4::glmer(
    formula = model_formula,
    data = gps_hfi_model_data_60,
    family = stats::binomial(link = "logit"),
    nAGQ = 1,
    control = control_glmer_60)}

hfi_models_60 <- stats::setNames(lapply(hfi_formulas_60,fit_hfi_glmer_60),names(hfi_formulas_60))

# Store models using their exact list names ----
global_common_hfi_model_ml_60 <- hfi_models_60[["1. Common global HFI slope"]]
global_correlated_hfi_model_ml_60 <- hfi_models_60[["2. Correlated individual global HFI slopes"]]
between_individual_hfi_model_ml_60 <- hfi_models_60[["3. Between-individual HFI effect"]]
within_common_hfi_model_ml_60 <- hfi_models_60[["4. Common within-individual HFI slope"]]
within_correlated_hfi_model_ml_60 <- hfi_models_60[["5. Correlated individual within-HFI slopes"]]
within_uncorrelated_hfi_model_ml_60 <- hfi_models_60[["6. Uncorrelated individual within-HFI slopes"]]

# 1.4 Compare the relevant random-effect structures ----
extract_lrt_p_value_60 <- function(lrt_table) {
  p_column <- grep("^Pr\\(",names(lrt_table),value = TRUE)
  if(length(p_column) != 1L || nrow(lrt_table) < 2L) return(NA_real_)
  as.numeric(lrt_table[2,p_column])
}

global_slope_lrt_60 <- stats::anova(
  global_common_hfi_model_ml_60,
  global_correlated_hfi_model_ml_60,
  test = "Chisq"
)

within_correlated_slope_lrt_60 <- stats::anova(
  within_common_hfi_model_ml_60,
  within_correlated_hfi_model_ml_60,
  test = "Chisq"
)

within_uncorrelated_slope_lrt_60 <- stats::anova(
  within_common_hfi_model_ml_60,
  within_uncorrelated_hfi_model_ml_60,
  test = "Chisq"
)

within_correlation_lrt_60 <- stats::anova(
  within_uncorrelated_hfi_model_ml_60,
  within_correlated_hfi_model_ml_60,
  test = "Chisq"
)

global_slope_lrt_p_60 <- extract_lrt_p_value_60(global_slope_lrt_60)
within_correlated_slope_lrt_p_60 <- extract_lrt_p_value_60(within_correlated_slope_lrt_60)
within_uncorrelated_slope_lrt_p_60 <- extract_lrt_p_value_60(within_uncorrelated_slope_lrt_60)
within_correlation_lrt_p_60 <- extract_lrt_p_value_60(within_correlation_lrt_60)

# 1.5 Define Q05-Q95 contrasts for global, between and within HFI ----
calculate_q95_q05_difference_60 <- function(x) {
  quantiles <- stats::quantile(
    x,probs = contrast_probabilities_60,
    na.rm = TRUE,names = FALSE)
  as.numeric(quantiles[2] - quantiles[1])
}

between_hfi_values_60 <- gps_hfi_model_data_60 %>%
  dplyr::distinct(individual_id,hfi_between_z)

hfi_contrast_differences_60 <- c(
  global = calculate_q95_q05_difference_60(gps_hfi_model_data_60$hfi_raw_z),
  between = calculate_q95_q05_difference_60(between_hfi_values_60$hfi_between_z),
  within = calculate_q95_q05_difference_60(gps_hfi_model_data_60$hfi_within_z)
)

# 1.6 Define each model question and focal parameter ----
hfi_model_metadata_60 <- tibble::tibble(
  model = names(hfi_models_60),
  question = c(
    "Is HFI associated with aerial persistence, assuming one common response?",
    "Do individuals differ in their global settlements slope, and is this slope correlated with baseline aerial persistence?",
    "Are individuals using areas with higher mean settlements more likely to remain aerial?",
    "Does an individual respond when encountered settlement differs from its usual HFI?",
    "Does the within-individual settlement response vary among individuals and correlate with baseline aerial persistence?",
    "Does the within-individual settlement response vary among individuals without intercept-slope correlation?"
  ),
  focal_term = c(
    "hfi_raw_z","hfi_raw_z","hfi_between_z",
    "hfi_within_z","hfi_within_z","hfi_within_z"
  ),
  random_slope_term = c(
    NA_character_,"hfi_raw_z",NA_character_,
    NA_character_,"hfi_within_z","hfi_within_z"
  ),
  focal_level = c(
    "global","global","between","within","within","within"
  ),
  random_slope_LRT_p = c(
    NA_real_,global_slope_lrt_p_60,NA_real_,NA_real_,
    within_correlated_slope_lrt_p_60,
    within_uncorrelated_slope_lrt_p_60
  ),
  correlation_LRT_p = c(
    NA_real_,NA_real_,NA_real_,NA_real_,
    within_correlation_lrt_p_60,NA_real_
  ),
  formula = unname(vapply(
    hfi_formulas_60,
    function(x) paste(deparse(x,width.cutoff = 500),collapse = " "),
    character(1)
  ))
)

# 1.7 Extract fixed effects and random-effect statistics ----
extract_model_convergence_60 <- function(fitted_model) {
  convergence_message <- fitted_model@optinfo$conv$lme4$messages
  optimizer_code <- fitted_model@optinfo$conv$opt
  is.null(convergence_message) &&
    (is.null(optimizer_code) || all(optimizer_code == 0))
}

extract_random_intercept_sd_60 <- function(fitted_model) {
  variance_table <- as.data.frame(lme4::VarCorr(fitted_model))
  value <- variance_table$sdcor[
    grepl("^individual_id",variance_table$grp) &
      variance_table$var1 == "(Intercept)" &
      is.na(variance_table$var2)
  ]
  if(length(value) == 0L) NA_real_ else value[[1]]
}

extract_random_slope_sd_60 <- function(fitted_model,random_slope_term) {
  if(is.na(random_slope_term)) return(NA_real_)
  variance_table <- as.data.frame(lme4::VarCorr(fitted_model))
  value <- variance_table$sdcor[
    grepl("^individual_id",variance_table$grp) &
      variance_table$var1 == random_slope_term &
      is.na(variance_table$var2)
  ]
  if(length(value) == 0L) NA_real_ else value[[1]]
}

extract_intercept_slope_correlation_60 <- function(fitted_model,random_slope_term) {
  if(is.na(random_slope_term)) return(NA_real_)
  variance_table <- as.data.frame(lme4::VarCorr(fitted_model))
  value <- variance_table$sdcor[
    grepl("^individual_id",variance_table$grp) &
      !is.na(variance_table$var2) &
      ((variance_table$var1 == "(Intercept)" &
          variance_table$var2 == random_slope_term) |
         (variance_table$var2 == "(Intercept)" &
            variance_table$var1 == random_slope_term))
  ]
  if(length(value) == 0L) NA_real_ else value[[1]]
}

hfi_model_results_60 <- dplyr::bind_rows(
  lapply(seq_len(nrow(hfi_model_metadata_60)),function(i) {
    model_information <- hfi_model_metadata_60[i,]
    model_name <- model_information$model[[1]]
    fitted_model <- hfi_models_60[[model_name]]
    coefficient_table <- summary(fitted_model)$coefficients
    focal_term <- model_information$focal_term[[1]]
    focal_level <- model_information$focal_level[[1]]
    random_slope_term <- model_information$random_slope_term[[1]]
    
    if(!focal_term %in% rownames(coefficient_table)) {
      stop("Term ",focal_term," is absent from model ",model_name,".")
    }
    
    estimate <- coefficient_table[focal_term,"Estimate"]
    standard_error <- coefficient_table[focal_term,"Std. Error"]
    confidence_low <- estimate - confidence_critical_value_60 * standard_error
    confidence_high <- estimate + confidence_critical_value_60 * standard_error
    contrast_difference <- hfi_contrast_differences_60[[focal_level]]
    model_loglik <- stats::logLik(fitted_model)
    
    tibble::tibble(
      model = model_name,
      question = model_information$question[[1]],
      formula = model_information$formula[[1]],
      focal_term = focal_term,
      n_observations = stats::nobs(fitted_model),
      n_individuals = dplyr::n_distinct(stats::model.frame(fitted_model)$individual_id),
      model_df = attr(model_loglik,"df"),
      AIC = stats::AIC(fitted_model),
      hfi_coefficient = estimate,
      hfi_confidence_low = confidence_low,
      hfi_confidence_high = confidence_high,
      Q95_Q05_log_odds = estimate * contrast_difference,
      Q95_Q05_log_odds_CI_low = confidence_low * contrast_difference,
      Q95_Q05_log_odds_CI_high = confidence_high * contrast_difference,
      random_intercept_SD = extract_random_intercept_sd_60(fitted_model),
      random_slope_SD = extract_random_slope_sd_60(fitted_model,random_slope_term),
      intercept_slope_correlation =
        extract_intercept_slope_correlation_60(fitted_model,random_slope_term),
      random_slope_LRT_p = model_information$random_slope_LRT_p[[1]],
      correlation_LRT_p = model_information$correlation_LRT_p[[1]],
      singular = lme4::isSingular(fitted_model,tol = 1e-4),
      converged = extract_model_convergence_60(fitted_model)
    )
  })
) %>%
  dplyr::mutate(delta_AIC = AIC - min(AIC))

# 1.8 Prepare the final comparison table ----
format_p_value_60 <- function(x) {
  ifelse(is.na(x),"--",format.pval(x,digits = 3,eps = 0.001))
}

extract_compact_model_statistics_60 <- function(i) {
  model_result <- hfi_model_results_60[i,]
  
  tibble::tibble(
    Statistic = c(
      "Model formula",
      "Focal HFI term",
      "Number of observations",
      "Number of individuals",
      "AIC",
      "Delta AIC",
      "Model df",
      "HFI coefficient",
      "HFI coefficient 95% CI",
      "Q05-Q95 HFI log-odds contrast",
      "Q05-Q95 HFI log-odds 95% CI",
      "Random-intercept SD",
      "Random-slope SD",
      "Intercept-slope correlation",
      "Random-slope LRT p-value",
      "Correlation LRT p-value",
      "Singular fit",
      "Converged"
    ),
    Value = c(
      model_result$formula,
      model_result$focal_term,
      as.character(model_result$n_observations),
      as.character(model_result$n_individuals),
      sprintf("%.1f",model_result$AIC),
      sprintf("%.1f",model_result$delta_AIC),
      sprintf("%.0f",model_result$model_df),
      sprintf("%.3f",model_result$hfi_coefficient),
      sprintf(
        "[%.3f, %.3f]",
        model_result$hfi_confidence_low,
        model_result$hfi_confidence_high
      ),
      sprintf("%.3f",model_result$Q95_Q05_log_odds),
      sprintf(
        "[%.3f, %.3f]",
        model_result$Q95_Q05_log_odds_CI_low,
        model_result$Q95_Q05_log_odds_CI_high
      ),
      ifelse(
        is.na(model_result$random_intercept_SD),
        "--",
        sprintf("%.3f",model_result$random_intercept_SD)
      ),
      ifelse(
        is.na(model_result$random_slope_SD),
        "--",
        sprintf("%.3f",model_result$random_slope_SD)
      ),
      ifelse(
        is.na(model_result$intercept_slope_correlation),
        "--",
        sprintf("%.3f",model_result$intercept_slope_correlation)
      ),
      format_p_value_60(model_result$random_slope_LRT_p),
      format_p_value_60(model_result$correlation_LRT_p),
      as.character(model_result$singular),
      as.character(model_result$converged)
    )
  )
}

model_statistics_list_60 <- stats::setNames(
  lapply(
    seq_len(nrow(hfi_model_results_60)),
    extract_compact_model_statistics_60
  ),
  hfi_model_results_60$model
)

hfi_model_comparison_table_60 <- dplyr::bind_rows(
  model_statistics_list_60,
  .id = "Model"
) %>%
  tidyr::pivot_wider(
    names_from = Model,
    values_from = Value
  )

# 1.9 Format the old-style HTML table ----
hfi_model_comparison_gt_60 <- hfi_model_comparison_table_60 %>%
  gt::gt() %>%
  gt::tab_header(
    title = gt::md("**GPS HFI individual-effect models**"),
    subtitle = paste0(
      "Binomial GLMMs; backbone + elevation + ruggedness; ",
      "HFI mean within 1000 m"
    )
  ) %>%
  gt::cols_label(
    Statistic = "Statistic"
  ) %>%
  gt::cols_align(
    align = "left",
    columns = gt::everything()
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_body(columns = Statistic)
  ) %>%
  gt::tab_style(
    style = gt::cell_text(
      font = "monospace",
      size = gt::px(11)
    ),
    locations = gt::cells_body(
      columns = -Statistic,
      rows = Statistic == "Model formula"
    )
  ) %>%
  gt::tab_style(
    style = gt::cell_fill(color = "grey95"),
    locations = gt::cells_body(
      rows = Statistic == "Model formula"
    )
  ) %>%
  gt::cols_width(
    Statistic ~ gt::px(240),
    gt::everything() ~ gt::px(340)
  ) %>%
  gt::tab_options(
    table.width = gt::px(2300),
    table.font.size = gt::px(12),
    heading.title.font.size = gt::px(18),
    heading.subtitle.font.size = gt::px(13),
    column_labels.font.weight = "bold",
    data_row.padding = gt::px(6)
  ) %>%
  gt::tab_source_note(
    source_note = gt::md(
      paste0(
        "**Random-slope tests:** Models 1–2 test variation in the global HFI slope; ",
        "Models 4–5 and 4–6 test variation in the within-individual HFI slope. ",
        "**Correlation test:** Model 6 versus Model 5 tests the intercept–slope correlation."
      )
    )
  )

# 1.10 Export and open the HTML table ----
hfi_model_comparison_file_60 <- file.path(
  results_directory_60,
  "GPS_HFI_six_individual_effect_models.html"
)

gt::gtsave(
  data = hfi_model_comparison_gt_60,
  filename = basename(hfi_model_comparison_file_60),
  path = dirname(hfi_model_comparison_file_60)
)

browseURL(normalizePath(hfi_model_comparison_file_60))

#------------------------------------------------------------------------------- STEP 2: compare GPS and ACC individual within-HFI responses ----
#' **Steps:**
#' (i) prepare the ACC aerial-versus-terrestrial dataset;
#' (ii) fit the retained uncorrelated individual within-HFI slope model;
#' (iii) extract and compare GPS and ACC individual HFI slopes;
#' (iv) identify individuals with positive estimated slopes in both datasets;
#' (v) describe feeding and resting behaviour among ACC individuals with
#'     positive aerial-versus-terrestrial HFI slopes.

# 2.1 Prepare the ACC aerial-versus-terrestrial dataset ----
acc_required_variables_60 <- unique(c(
  "individual_id",
  "transition_destination",
  "remain_aerial",
  backbone_terms_60,
  environmental_terms_60,
  hfi_variable_60
))

acc_hfi_model_data_60 <- GE_acc_weighted %>%
  dplyr::select(dplyr::all_of(acc_required_variables_60)) %>%
  tidyr::drop_na() %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(
        setdiff(
          acc_required_variables_60,
          c("individual_id","transition_destination")
        )
      ),
      is.finite
    )
  ) %>%
  dplyr::mutate(
    individual_id = droplevels(factor(individual_id)),
    hfi_raw_z = .data[[hfi_variable_60]]
  )

# 2.2 Fit model 2 to the ACC dataset ----
acc_global_correlated_formula_60 <-
  hfi_formulas_60[["2. Correlated individual global HFI slopes"]]

acc_global_correlated_hfi_model_ml_60 <- lme4::glmer(
  formula = acc_global_correlated_formula_60,
  data = acc_hfi_model_data_60,
  family = stats::binomial(link = "logit"),
  nAGQ = 1,
  control = control_glmer_60
)

# 2.3 Extract GPS and ACC individual global HFI slopes ----
extract_individual_global_hfi_slopes_60 <- function(fitted_model,dataset_name) {
  population_slope <- unname(lme4::fixef(fitted_model)[["hfi_raw_z"]])
  individual_random_effects <- lme4::ranef(fitted_model)$individual_id
  
  if(length(population_slope) != 1L || !is.finite(population_slope)) {
    stop("The fixed global HFI slope was not found for ",dataset_name,".")
  }
  
  if(!"hfi_raw_z" %in% colnames(individual_random_effects)) {
    stop("The individual global HFI random slope was not found for ",dataset_name,".")
  }
  
  tibble::tibble(
    dataset = dataset_name,
    individual_id = rownames(individual_random_effects),
    population_hfi_slope = population_slope,
    individual_slope_deviation = individual_random_effects[,"hfi_raw_z"],
    individual_hfi_slope = population_slope + individual_slope_deviation,
    positive_estimated_slope = individual_hfi_slope > 0
  ) %>%
    dplyr::arrange(dplyr::desc(individual_hfi_slope))
}

gps_global_correlated_hfi_model_ml_60 <-
  hfi_models_60[["2. Correlated individual global HFI slopes"]]

gps_individual_hfi_slopes_60 <- extract_individual_global_hfi_slopes_60(
  fitted_model = gps_global_correlated_hfi_model_ml_60,
  dataset_name = "GPS"
)

acc_individual_hfi_slopes_60 <- extract_individual_global_hfi_slopes_60(
  fitted_model = acc_global_correlated_hfi_model_ml_60,
  dataset_name = "ACC"
)

# 2.4 Compare individual names and slope signs between GPS and ACC ----
individual_hfi_slope_comparison_60 <- dplyr::inner_join(
  gps_individual_hfi_slopes_60 %>%
    dplyr::select(
      individual_id,
      gps_hfi_slope = individual_hfi_slope
    ),
  acc_individual_hfi_slopes_60 %>%
    dplyr::select(
      individual_id,
      acc_hfi_slope = individual_hfi_slope
    ),
  by = "individual_id"
) %>%
  dplyr::mutate(
    gps_positive_slope = gps_hfi_slope > 0,
    acc_positive_slope = acc_hfi_slope > 0,
    gps_negative_slope = gps_hfi_slope < 0,
    acc_negative_slope = acc_hfi_slope < 0,
    slope_comparison = dplyr::case_when(
      gps_negative_slope & acc_negative_slope ~ "Negative in both",
      gps_negative_slope & acc_positive_slope ~ "Negative in GPS only",
      gps_positive_slope & acc_negative_slope ~ "Negative in ACC only",
      gps_positive_slope & acc_positive_slope ~ "Positive in both",
      TRUE ~ "Slope equal to zero"
    )
  )

# Retain and display only individuals with at least one negative slope ----
negative_individual_hfi_slopes_60 <- individual_hfi_slope_comparison_60 %>%
  dplyr::filter(
    gps_negative_slope | acc_negative_slope
  ) %>%
  dplyr::select(
    individual_id,
    gps_hfi_slope,
    acc_hfi_slope,
    slope_comparison
  ) %>%
  dplyr::arrange(
    factor(
      slope_comparison,
      levels = c(
        "Negative in both",
        "Negative in GPS only",
        "Negative in ACC only",
        "Slope equal to zero"
      )
    ),
    gps_hfi_slope,
    acc_hfi_slope
  )

print(
  negative_individual_hfi_slopes_60,
  n = Inf,
  width = Inf
)


# for 20 minutes dataset 
# individual_id            gps_hfi_slope acc_hfi_slope slope_comparison    
# Flüela2 21 (eobs 7043)         -0.0211       0.00936 Negative in GPS only
# Dischma1 19 (eobs 7006)         0.0375      -0.0107  Negative in ACC only
# Sinestra1 19 (eobs 7003)        0.0515      -0.00986 Negative in ACC only
# Ettenberg22 (eobs 10539)        0.0804      -0.00718 Negative in ACC only
# Fahrntal19 (eobs 7014)          0.0854      -0.0321  Negative in ACC only

# Summary of positive and negative slope agreement ----
individual_hfi_slope_summary_60 <- individual_hfi_slope_comparison_60 %>%
  dplyr::summarise(
    n_shared_individuals = dplyr::n(),
    n_positive_in_both = sum(
      gps_positive_slope & acc_positive_slope
    ),
    n_negative_in_both = sum(
      gps_negative_slope & acc_negative_slope
    ),
    n_negative_in_GPS_only = sum(
      gps_negative_slope & acc_positive_slope
    ),
    n_negative_in_ACC_only = sum(
      gps_positive_slope & acc_negative_slope
    ),
    n_negative_in_GPS = sum(gps_negative_slope),
    n_negative_in_ACC = sum(acc_negative_slope)
  )

print(individual_hfi_slope_summary_60)
# n_shared_individuals n_positive_in_both n_negative_in_both n_negative_in_GPS_only n_negative_in_ACC_only n_negative_in_GPS
#    55                 50                  0                      1                      4                 1


#------------------------------------------------------------------------------- STEP 3: natal experience and individual variation in aerial persistence (GPS only) ----
#' **Steps:**
#' (i) two-step approach, as for the iSSF: one binomial GLM per individual (backbone +
#'     elevation + ruggedness + HFI) gives the individual HFI slope on the probability of
#'     remaining aerial; the individual slopes are then regressed on the 90th percentile of
#'     settlement density in the natal territory (natal_q90_NT), with territory as a random
#'     effect (model A without predictor, model B with natal_q90_NT);
#' (ii) the three familiarity covariates (position_NT, above_natal_NT, excess_NT) are defined
#'     for each location and not for each individual: they cannot be predictors of the
#'     individual slopes. They are added one at a time to the retained model 2
#'     (hfi_raw_z + (1 + hfi_raw_z | individual_id)) and compared with the same covariate built
#'     on a reference common to all individuals (natal, common, both).
#' **Expected under the hypothesis:** the HFI slope is positive (eagles remain aerial over
#' settlements). Individuals from more built natal territories should have a less positive
#' slope (negative natal coefficient in step (i)), and unfamiliar settlement levels should
#' increase aerial persistence (positive familiarity coefficients in step (ii)).
#' Requires the objects of STEP 1 (parameters, 1.1, 1.2 and the extraction functions of 1.7).

# Parameters
natal_descriptor_60 <- "natal_q90_NT"
familiarity_covariates_60 <- c("position_NT","above_natal_NT","excess_NT")
common_excess_probability_60 <- 0.95
individual_terms_60 <- c(backbone_terms_60,environmental_terms_60,"hfi_raw_z")

# 3.1 Prepare the GPS dataset with the familiarity covariates ----
controle_missing_familiarity_60 <- setdiff(c(paste0(familiarity_covariates_60,"_z"),"settlement_density"),names(GE_gps_weighted))
if(length(controle_missing_familiarity_60) > 0L) stop("Missing columns: ",paste(controle_missing_familiarity_60,collapse = ", "))

gps_familiarity_data_60 <- GE_gps_weighted %>%
  dplyr::select(dplyr::all_of(unique(c(control_required_variables_60,"settlement_density",paste0(familiarity_covariates_60,"_z"))))) %>%
  tidyr::drop_na() %>%
  dplyr::filter(dplyr::if_all(-individual_id,is.finite)) %>%
  dplyr::mutate(
    individual_id = droplevels(factor(individual_id)),
    hfi_raw_z = .data[[hfi_variable_60]],
    key = gsub("\\s+"," ",trimws(as.character(individual_id))))

# 3.2 Natal descriptor per individual (same values as in the iSSF) ----
# 3.2 Natal descriptor per individual (same polygons, raster and definition as in the iSSF) ----
#' Computed here from the natal polygons, so that STEP 3 does not depend on the columns
#' stored in the iSSF annotated dataset. territory_id is built exactly as in Extract_covariates.R.
settlement_density_raster_60 <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/settlement_density_1km2_100m.tif")
terra::crs(settlement_density_raster_60) <- "EPSG:3035"
id_lookup_path_60 <- "/Users/louisefaure/Desktop/dossier sans titre/donnees aigles gps burst/gps_bursts_raw_move2.rds"
natal_polygon_path_60 <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/natal_polygons/natal_polygon_edges.csv"

id_lookup_60 <- move2::mt_track_data(readRDS(id_lookup_path_60)) %>%
  dplyr::transmute(id = as.character(individual_id),key = gsub("\\s+"," ",trimws(as.character(individual_local_identifier))))

natal_polygons_60 <- readr::read_csv(natal_polygon_path_60,show_col_types = FALSE) %>%
  dplyr::transmute(x,y,id = as.character(id)) %>%
  dplyr::distinct(id,x,y) %>%
  sf::st_as_sf(coords = c("x","y"),crs = 4326) %>%
  sf::st_transform(3035) %>%
  dplyr::group_by(id) %>%
  dplyr::filter(dplyr::n() >= 3) %>%
  dplyr::summarise(do_union = FALSE) %>%
  sf::st_cast("LINESTRING") %>%
  sf::st_cast("POLYGON") %>%
  sf::st_make_valid() %>%
  sf::st_collection_extract("POLYGON") %>%
  dplyr::inner_join(id_lookup_60,by = "id") %>%
  dplyr::mutate(territory_id = sapply(sf::st_equals_exact(.,.,par = 1),min))

weighted_quantile_60 <- function(v,w,p) {o <- order(v); v[o][which(cumsum(w[o]) / sum(w) >= p)[1]]}

natal_descriptors_60 <- terra::extract(settlement_density_raster_60,terra::vect(natal_polygons_60),exact = TRUE) %>%
  stats::setNames(c("row","value","fraction")) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(
    key = natal_polygons_60$key[row],
    territory_id = natal_polygons_60$territory_id[row],
    value = round(value,6)) %>%
  dplyr::group_by(key,territory_id) %>%
  dplyr::summarise(
    natal_q90_NT = weighted_quantile_60(value,fraction,0.90),
    natal_prop_built_NT = stats::weighted.mean(value > 0,fraction),
    natal_mean_NT = stats::weighted.mean(value,fraction),
    .groups = "drop") %>%
  dplyr::mutate(natal_value = .data[[natal_descriptor_60]])

# 3.3 Step 1: one binomial GLM per individual ----
individual_formula_60 <- stats::as.formula(paste("remain_aerial ~",paste(individual_terms_60,collapse = " + ")))

individual_data_60 <- split(gps_familiarity_data_60,gps_familiarity_data_60$individual_id,drop = TRUE)
individual_models_60 <- lapply(individual_data_60,function(d) {
  tryCatch(stats::glm(individual_formula_60,family = stats::binomial(link = "logit"),data = d),error = function(e) NULL)})

individual_slopes_60 <- dplyr::bind_rows(lapply(names(individual_models_60),function(ind) {
  fitted_model <- individual_models_60[[ind]]
  d <- individual_data_60[[ind]]
  coefficient_table <- if(is.null(fitted_model)) NULL else summary(fitted_model)$coefficients
  has_hfi <- !is.null(coefficient_table) && "hfi_raw_z" %in% rownames(coefficient_table)
  tibble::tibble(
    individual_id = ind,
    key = d$key[1],
    n_transitions = nrow(d),
    n_landings = sum(d$remain_aerial == 0L),
    converged = has_hfi && isTRUE(fitted_model$converged),
    slope = if(has_hfi) unname(coefficient_table["hfi_raw_z","Estimate"]) else NA_real_,
    slope_se = if(has_hfi) unname(coefficient_table["hfi_raw_z","Std. Error"]) else NA_real_)
})) %>%
  dplyr::left_join(natal_descriptors_60,by = "key")

controle_individual_fits_60 <- individual_slopes_60 %>%
  dplyr::summarise(
    n_individuals = dplyr::n(),
    n_converged = sum(converged),
    n_missing_natal = sum(is.na(natal_value)),
    median_slope_se = stats::median(slope_se,na.rm = TRUE),
    max_slope_se = max(slope_se,na.rm = TRUE))

controle_unstable_individuals_60 <- individual_slopes_60 %>%
  dplyr::filter(!converged | slope_se > 5 * stats::median(slope_se,na.rm = TRUE)) %>%
  dplyr::arrange(dplyr::desc(slope_se))

# 3.4 Step 2: regress the individual slopes on the natal descriptor ----
#' slope_i = X_i beta + u_territory + u_individual + e_i, with e_i ~ N(0, slope_se_i^2) known.
#' Maximum likelihood fit, equivalent to a random-effects meta-regression (metafor::rma.mv, ML).
two_step_data_60 <- individual_slopes_60 %>%
  dplyr::filter(converged,is.finite(slope),is.finite(slope_se),!is.na(natal_value),!is.na(territory_id)) %>%
  dplyr::mutate(
    natal_value_z = as.numeric(scale(natal_value)),
    territory_id = factor(territory_id))

fit_meta_regression_60 <- function(moderator_formula,data) {
  X <- stats::model.matrix(moderator_formula,data)
  y <- data$slope
  Z <- stats::model.matrix(~ 0 + territory_id,data)
  sigma_matrix <- function(log_tau2) diag(data$slope_se^2 + exp(log_tau2[2]),nrow(data)) + exp(log_tau2[1]) * tcrossprod(Z)
  gls <- function(S) {
    XtSi <- crossprod(X,solve(S))
    cov_beta <- solve(XtSi %*% X)
    list(beta = stats::setNames(drop(cov_beta %*% XtSi %*% y),colnames(X)),cov_beta = cov_beta)}
  negative_loglik <- function(log_tau2) {
    S <- sigma_matrix(log_tau2)
    r <- y - drop(X %*% gls(S)$beta)
    0.5 * (as.numeric(determinant(S,logarithm = TRUE)$modulus) + drop(crossprod(r,solve(S,r))) + length(y) * log(2 * pi))}
  opt <- stats::optim(c(log(0.1),log(0.1)),negative_loglik,method = "L-BFGS-B",lower = c(-20,-20),upper = c(5,5))
  g <- gls(sigma_matrix(opt$par))
  list(
    beta = g$beta,
    se = stats::setNames(sqrt(diag(g$cov_beta)),colnames(X)),
    tau2 = stats::setNames(exp(opt$par),c("territory","individual")),
    logLik = -opt$value,
    AIC = 2 * opt$value + 2 * (ncol(X) + 2),
    convergence = opt$convergence)
}

two_step_models_60 <- list(
  A_no_predictor = fit_meta_regression_60(~ 1,two_step_data_60),
  B_natal = fit_meta_regression_60(~ natal_value_z,two_step_data_60))

# heterogeneity (Cochran's Q): do individual slopes vary more than their sampling error?
fixed_effect_A_60 <- stats::lm(slope ~ 1,data = two_step_data_60,weights = 1 / slope_se^2)
heterogeneity_Q_60 <- sum(stats::weights(fixed_effect_A_60) * stats::residuals(fixed_effect_A_60)^2)

# 3.5 Two-step results ----
#' Expected result 1: heterogeneity_p < 0.05 (individuals differ more than their imprecision predicts)
#' Expected result 2: natal_coefficient < 0 with a CI excluding 0 (less aerial persistence over
#' settlements for individuals from more built natal territories)
#' Expected result 3: variance_explained_percent clearly above 0, delta_AIC_B_minus_A < -2
model_A_60 <- two_step_models_60$A_no_predictor
model_B_60 <- two_step_models_60$B_natal

two_step_results_60 <- tibble::tibble(
  natal_descriptor = natal_descriptor_60,
  n_individuals = nrow(two_step_data_60),
  n_territories = dplyr::n_distinct(two_step_data_60$territory_id),
  optimizer_ok = model_A_60$convergence == 0 & model_B_60$convergence == 0,
  # expected result 1
  mean_slope = unname(model_A_60$beta["(Intercept)"]),
  between_individual_sd = sqrt(sum(model_A_60$tau2)),
  prop_individuals_positive_slope = 1 - stats::pnorm(0,mean = mean_slope,sd = between_individual_sd),
  heterogeneity_Q = heterogeneity_Q_60,
  heterogeneity_p = stats::pchisq(heterogeneity_Q_60,df = nrow(two_step_data_60) - 1,lower.tail = FALSE),
  # expected result 2
  natal_coefficient = unname(model_B_60$beta["natal_value_z"]),
  natal_CI_95 = sprintf("[%.3f; %.3f]",natal_coefficient - confidence_critical_value_60 * model_B_60$se["natal_value_z"],natal_coefficient + confidence_critical_value_60 * model_B_60$se["natal_value_z"]),
  natal_p = 2 * stats::pnorm(-abs(natal_coefficient / model_B_60$se["natal_value_z"])),
  # expected result 3
  variance_explained_percent = 100 * (1 - sum(model_B_60$tau2) / sum(model_A_60$tau2)),
  delta_AIC_B_minus_A = model_B_60$AIC - model_A_60$AIC,
  LRT_p = stats::pchisq(2 * (model_B_60$logLik - model_A_60$logLik),df = 1,lower.tail = FALSE))

print(two_step_results_60,width = Inf)

# 3.6 Familiarity covariates in the retained population model (model 2) ----
#' Each familiarity covariate is compared with the same covariate built on a reference common
#' to all individuals (all locations of the dataset): mid-rank of settlement density for the
#' position, the same rank above the median for above_natal, density above the 95th percentile
#' for the excess. The natal version is supported only if the "_both" model has a lower AIC than
#' the "_common" model and the natal coefficient remains distinct from 0.
reference_density_60 <- sort(round(gps_familiarity_data_60$settlement_density,6))
common_threshold_60 <- stats::quantile(gps_familiarity_data_60$settlement_density,common_excess_probability_60,names = FALSE)

gps_familiarity_data_60 <- gps_familiarity_data_60 %>%
  dplyr::mutate(
    position_common = (findInterval(round(settlement_density,6),reference_density_60) + findInterval(round(settlement_density,6),reference_density_60,left.open = TRUE)) / (2 * length(reference_density_60)),
    above_common = 2 * pmax(position_common - 0.5,0),
    excess_common = pmax(settlement_density - common_threshold_60,0)) %>%
  dplyr::mutate(dplyr::across(c(position_common,above_common,excess_common),~ as.numeric(scale(.x)),.names = "{.col}_z"))

base_familiarity_formula_60 <- hfi_formulas_60[["2. Correlated individual global HFI slopes"]]

familiarity_formulas_60 <- list(
  base_model = base_familiarity_formula_60,
  position_natal = stats::update(base_familiarity_formula_60,. ~ . + position_NT_z),
  position_common = stats::update(base_familiarity_formula_60,. ~ . + position_common_z),
  position_both = stats::update(base_familiarity_formula_60,. ~ . + position_common_z + position_NT_z),
  above_natal = stats::update(base_familiarity_formula_60,. ~ . + above_natal_NT_z),
  above_common = stats::update(base_familiarity_formula_60,. ~ . + above_common_z),
  above_both = stats::update(base_familiarity_formula_60,. ~ . + above_common_z + above_natal_NT_z),
  excess_natal = stats::update(base_familiarity_formula_60,. ~ . + excess_NT_z),
  excess_common = stats::update(base_familiarity_formula_60,. ~ . + excess_common_z),
  excess_both = stats::update(base_familiarity_formula_60,. ~ . + excess_common_z + excess_NT_z))

fit_familiarity_glmer_60 <- function(model_formula) {
  lme4::glmer(
    formula = model_formula,
    data = gps_familiarity_data_60,
    family = stats::binomial(link = "logit"),
    nAGQ = 1,
    control = control_glmer_60)}

familiarity_models_60 <- stats::setNames(lapply(familiarity_formulas_60,fit_familiarity_glmer_60),names(familiarity_formulas_60))

# 3.7 Compare the familiarity models ----
familiarity_natal_terms_60 <- c(position = "position_NT_z",above = "above_natal_NT_z",excess = "excess_NT_z")

familiarity_model_comparison_60 <- dplyr::bind_rows(lapply(names(familiarity_models_60),function(model_name) {
  fitted_model <- familiarity_models_60[[model_name]]
  coefficient_table <- summary(fitted_model)$coefficients
  covariate <- if(model_name == "base_model") "none" else sub("_.*","",model_name)
  natal_term <- if(covariate == "none") NA_character_ else familiarity_natal_terms_60[[covariate]]
  has_natal <- !is.na(natal_term) && natal_term %in% rownames(coefficient_table)
  tibble::tibble(
    model = model_name,
    covariate = covariate,
    AIC = stats::AIC(fitted_model),
    natal_coefficient = if(has_natal) unname(coefficient_table[natal_term,"Estimate"]) else NA_real_,
    natal_standard_error = if(has_natal) unname(coefficient_table[natal_term,"Std. Error"]) else NA_real_,
    hfi_coefficient = unname(coefficient_table["hfi_raw_z","Estimate"]),
    random_slope_SD = extract_random_slope_sd_60(fitted_model,"hfi_raw_z"),
    singular = lme4::isSingular(fitted_model,tol = 1e-4),
    converged = extract_model_convergence_60(fitted_model))
})) %>%
  dplyr::mutate(
    delta_AIC_from_base = AIC - AIC[model == "base_model"],
    natal_CI_95 = ifelse(is.na(natal_coefficient),"--",sprintf("%.3f [%.3f; %.3f]",natal_coefficient,natal_coefficient - confidence_critical_value_60 * natal_standard_error,natal_coefficient + confidence_critical_value_60 * natal_standard_error))) %>%
  dplyr::group_by(covariate) %>%
  dplyr::mutate(delta_AIC_within_covariate = AIC - min(AIC)) %>%
  dplyr::ungroup() %>%
  dplyr::select(model,covariate,AIC,delta_AIC_from_base,delta_AIC_within_covariate,natal_coefficient,natal_CI_95,hfi_coefficient,random_slope_SD,singular,converged)

print(familiarity_model_comparison_60,n = Inf,width = Inf)

# 3.8 Export the two tables ----
utils::write.csv(two_step_results_60,file.path(results_directory_60,"GPS_persistence_two_step_natal.csv"),row.names = FALSE)
utils::write.csv(familiarity_model_comparison_60,file.path(results_directory_60,"GPS_persistence_familiarity_models.csv"),row.names = FALSE)


# 3.9 Natal effect as an effect size ----
#' Change in the Q05-Q95 HFI contrast (log-odds of remaining aerial) between an eagle from a
#' low-built natal territory (q10 of the natal descriptor) and one from a high-built territory (q90).
natal_z_contrast_60 <- stats::quantile(two_step_data_60$natal_value_z,c(0.10,0.90),names = FALSE)
hfi_difference_60 <- hfi_contrast_differences_60[["global"]]

natal_effect_size_60 <- tibble::tibble(
  natal_territory = c("low built (q10)","high built (q90)"),
  natal_q90_NT = stats::quantile(two_step_data_60$natal_value,c(0.10,0.90),names = FALSE),
  hfi_slope = unname(model_B_60$beta["(Intercept)"] + model_B_60$beta["natal_value_z"] * natal_z_contrast_60),
  log_odds_Q95_Q05_hfi = hfi_slope * hfi_difference_60)

natal_effect_difference_60 <- tibble::tibble(
  difference_log_odds = unname(model_B_60$beta["natal_value_z"]) * diff(natal_z_contrast_60) * hfi_difference_60,
  difference_se = unname(model_B_60$se["natal_value_z"]) * diff(natal_z_contrast_60) * hfi_difference_60,
  CI_95 = sprintf("[%.3f; %.3f]",difference_log_odds - confidence_critical_value_60 * difference_se,difference_log_odds + confidence_critical_value_60 * difference_se),
  relative_to_mean_contrast = difference_log_odds / (unname(model_A_60$beta["(Intercept)"]) * hfi_difference_60))

print(natal_effect_size_60); print(natal_effect_difference_60)
