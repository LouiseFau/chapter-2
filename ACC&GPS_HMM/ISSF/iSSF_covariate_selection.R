#'-------------------------------------------------------------------------------
#' Title: iSSF covariates selection ----
#' Authors : Louise Faure
#' Date : 31.07.26
#' **Info:** this script follow the Data_processing_annotation.R script where I 
#' generated random step based on a gamma and uniform distribution for step lenght
#' and turning angle, and extracted the environmental values below each data point.  
#' **Purpose:** 
#' (1) prepare the dataset and determine which familiarity index best explained
#' the variance within the dataset. 
#' (2) fit several models that group covariates based on biological interpretation 
#' (3) select the model with lowest AIC, few covariates, RSS for q05 and q95, CI 
#' of the RSS, standard error for HFI
#' (4) for the selected model control the VIF and Pearson correlation coefficient
#' ------------------------------------------------------------------------------

# Libraries ----
library(tidyverse)
library(glmmTMB)
library(corrr)
library(gt)

# Parameters ----
minimum_landings_60 <- 30L
fixed_stratum_sd_60 <- 1e3
model_control_60 <- glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000))

# Colors ----
clr <- oce::oceColorsPalette(100)[9] #was [2] before
clr_light <- oce::oceColorsPalette(100)[10]
clr2 <- oce::oceColorsPalette(100)[80]

# Golden eagle dataset ----
annotated_data <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/issf_generated_observed_location_annotated(2).rds")

#------------------------------------------------------------------------------- STEP 1: dataset preparation and model fitting ----
#' **Steps:**
#' (i) retain complete choice sets and individuals with at least 30 landings;
#' (ii) standardize environmental covariates and prepare movement terms;

# 1.1 Prepare complete choice sets and retain supported individuals ----
environmental_covariates_60 <- c("settlement_density","elevation_100m","ruggedness_100m","prop_forest_5cells","prop_low_vegetation_5cells","distance_to_ridgeline_100m")
required_columns_60 <- c("used","stratum","individual.local.identifier","step_length_km","turning_angle_rad",environmental_covariates_60)
controle_missing_columns_60 <- setdiff(required_columns_60,names(annotated_data))
if(length(controle_missing_columns_60) > 0L) stop("Missing columns: ",paste(controle_missing_columns_60,collapse = ", "))

data_complete_60 <- annotated_data %>%
  dplyr::mutate(used = as.integer(used),individual.local.identifier = as.character(individual.local.identifier),stratum = as.character(stratum)) %>%
  dplyr::filter(!is.na(individual.local.identifier),!is.na(stratum),used %in% c(0L,1L),step_length_km > 0,dplyr::if_all(dplyr::all_of(c("step_length_km","turning_angle_rad",environmental_covariates_60)),~ !is.na(.x) & is.finite(.x))) %>%
  dplyr::group_by(individual.local.identifier,stratum) %>%
  dplyr::filter(sum(used == 1L) == 1L,sum(used == 0L) >= 1L) %>%
  dplyr::ungroup()

controle_individual_support_60 <- data_complete_60 %>%
  dplyr::filter(used == 1L) %>%
  dplyr::count(individual.local.identifier,name = "n_landings") %>%
  dplyr::mutate(retained = n_landings >= minimum_landings_60) %>%
  dplyr::arrange(n_landings)

retained_individuals_60 <- controle_individual_support_60 %>%
  dplyr::filter(retained) %>%
  dplyr::pull(individual.local.identifier)

data_model_60 <- data_complete_60 %>%
  dplyr::filter(individual.local.identifier %in% retained_individuals_60) %>%
  dplyr::mutate(animal_ID = factor(individual.local.identifier),stratum_ID = interaction(individual.local.identifier,stratum,drop = TRUE,lex.order = TRUE))

# 1.2 Standardize environmental and familiarity covariates and prepare movement terms ----
familiarity_covariates_60 <- c("position_NT","above_natal_NT","excess_NT")
controle_missing_familiarity_60 <- setdiff(familiarity_covariates_60,names(data_model_60))
if(length(controle_missing_familiarity_60) > 0L) stop("Missing columns: ",paste(controle_missing_familiarity_60,collapse = ", "))

standardized_covariates_60 <- c(environmental_covariates_60,familiarity_covariates_60)

standardization_parameters_60 <- tibble::tibble(
  variable = standardized_covariates_60,
  center = vapply(standardized_covariates_60,\(x) mean(data_model_60[[x]],na.rm = TRUE),numeric(1)),
  scale = vapply(standardized_covariates_60,\(x) stats::sd(data_model_60[[x]],na.rm = TRUE),numeric(1))
)

data_model_60 <- data_model_60 %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(standardized_covariates_60),~ as.numeric((.x - mean(.x,na.rm = TRUE)) / stats::sd(.x,na.rm = TRUE)),.names = "{.col}_z"),
                log_step_length_km = log(step_length_km),
                cos_turning_angle = cos(turning_angle_rad))


#------------------------------------------------------------------------------- STEP 2 : fit several biological models to choose environmental covariates ----
#' **Steps:**
#' (i) fit candidate biological iSSF models on the same analytical dataset;
#' (ii) define the q05-q95 HFI contrast from available destinations;
#' (iii) extract model formula, AIC and population HFI estimates;
#' (iv) calculate population-level relative selection strength and its CI;
#' (v) display and save the comparison table, highlighting the best-AIC model.
#' (vi) validate with RMSE and export best model

# Parameters ----
rss_confidence_level_60 <- 0.95
issf_comparison_pdf_60 <- "issf_model_comparison_60.pdf"
issf_comparison_csv_60 <- "issf_model_comparison_60.csv"
best_model_fill_60 <- "#D9EAD3"

# 2.1 Define candidate biological models ----
formula_null_60 <- used ~ -1 + settlement_density_z + step_length_km + log_step_length_km + cos_turning_angle + (1 | stratum_ID) + (0 + settlement_density_z | animal_ID)

issf_formulas_60 <- list(
  null = formula_null_60,
  ia_elevation = update(formula_null_60,. ~ . + elevation_100m_z),
  ib_ruggedness = update(formula_null_60,. ~ . + ruggedness_100m_z),
  ic_elevation_ruggedness = update(formula_null_60,. ~ . + elevation_100m_z + ruggedness_100m_z),
  id_elevation_forest = update(formula_null_60,. ~ . + elevation_100m_z + prop_forest_5cells_z),
  ie_elevation_open_habitat = update(formula_null_60,. ~ . + elevation_100m_z + prop_low_vegetation_5cells_z),
  if_elevation_open_habitat_ruggedness = update(formula_null_60,. ~ . + elevation_100m_z + prop_low_vegetation_5cells_z + ruggedness_100m_z),
  ig_elevation_ridgeline = update(formula_null_60,. ~ . + elevation_100m_z + distance_to_ridgeline_100m_z))

fit_issf_model_60 <- function(model_formula,data){
  glmmTMB::glmmTMB(formula = model_formula,family = poisson(link = "log"),data = data,map = list(theta = factor(c(NA,1L))),start = list(theta = c(log(fixed_stratum_sd_60),0)),control = model_control_60)}

issf_models_60 <- purrr::map(issf_formulas_60,fit_issf_model_60,data = data_model_60)

# 2.2 Define the q05-q95 HFI contrast from available destinations ----
hfi_standardization_60 <- standardization_parameters_60 %>% dplyr::filter(variable == "settlement_density")

controle_hfi_quantiles_60 <- data_model_60 %>%
  dplyr::filter(used == 0L) %>%
  dplyr::summarise(q05_hfi = as.numeric(stats::quantile(settlement_density,0.05,na.rm = TRUE)),q95_hfi = as.numeric(stats::quantile(settlement_density,0.95,na.rm = TRUE))) %>%
  dplyr::mutate(q05_hfi_z = (q05_hfi - hfi_standardization_60$center) / hfi_standardization_60$scale,q95_hfi_z = (q95_hfi - hfi_standardization_60$center) / hfi_standardization_60$scale,hfi_q95_q05_difference_z = q95_hfi_z - q05_hfi_z)

# 2.3 Extract model information and calculate HFI RSS ----
extract_issf_information_60 <- function(model_object,model_name){
  coefficient_table <- summary(model_object)$coefficients$cond
  if(!"settlement_density_z" %in% rownames(coefficient_table)) stop("HFI coefficient not found in model: ",model_name)
  hfi_coefficient <- unname(coefficient_table["settlement_density_z","Estimate"])
  hfi_standard_error <- unname(coefficient_table["settlement_density_z","Std. Error"])
  hfi_difference_z <- controle_hfi_quantiles_60$hfi_q95_q05_difference_z
  confidence_multiplier <- stats::qnorm(1 - (1 - rss_confidence_level_60) / 2)
  log_rss <- hfi_coefficient * hfi_difference_z
  log_rss_standard_error <- hfi_standard_error * abs(hfi_difference_z)
  tibble::tibble(model = model_name,AIC = stats::AIC(model_object),hfi_coefficient = hfi_coefficient,hfi_standard_error = hfi_standard_error,RSS_q95_vs_q05 = exp(log_rss),RSS_confidence_low = exp(log_rss - confidence_multiplier * log_rss_standard_error),RSS_confidence_high = exp(log_rss + confidence_multiplier * log_rss_standard_error))
}

issf_model_comparison_60 <- purrr::imap_dfr(issf_models_60,~ extract_issf_information_60(model_object = .x,model_name = .y)) %>%
  dplyr::mutate(delta_AIC = AIC - min(AIC),best_AIC = AIC == min(AIC),RSS_CI_95 = sprintf("%.3f [%.3f; %.3f]",RSS_q95_vs_q05,RSS_confidence_low,RSS_confidence_high)) %>%
  dplyr::arrange(AIC) %>%
  dplyr::select(model,AIC,delta_AIC,hfi_coefficient,hfi_standard_error,RSS_q95_vs_q05,RSS_confidence_low,RSS_confidence_high,RSS_CI_95,best_AIC)

print(issf_model_comparison_60,n = Inf) # best model is elevation_open_habitat_ruggedness, hfi coeff = -2.12


# 2.4 remove one by one the variables from the best model and check their influence on the hfi_coefficient  ----
best_model_name_60 <- issf_model_comparison_60 %>%
  dplyr::slice_min(AIC,n = 1,with_ties = FALSE) %>%
  dplyr::pull(model)

best_model_60 <- issf_models_60[[best_model_name_60]]
best_formula_60 <- stats::formula(best_model_60)

environmental_terms_to_remove_60 <- c(elevation = "elevation_100m_z",ruggedness = "ruggedness_100m_z",open_habitat = "prop_low_vegetation_5cells_z")

best_model_terms_60 <- attr(stats::terms(best_formula_60),"term.labels")
terms_to_remove_60 <- environmental_terms_to_remove_60[environmental_terms_to_remove_60 %in% best_model_terms_60]

reduced_formulas_60 <- purrr::imap(
  terms_to_remove_60,
  ~ stats::update.formula(best_formula_60,stats::as.formula(paste(". ~ . -",.x))))
names(reduced_formulas_60) <- paste0("without_",names(terms_to_remove_60))

issf_hfi_sensitivity_models_60 <- c(
  list(full_model = best_model_60),
  purrr::map(reduced_formulas_60,fit_issf_model_60,data = data_model_60))

sensitivity_metadata_60 <- tibble::tibble(
  model = names(issf_hfi_sensitivity_models_60),
  removed_covariate = c("none",names(terms_to_remove_60)))

issf_hfi_sensitivity_60 <- purrr::imap_dfr(issf_hfi_sensitivity_models_60,function(model_object,model_name) {
  coefficient_table <- summary(model_object)$coefficients$cond
  tibble::tibble(
    model = model_name,
    AIC = stats::AIC(model_object),
    hfi_coefficient = unname(coefficient_table["settlement_density_z","Estimate"]),
    hfi_standard_error = unname(coefficient_table["settlement_density_z","Std. Error"]))
}) %>% dplyr::left_join(sensitivity_metadata_60,by = "model")

full_AIC_60 <- issf_hfi_sensitivity_60$AIC[issf_hfi_sensitivity_60$model == "full_model"]
full_hfi_coefficient_60 <- issf_hfi_sensitivity_60$hfi_coefficient[issf_hfi_sensitivity_60$model == "full_model"]

issf_hfi_sensitivity_60 <- issf_hfi_sensitivity_60 %>%
  dplyr::mutate(
    delta_AIC_from_full = AIC - full_AIC_60,
    hfi_coefficient_difference = hfi_coefficient - full_hfi_coefficient_60,
    hfi_magnitude_change_percent = 100 * (abs(hfi_coefficient) - abs(full_hfi_coefficient_60)) / abs(full_hfi_coefficient_60),
    hfi_CI_low = hfi_coefficient - 1.96 * hfi_standard_error,
    hfi_CI_high = hfi_coefficient + 1.96 * hfi_standard_error) %>%
  dplyr::arrange(match(model,c("full_model",paste0("without_",names(terms_to_remove_60))))) %>%
  dplyr::select(model,removed_covariate,AIC,delta_AIC_from_full,hfi_coefficient,hfi_standard_error,
                hfi_coefficient_difference,hfi_magnitude_change_percent,hfi_CI_low,hfi_CI_high)

print(issf_hfi_sensitivity_60,n = Inf)

# 2.5 control pearson correlation and vif in the best model ----
data_model_60 %>%
  dplyr::select(c("settlement_density","elevation_100m","ruggedness_100m",
                  "prop_low_vegetation_5cells","step_length_km")) %>%
  corrr::correlate()

performance::check_collinearity(best_model_60)

# 2.6 Print model statistical summary ----
summary(best_model_60)

# 2.7.1 Extract coefficient estimates and confidence intervals ----
confint(best_model_60)

# 2.7.2 Extract individual-specific random effects ----
ranef(best_model_60)[[1]]$animal_ID

# 3.2 Model validation obtained by calculating the RMSE ----
performance::performance_rmse(best_model_60) # 0.13

# Export best model 
saveRDS(list(model = best_model_60, model_name = best_model_name_60, formula = best_formula_60, standardization = standardization_parameters_60, hfi_quantiles = controle_hfi_quantiles_60), file = "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/ISSF/issf_best_model_60.rds")

# VISUALISATION n°1 : PLOT coefficient estimates -------------------------------------------------------------------------------------------------------- 
graph_60 <- confint(best_model_60) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Factor") %>%
  dplyr::filter(!grepl("^Std.Dev",Factor))

colnames(graph_60)[c(2,3)] <- c("Lower","Upper")

coefficient_labels_60 <- c(
  settlement_density_z = "HFI",
  step_length_km = "Step length",
  log_step_length_km = "Log step length",
  cos_turning_angle = "Cosine turning angle",
  elevation_100m_z = "Elevation",
  prop_low_vegetation_5cells_z = "Low vegetation",
  ruggedness_100m_z = "Ruggedness"
)

graph_60 <- graph_60 %>%
  dplyr::mutate(Factor = factor(Factor,levels = rev(names(coefficient_labels_60))))

coefs_60 <- ggplot2::ggplot(
  graph_60,
  ggplot2::aes(x = Estimate,y = Factor)
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray",
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    color = clr,
    size = 2.4
  ) +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = Lower,xmax = Upper),
    color = clr,
    linewidth = 0.8
  ) +
  ggplot2::scale_y_discrete(
    labels = coefficient_labels_60
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    text = ggplot2::element_text(
      family = "Baskerville",
      size = 14
    ),
    axis.text.x = ggplot2::element_text(
      family = "Baskerville",
      size = 14
    ),
    axis.text.y = ggplot2::element_text(
      family = "Baskerville",
      size = 15
    )
  )

print(coefs_60)

ggsave(coefs_60, filename = "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/ISSF/coeffs.svg", 
       width = 7, height = 2, dpi = 400)




# VISUALISATION n°2 : individual-specific coefficients ------------------------------------------------------------------ 
population_hfi_60 <- unname(glmmTMB::fixef(best_model_60)$cond["settlement_density_z"])

individual_hfi_coefficient_60 <- glmmTMB::ranef(best_model_60,condVar = TRUE) %>%
  as.data.frame() %>%
  dplyr::filter(grpvar == "animal_ID",term == "settlement_density_z") %>%
  dplyr::mutate(individual_id = grp,
                individual_coefficient = population_hfi_60 + condval,
                CI_low = population_hfi_60 + condval - 1.96 * condsd,
                CI_high = population_hfi_60 + condval + 1.96 * condsd) %>%
  dplyr::arrange(individual_coefficient) %>%
  dplyr::mutate(individual_id = factor(individual_id,levels = individual_id))

individual_hfi_coefficient_plot_60 <- ggplot2::ggplot(
  individual_hfi_coefficient_60,
  ggplot2::aes(x = individual_coefficient,y = individual_id)
) +
  ggplot2::geom_vline(xintercept = 0,color = "black",linewidth = 0.5) +
  ggplot2::geom_vline(xintercept = population_hfi_60,linetype = "dashed",color = "gray",linewidth = 0.6) +
  ggplot2::geom_linerange(ggplot2::aes(xmin = CI_low,xmax = CI_high),linewidth = 0.6) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::labs(x = "Individual HFI coefficient",y = "",
                subtitle = "Dashed line: population HFI coefficient; solid line: no HFI effect") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 8))

print(individual_hfi_coefficient_plot_60)

ggsave(plot = individual_hfi_coefficient_plot_60, filename = "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/ISSF/ind_coefs.pdf", 
       width = 7, height = 9, dpi = 300)

# VISUALISATION n°3: individual-specific HFI RSS ---------------------------------------------------------------------
hfi_difference_60 <- controle_hfi_quantiles_60$hfi_q95_q05_difference_z

individual_hfi_rss_60 <- individual_hfi_coefficient_60 %>%
  dplyr::mutate(
    RSS_Q95_vs_Q05 = exp(individual_coefficient * hfi_difference_60),
    RSS_CI_low = exp(CI_low * hfi_difference_60),
    RSS_CI_high = exp(CI_high * hfi_difference_60)
  )

population_hfi_rss_60 <- exp(population_hfi_60 * hfi_difference_60)

individual_hfi_rss_plot_60 <- ggplot2::ggplot(
  individual_hfi_rss_60,
  ggplot2::aes(x = RSS_Q95_vs_Q05,y = individual_id)
) +
  ggplot2::geom_vline(xintercept = 1,color = "black",linewidth = 0.5) +
  ggplot2::geom_vline(xintercept = population_hfi_rss_60,
                      linetype = "dashed",color = "gray",linewidth = 0.6) +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = RSS_CI_low,xmax = RSS_CI_high),
    linewidth = 0.6
  ) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    x = "Individual RSS: HFI Q95 versus Q05",
    y = "",
    subtitle = "Dashed line: population RSS; solid line: RSS = 1"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 8))

print(individual_hfi_rss_plot_60)
ggsave(plot = individual_hfi_rss_plot_60, filename = "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/ISSF/ind_rss.pdf", 
       width = 7, height = 9, dpi = 300)


# 3.3 Control coefficient stability by removing individuals ----
#' **Hypothesis**: individuals with more landing events have a different behavior
#' from the other individuals. To test this assumption, we can remove individual 
#' with more landing event. 
