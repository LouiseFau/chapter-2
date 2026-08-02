#'-------------------------------------------------------------------------------
#' Title: iSSF covariates selection ----
#' Authors : Louise Faure
#' Date : 31.07.26
#' **Info:** this script follow the Data_processing_annotation.R script where I 
#' generated random step based on a gamma and uniform distribution for step lenght
#' and turning angle, and extracted the values below each data point.  
#' **Purpose:** 
#' (1) prepare the dataset and fit several models that group covariates based on 
#' biological interpretation 
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

# Golden eagle dataset ----
annotated_data <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/issf_generated_observed_location_annotated(2).rds")

#------------------------------------------------------------------------------- STEP 1: dataset preparation and model fitting ----
#' **Steps:**
#' (i) retain complete choice sets and individuals with at least 30 landings;
#' (ii) standardize environmental covariates and prepare movement terms;
#' (iii) fit candidate biological iSSF models on the same analytical dataset;
#' (iv) compile model-fit controls and HFI coefficient summaries.

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

# 1.2 Standardize environmental covariates and prepare movement terms ----
standardization_parameters_60 <- tibble::tibble(
  variable = environmental_covariates_60,
  center = vapply(environmental_covariates_60,\(x) mean(data_model_60[[x]]),numeric(1)),
  scale = vapply(environmental_covariates_60,\(x) stats::sd(data_model_60[[x]]),numeric(1))
)

data_model_60 <- data_model_60 %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(environmental_covariates_60),~ as.numeric(scale(.x)),.names = "{.col}_z"),
                log_step_length_km = log(step_length_km),
                cos_turning_angle = cos(turning_angle_rad))

# 1.3 Define candidate biological models ----
formula_null_60 <- used ~ -1 + settlement_density_z + step_length_km + log_step_length_km + cos_turning_angle + (1 | stratum_ID) + (0 + settlement_density_z | animal_ID)

issf_formulas_60 <- list(
  null = formula_null_60,
  ia_elevation = update(formula_null_60,. ~ . + elevation_100m_z),
  ib_ruggedness = update(formula_null_60,. ~ . + ruggedness_100m_z),
  ic_elevation_ruggedness = update(formula_null_60,. ~ . + elevation_100m_z + ruggedness_100m_z),
  id_elevation_forest = update(formula_null_60,. ~ . + elevation_100m_z + prop_forest_5cells_z),
  ie_elevation_open_habitat = update(formula_null_60,. ~ . + elevation_100m_z + prop_low_vegetation_5cells_z),
  if_elevation_open_habitat_ruggedness = update(formula_null_60,. ~ . + elevation_100m_z + prop_low_vegetation_5cells_z + ruggedness_100m_z),
  ig_elevation_ridgeline = update(formula_null_60,. ~ . + elevation_100m_z + distance_to_ridgeline_100m_z)
)

# 1.4 Fit all candidate models ----
fit_issf_model_60 <- function(model_formula,data){
  glmmTMB::glmmTMB(formula = model_formula,family = poisson(link = "log"),data = data,map = list(theta = factor(c(NA,1L))),start = list(theta = c(log(fixed_stratum_sd_60),0)),control = model_control_60)
}

issf_models_60 <- purrr::map(issf_formulas_60,fit_issf_model_60,data = data_model_60)


#------------------------------------------------------------------------------- STEP 2: extract and compare candidate-model information ----
#' **Steps:**
#' (i) define the q05-q95 HFI contrast from available destinations;
#' (ii) extract model formula, AIC and population HFI estimates;
#' (iii) calculate population-level relative selection strength and its CI;
#' (iv) display and save the comparison table, highlighting the best-AIC model.

# Parameters ----
rss_confidence_level_60 <- 0.95
issf_comparison_pdf_60 <- "issf_model_comparison_60.pdf"
issf_comparison_csv_60 <- "issf_model_comparison_60.csv"
best_model_fill_60 <- "#D9EAD3"

# 2.1 Define the q05-q95 HFI contrast from available destinations ----
hfi_standardization_60 <- standardization_parameters_60 %>% dplyr::filter(variable == "settlement_density")

controle_hfi_quantiles_60 <- data_model_60 %>%
  dplyr::filter(used == 0L) %>%
  dplyr::summarise(q05_hfi = as.numeric(stats::quantile(settlement_density,0.05,na.rm = TRUE)),q95_hfi = as.numeric(stats::quantile(settlement_density,0.95,na.rm = TRUE))) %>%
  dplyr::mutate(q05_hfi_z = (q05_hfi - hfi_standardization_60$center) / hfi_standardization_60$scale,q95_hfi_z = (q95_hfi - hfi_standardization_60$center) / hfi_standardization_60$scale,hfi_q95_q05_difference_z = q95_hfi_z - q05_hfi_z)

# 2.2 Extract model information and calculate HFI RSS ----
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

print(issf_model_comparison_60,n = Inf)


# 2.3 Create and export the model-comparison table ----
issf_model_comparison_display_60 <- issf_model_comparison_60 %>%
  dplyr::select(model,AIC,delta_AIC,hfi_coefficient,hfi_standard_error,RSS_q95_vs_q05,RSS_CI_95,best_AIC)

issf_model_comparison_table_60 <- issf_model_comparison_display_60 %>%
  gt::gt(rowname_col = "model") %>%
  gt::tab_header(title = gt::md("**Landing iSSF candidate-model comparison**"),subtitle = paste0("Population HFI RSS compares available-landscape q95 with q05: ",round(controle_hfi_quantiles_60$q05_hfi,3)," to ",round(controle_hfi_quantiles_60$q95_hfi,3))) %>%
  gt::cols_label(AIC = "AIC",delta_AIC = "ΔAIC",hfi_coefficient = "HFI coefficient",hfi_standard_error = "HFI SE",RSS_q95_vs_q05 = "RSS q95:q05",RSS_CI_95 = "RSS [95% CI]") %>%
  gt::fmt_number(columns = c(AIC,delta_AIC),decimals = 1) %>%
  gt::fmt_number(columns = c(hfi_coefficient,hfi_standard_error,RSS_q95_vs_q05),decimals = 3) %>%
  gt::cols_hide(columns = best_AIC) %>%
  gt::tab_style(style = list(gt::cell_fill(color = best_model_fill_60),gt::cell_text(weight = "bold")),locations = gt::cells_body(rows = best_AIC)) %>%
  gt::cols_width(AIC ~ gt::px(80),delta_AIC ~ gt::px(70),hfi_coefficient ~ gt::px(105),hfi_standard_error ~ gt::px(85),RSS_q95_vs_q05 ~ gt::px(100),RSS_CI_95 ~ gt::px(145)) %>%
  gt::tab_options(table.font.size = gt::px(9),data_row.padding = gt::px(4),heading.title.font.size = gt::px(15),heading.subtitle.font.size = gt::px(10),table.width = gt::pct(100))

issf_model_comparison_table_60
gt::gtsave(issf_model_comparison_table_60,filename = issf_comparison_pdf_60,vwidth = 1800,vheight = 1100,expand = 10)
readr::write_csv(issf_model_comparison_60 %>% dplyr::select(-best_AIC),issf_comparison_csv_60)
