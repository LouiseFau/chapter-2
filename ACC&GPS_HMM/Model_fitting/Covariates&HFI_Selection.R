#' -----------------------------------------------------------------------------
# Title: Covariate and HFI variables selection ----
#' Authors : Louise Faure
#' Date : 20.07.26
#' 
#' **Info:** this script follow the Transition_Matrix_Preparation.R script which 
#' standardize covariates and prepare the transition matrix using aerial state 
#' as the reference. This script can be applied to both GPS and ACC data. 
#' 
#' **Main steps:**
#' (1) define the behavioural backbone model using gps dataset:
#'     (ii) retain the diel cosinor and individual random intercept;
#'     (iii) compare biologically plausible linear, quadratic or smooth forms for
#'          age and elapsed aerial duration for the acc and gps models 
#'     (iv) retain the most parsimonious supported structure.
#'
#' (2) identify potential confounders and HFI metrics for gps and acc datasets
#'     (i) inspect correlations using Pearson correlation coefficient and remove
#'     highly correlated variables
#'     (ii) identify HFI metrics supports to the data 
#'
#' (4) select model covariates
#'     (i) define a group of biologically meaningful groups 
#'     (ii) fit models for both acc and gps data and calculate the logg odds 
#'     between q05 and q95.
#'     (iii) compile in two tables (one for gps and the other for acc): the 
#'     degree of freedom, the AIC, the delta AIC, hfi coefficient, q95-Q05 log-odds
#'     (iv) compare aic for population density and settlement density
#'     (iv) check fof VIF in the selected model
#'
#' (5) print summary results for both dataset (coefficient HFI, log odds, and 
#' log odd CI) 
#' -----------------------------------------------------------------------------


# library
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(ggplot2)
library(lme4)
library(gt)         # for ploting results
library(htmltools)  # for ploting results


# acc and gps data weighted and standardized data ----
GE_acc_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/acc_weighted(2).rds")
GE_gps_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/gps_weighted(2).rds")

#------------------------------------------------------------------------------- STEP 1: define behavioural backbone models ----
#' **Steps:**
#' (i) prepare three modelling datasets:
#'     a) GPS aerial versus terrestrial;
#'     b) ACC aerial versus terrestrial;
#'     c) ACC feeding versus resting conditional on a terrestrial transition;
#' (ii) retain the diel cosinor and individual random intercept;
#' (iii) compare linear or quadratic forms for age and elapsed aerial duration;
#' (iv) retain the simplest supported structure when delta AIC <= 4.

# Parameters ----
delta_aic_threshold_60 <- 4
duration_terms_60 <- list(linear = "duration_z",quadratic = c("duration_z","duration_z2"))
control_glmer_60 <- lme4::glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))

# 1.1 Prepare modelling datasets ----
add_process_individual_weights_60 <- function(data) {
  data %>%
    dplyr::select(-dplyr::any_of(c("backbone_n_transitions","backbone_weight_raw","backbone_weight"))) %>%
    dplyr::group_by(individual_id) %>%
    dplyr::mutate(
      backbone_n_transitions = dplyr::n(),
      backbone_weight_raw = 1 / backbone_n_transitions) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      backbone_weight = backbone_weight_raw / mean(backbone_weight_raw))}

prepare_backbone_dataset_60 <- function(data,response_name) {required_variables <- c("individual_id",response_name,"cos_diel_c","sin_diel_c",
    "duration_z")
  data %>%
    dplyr::mutate(
      individual_id = factor(individual_id),
      duration_z2 = duration_z^2) %>%
    tidyr::drop_na(dplyr::all_of(required_variables)) %>%
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(c(response_name,"cos_diel_c","sin_diel_c","duration_z")),
        is.finite)) %>%
    dplyr::mutate(individual_id = droplevels(individual_id)) %>%
    add_process_individual_weights_60()}

gps_backbone_data_60 <- GE_gps_hfi %>% dplyr::mutate(remain_aerial = as.integer(remain_aerial)) %>%
  prepare_backbone_dataset_60("remain_aerial")

acc_aerial_backbone_data_60 <- GE_acc_hfi %>% dplyr::mutate(remain_aerial = as.integer(transition_destination == "aerial")) %>%
  prepare_backbone_dataset_60("remain_aerial")

acc_feeding_resting_backbone_data_60 <- GE_acc_hfi %>% dplyr::filter(transition_destination %in% c("resting","feeding")) %>%
  dplyr::mutate(feeding_vs_resting = as.integer(transition_destination == "feeding")) %>%
  prepare_backbone_dataset_60("feeding_vs_resting")

backbone_process_data_60 <- list(
  GPS_aerial_vs_terrestrial = gps_backbone_data_60,
  ACC_aerial_vs_terrestrial = acc_aerial_backbone_data_60,
  ACC_feeding_vs_resting = acc_feeding_resting_backbone_data_60)

backbone_process_response_60 <- c(GPS_aerial_vs_terrestrial = "remain_aerial",ACC_aerial_vs_terrestrial = "remain_aerial",ACC_feeding_vs_resting = "feeding_vs_resting")

# 1.2 Define candidate models ----
backbone_model_metadata_60 <- tidyr::expand_grid(
  process = names(backbone_process_data_60),
  duration_form = names(duration_terms_60)
) %>%
  dplyr::mutate(
    response = unname(backbone_process_response_60[process]),
    duration_complexity = dplyr::recode(
      duration_form,
      linear = 1L,
      quadratic = 2L
    ),
    model_id = paste(process,duration_form,sep = "__")
  )

candidate_backbone_formulas_60 <- stats::setNames(
  lapply(seq_len(nrow(backbone_model_metadata_60)),function(i) {
    model_terms <- c(
      "cos_diel_c",
      "sin_diel_c",
      duration_terms_60[[backbone_model_metadata_60$duration_form[[i]]]],
      "(1 | individual_id)"
    )
    
    stats::reformulate(
      termlabels = model_terms,
      response = backbone_model_metadata_60$response[[i]]
    )
  }),
  backbone_model_metadata_60$model_id
)

# 1.3 Fit linear and quadratic models ----
fit_backbone_model_60 <- function(model_id) {
  model_information <- backbone_model_metadata_60 %>%
    dplyr::filter(.data$model_id == .env$model_id)
  
  lme4::glmer(
    formula = candidate_backbone_formulas_60[[model_id]],
    data = backbone_process_data_60[[model_information$process[[1]]]],
    weights = backbone_weight,
    family = stats::binomial(link = "logit"),
    nAGQ = 1,
    control = control_glmer_60
  )
}

backbone_models_60 <- stats::setNames(lapply(backbone_model_metadata_60$model_id,fit_backbone_model_60),backbone_model_metadata_60$model_id)

# 1.4 Compare duration structures ----
backbone_model_comparison_60 <- dplyr::bind_rows(
  lapply(names(backbone_models_60),function(model_id) {
    fitted_model <- backbone_models_60[[model_id]]
    
    model_information <- backbone_model_metadata_60 %>%
      dplyr::filter(.data$model_id == .env$model_id)
    
    model_loglik <- stats::logLik(fitted_model)
    convergence_message <- fitted_model@optinfo$conv$lme4$messages
    optimizer_code <- fitted_model@optinfo$conv$opt
    
    tibble::tibble(
      process = model_information$process[[1]],
      model_id = model_id,
      duration_form = model_information$duration_form[[1]],
      duration_complexity = model_information$duration_complexity[[1]],
      n_observations = stats::nobs(fitted_model),
      model_df = attr(model_loglik,"df"),
      log_likelihood = as.numeric(model_loglik),
      AIC = stats::AIC(fitted_model),
      converged = is.null(convergence_message) &&
        (is.null(optimizer_code) || all(optimizer_code == 0)),
      singular = lme4::isSingular(fitted_model,tol = 1e-4),
      convergence_message = if(is.null(convergence_message)) {
        NA_character_
      } else {
        paste(convergence_message,collapse = "; ")
      }
    )
  })
) %>%
  dplyr::group_by(process) %>%
  dplyr::arrange(AIC,.by_group = TRUE) %>%
  dplyr::mutate(
    delta_AIC = AIC - min(AIC),
    akaike_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)),
    supported = delta_AIC <= delta_aic_threshold_60
  ) %>%
  dplyr::ungroup()

# 1.5 Select the simplest supported duration structure ----
selected_backbone_information_60 <- backbone_model_comparison_60 %>%
  dplyr::filter(converged,supported) %>%
  dplyr::group_by(process) %>%
  dplyr::arrange(
    singular,
    duration_complexity,
    AIC,
    .by_group = TRUE
  ) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

backbone_model_comparison_60 <- backbone_model_comparison_60 %>%
  dplyr::left_join(
    selected_backbone_information_60 %>%
      dplyr::transmute(process,model_id,selected = TRUE),
    by = c("process","model_id")
  ) %>%
  dplyr::mutate(selected = dplyr::coalesce(selected,FALSE))

# 1.6 Store selected models and formulas ----
selected_backbone_models_60 <- stats::setNames(
  lapply(names(backbone_process_data_60),function(process_name) {
    selected_model_id <- selected_backbone_information_60$model_id[
      selected_backbone_information_60$process == process_name
    ]
    
    backbone_models_60[[selected_model_id]]
  }),
  names(backbone_process_data_60))

# 1.7 Selection of duration form ----
print(selected_backbone_information_60 %>% dplyr::select(process,duration_form,AIC,delta_AIC,akaike_weight,singular),n = Inf,width = Inf)

# process                   duration_form    AIC delta_AIC akaike_weight singular
# ACC_aerial_vs_terrestrial linear        12927.      2.75         0.202 FALSE   
# ACC_feeding_vs_resting    linear         3683.      0            0.665 FALSE   
# GPS_aerial_vs_terrestrial quadratic     14046.      0            0.978 FALSE 

selected_backbone_formulas_60 <- list(
  GPS_aerial_vs_terrestrial =
    remain_aerial ~ cos_diel_c + sin_diel_c + duration_z + (1 | individual_id),
  ACC_aerial_vs_terrestrial =
    remain_aerial ~ cos_diel_c + sin_diel_c + duration_z + (1 | individual_id),
  ACC_feeding_vs_resting =
    feeding_vs_resting ~ cos_diel_c + sin_diel_c + duration_z + (1 | individual_id))



#------------------------------------------------------------------------------- STEP 2: identify potential confounders and HFI metrics ----
#' **Steps:**
#' (i) inspect Pearson correlations separately for the three processes;
#' (ii) remove redundant environmental covariates using biological justification;
#' (iii) evaluate absolute and relative HFI support in the three processes;
#' (iv) retain HFI metrics supported across GPS, ACC aerial-versus-terrestrial
#'      and ACC feeding-versus-resting.

# Parameters ----
correlation_threshold_60 <- 0.70
environmental_covariates_60 <- c("elevation_100m","ruggedness_100m","slope_100m",
                                 "distance_to_ridgeline_100m","prop_forest_5cells",
                                 "prop_low_vegetation_5cells","prop_rocky_terrain_5cells", 
                                 "population_density", "settlement_density")
hfi_classes_60 <- c(0,0.10,0.20,0.60,0.80,Inf)
hfi_labels_60 <- c("0–0.10","0.10–0.20","0.20–0.60","0.60–0.80",">0.80")
process_labels_60 <- c(
  GPS_aerial_vs_terrestrial = "GPS: aerial vs terrestrial",
  ACC_aerial_vs_terrestrial = "ACC: aerial vs terrestrial",
  ACC_feeding_vs_resting = "ACC: feeding vs resting")

step2_process_data_60 <- backbone_process_data_60

# 2.1 Inspect correlations among environmental covariates ----
correlation_matrix_list_60 <- stats::setNames(lapply(step2_process_data_60,function(data) {
  matrix_60 <- stats::cor(data[,environmental_covariates_60],method = "pearson",use = "pairwise.complete.obs")
  dimnames(matrix_60) <- list(environmental_covariates_60,environmental_covariates_60)
  matrix_60
}),names(step2_process_data_60))

compile_lower_correlation_60 <- function(matrix_60,process_name) {
  indices <- which(lower.tri(matrix_60),arr.ind = TRUE)
  tibble::tibble(process = process_name,variable_1 = rownames(matrix_60)[indices[,1]],
                 variable_2 = colnames(matrix_60)[indices[,2]],pearson_r = matrix_60[indices])
}

correlation_summary_60 <- dplyr::bind_rows(lapply(names(correlation_matrix_list_60),function(process_name) {
  compile_lower_correlation_60(correlation_matrix_list_60[[process_name]],process_name)
})) %>%
  dplyr::mutate(process_label = unname(process_labels_60[process])) %>%
  dplyr::arrange(process,dplyr::desc(abs(pearson_r)))

invisible(lapply(names(correlation_matrix_list_60),function(process_name) {
  corrplot::corrplot(
    correlation_matrix_list_60[[process_name]],method = "color",type = "lower",
    order = "original",col = corrplot::COL2("RdBu",200),col.lim = c(-1,1),
    diag = FALSE,addCoef.col = "black",number.digits = 2,number.cex = 0.6,
    tl.col = "black",tl.cex = 0.6,tl.srt = 45,cl.cex = 0.8,
    addgrid.col = "white",
    title = paste("Pearson correlations —",process_labels_60[[process_name]]),
    mar = c(0,0,3,0))
}))

# 2.2 Retain non-redundant environmental covariates ----
excluded_environmental_covariates_60 <- tibble::tibble(
  variable = c("slope_100m","prop_rocky_terrain_5cells"),
  reason = c("Redundant with ruggedness and less representative of topographic complexity",
             "Associated with elevation and of limited additional biological interpretability"))

selected_environmental_covariates_60 <- setdiff(environmental_covariates_60,excluded_environmental_covariates_60$variable)

# 2.3 Evaluate support along relative human-pressure gradients ----
#' Population and settlement values are independently rescaled between their
#' empirical Q05 and Q95 within each behavioural process. The interval is then
#' divided into ten equal-width relative classes.
human_pressure_variables_60 <- c("population_density","settlement_density")
human_pressure_labels_60 <- c(population_density = "Population density",settlement_density = "Settlement density")
relative_class_labels_60 <- c("0–10%","10–20%","20–30%","30–40%","40–50%","50–60%","60–70%","70–80%","80–90%","90–100%")

human_pressure_long_60 <- dplyr::bind_rows(lapply(names(step2_process_data_60),function(process_name) {
  step2_process_data_60[[process_name]] %>%
    dplyr::select(
      individual_id,
      burst_id,
      transition_destination,
      dplyr::all_of(human_pressure_variables_60)) %>%
    dplyr::mutate(
      process = process_name,
      process_label = unname(process_labels_60[process_name]),
      burst_uid = interaction(individual_id,burst_id,drop = TRUE,lex.order = TRUE)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(human_pressure_variables_60),
      names_to = "human_variable",
      values_to = "human_value"
    )})) %>%dplyr::filter(is.finite(human_value))

# 2.4 Calculate empirical Q05 and Q95 ----
human_pressure_thresholds_60 <- human_pressure_long_60 %>%
  dplyr::group_by(process,process_label,human_variable) %>%
  dplyr::summarise(
    empirical_q05 = as.numeric(stats::quantile(human_value,0.05,names = FALSE)),
    empirical_q95 = as.numeric(stats::quantile(human_value,0.95,names = FALSE)),
    .groups = "drop")

# 2.5 Divide each empirical Q05-Q95 gradient into ten relative classes ----
human_pressure_class_support_60 <- human_pressure_long_60 %>%
  dplyr::left_join(
    human_pressure_thresholds_60,
    by = c("process","process_label","human_variable")) %>%
  dplyr::filter(
    human_value >= empirical_q05,
    human_value <= empirical_q95,
    empirical_q95 > empirical_q05) %>%
  dplyr::mutate(
    relative_position = (human_value - empirical_q05) / (empirical_q95 - empirical_q05),
    relative_class_number = pmin(10L,floor(relative_position * 10) + 1L),
    relative_class = factor(
      relative_class_labels_60[relative_class_number],
      levels = relative_class_labels_60)
  ) %>%
  dplyr::group_by(
    process,
    process_label,
    human_variable,
    relative_class
  ) %>%
  dplyr::summarise(
    n_transitions = dplyr::n(),
    n_individuals = dplyr::n_distinct(individual_id),
    n_bursts = dplyr::n_distinct(burst_uid),
    .groups = "drop"
  ) %>%
  dplyr::group_by(process,process_label,human_variable) %>%
  tidyr::complete(
    relative_class = factor(
      relative_class_labels_60,
      levels = relative_class_labels_60
    ),
    fill = list(
      n_transitions = 0L,
      n_individuals = 0L,
      n_bursts = 0L
    )) %>% dplyr::ungroup()

# Prepare plot data in the same format as visualisation n°1 ----
human_pressure_class_plot_data_60 <- human_pressure_class_support_60 %>%
  dplyr::select(
    process_label,
    human_variable,
    relative_class,
    n_transitions,
    n_individuals,
    n_bursts
  ) %>%
  tidyr::pivot_longer(
    cols = c(n_transitions,n_individuals,n_bursts),
    names_to = "support_measure",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    support_measure = dplyr::recode(
      support_measure,
      n_transitions = "Transitions",
      n_individuals = "Individuals",
      n_bursts = "Bursts"
    ),
    support_measure = factor(
      support_measure,
      levels = c("Transitions","Individuals","Bursts")
    ),
    human_variable = dplyr::recode(
      human_variable,
      population_density = "Population density",
      settlement_density = "Settlement density"
    ),
    human_variable = factor(
      human_variable,
      levels = c("Settlement density","Population density")
    ),
    relative_class = factor(
      relative_class,
      levels = relative_class_labels_60
    )
  ) %>%
  dplyr::group_by(
    process_label,
    human_variable,
    support_measure
  ) %>%
  dplyr::mutate(
    relative_support = if(max(count,na.rm = TRUE) > 0) {
      count / max(count,na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  dplyr::ungroup()

# VISUALISATION: support across relative Q05-Q95 classes ----
human_pressure_class_support_plot_60 <- ggplot2::ggplot(
  human_pressure_class_plot_data_60,
  ggplot2::aes(
    x = relative_class,
    y = human_variable,
    fill = relative_support)) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = count),
    size = 2.8) +
  ggplot2::facet_grid(
    support_measure ~ process_label) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "darkseagreen4",
    limits = c(0,1),
    breaks = c(0,0.25,0.50,0.75,1),
    labels = scales::percent_format(accuracy = 1),name = "Relative\nsupport"
  ) +
  ggplot2::labs(
    x = "Relative class between empirical Q05 and Q95",
    y = "Human-pressure variable",
    title = "Empirical support across relative human-pressure classes",
    subtitle = "Numbers show absolute counts; colour is relative within each variable, measure and process") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "right")

print(human_pressure_class_support_plot_60)




# ------------------------------------------------------------------------------ STEP 4: select model covariates ----
#' **Steps:**
#' (i) define six biologically meaningful environmental adjustment structures;
#' 
#'  (ia) backbone + building + elevation : elevation is assumed to be the main parameter
#'  cofounding of building location 
#'  (ib) backbone + building + ruggedness : humans does not build in highly rugged 
#'  terrain
#'  (ic) backbone + building + elevation + ruggedness : elevation and ruggedness are 
#'  assumed to influence human settlements and therefore may be cofounders. 
#'  (id) backbone + building + elevation + proportion of forest : elevation controls
#'  for settlements position and proportion of forest may influence landing behaviors
#'  (ie) backbone + building + elevation + open habitat : open habitat may positively 
#'  influence flight behaviors by favoring uplifts. 
#'  (if) backbone + building + elevation + open habitat + ruggedness : ruggedness and 
#'  open habitats influence both uplifts and therefore the probability of 
#'  remaining aerial.
#'  (ig) backbone + building + elevation + distance to ridgeline : version 
#'  controlling only for orographic uplifts which are more used at an early age 
#'  (cf. Nourani et al., Elife, 2024)
#'  
#' (ii) calculate one common raw Q05-Q95 HFI contrast supported by all three processes;
#' (iii) fit candidate weighted binomial GLMMs for:
#'       a) GPS aerial versus terrestrial;
#'       b) ACC aerial versus terrestrial;
#'       c) ACC feeding versus resting;
#' (iv) compile AIC, delta AIC, HFI coefficients and Q05-Q95 log-odds contrasts;
#' (v) save the model-comparison tables as a two-page A3 landscape PDF.

# Parameters ----
selected_hfi_candidates_60 <- c( "population_density","settlement_density")
hfi_metric_labels_60 <- c(population_density = "Population density",settlement_density = "Settlement density")
model_delta_aic_threshold_60 <- 4

process_labels_60 <- c(GPS_aerial_vs_terrestrial = "GPS: aerial vs terrestrial",
  ACC_aerial_vs_terrestrial = "ACC: aerial vs terrestrial",
  ACC_feeding_vs_resting = "ACC: feeding vs resting")

model_structure_labels_60 <- c(
  elevation = "Elevation",
  ruggedness = "Ruggedness",
  elevation_ruggedness = "Elevation + ruggedness",
  forest = "Elevation + forest",
  open_habitat = "Elevation + open habitat",
  open_habitat_ruggedness = "Elevation + open habitat + ruggedness",
  ridgeline = "Elevation + distance to ridgeline")

hfi_metric_labels_60 <- c( population_density = "Population density",settlement_density = "Settlement density")

# 4.1 Define environmental adjustment structures ----
model_covariate_structures_60 <- list(
  elevation = "elevation_100m_z",
  ruggedness = "ruggedness_100m_z",
  elevation_ruggedness = c("elevation_100m_z","ruggedness_100m_z"),
  forest = c("elevation_100m_z","prop_forest_5cells_z"),
  open_habitat = c("elevation_100m_z","prop_low_vegetation_5cells_z"),
  open_habitat_ruggedness = c("elevation_100m_z","prop_low_vegetation_5cells_z","ruggedness_100m_z"),
  ridgeline = c("elevation_100m_z","distance_to_ridgeline_100m_z"))

model_covariate_metadata_60 <- tidyr::expand_grid(
  process = names(backbone_process_data_60),
  hfi_variable = selected_hfi_candidates_60,
  model_structure = names(model_covariate_structures_60)) %>%
  dplyr::mutate(
    hfi_term = paste0(hfi_variable,"_z"),
    model_id = paste(process,hfi_variable,model_structure,sep = "__"))

# 4.2 Calculate process-specific raw Q05 and Q95 ----
selected_hfi_candidates_60 <- c("population_density","settlement_density")

hfi_thresholds_60 <- dplyr::bind_rows(
  lapply(names(backbone_process_data_60),function(process_name) {
    data_process <- backbone_process_data_60[[process_name]]
    
    dplyr::bind_rows(
      lapply(selected_hfi_candidates_60,function(hfi_variable) {
        values <- data_process[[hfi_variable]]
        values <- values[is.finite(values)]
        
        tibble::tibble(
          process = process_name,
          hfi_variable = hfi_variable,
          hfi_q05 = as.numeric(stats::quantile(values,0.05,names = FALSE)),
          hfi_q95 = as.numeric(stats::quantile(values,0.95,names = FALSE)))})
    )
  })
)

# 4.2 Define a common raw HFI Q05-Q95 contrast ----
# for defining q05 and q95, we take the smallest value of q95 and the highest 
# value of q05 for ACC and GPS. This allow to compare the contrast on the same
# scale
common_hfi_thresholds_60 <- hfi_thresholds_60 %>%
  dplyr::filter(as.character(hfi_variable) %in% selected_hfi_candidates_60) %>%
  dplyr::mutate(hfi_variable = as.character(hfi_variable)) %>%
  dplyr::group_by(hfi_variable) %>%
  dplyr::summarise(
    hfi_q05_raw = max(hfi_q05,na.rm = TRUE),
    hfi_q95_raw = min(hfi_q95,na.rm = TRUE),
    .groups = "drop")

derive_hfi_scaling_60 <- function(data,process_name,hfi_variable) {
  hfi_term <- paste0(hfi_variable,"_z")
  raw_value <- data[[hfi_variable]]
  standardized_value <- data[[hfi_term]]
  scale_value <- stats::cov(raw_value,standardized_value) / stats::var(standardized_value)
  center_value <- mean(raw_value) - scale_value * mean(standardized_value)
  tibble::tibble(process = process_name,hfi_variable = hfi_variable,
                 center = center_value,scale = scale_value)
}

hfi_process_scaling_60 <- dplyr::bind_rows(lapply(names(backbone_process_data_60),function(process_name) {
  dplyr::bind_rows(lapply(selected_hfi_candidates_60,function(hfi_variable) {
    derive_hfi_scaling_60(backbone_process_data_60[[process_name]],process_name,hfi_variable)
  }))
}))

hfi_contrast_scaling_60 <- tidyr::expand_grid(
  process = names(backbone_process_data_60),
  hfi_variable = selected_hfi_candidates_60) %>%
  dplyr::left_join(common_hfi_thresholds_60,by = "hfi_variable") %>%
  dplyr::left_join(hfi_process_scaling_60,by = c("process","hfi_variable")) %>%
  dplyr::mutate(
    hfi_q05_z = (hfi_q05_raw - center) / scale,
    hfi_q95_z = (hfi_q95_raw - center) / scale,
    hfi_q95_q05_difference_z = hfi_q95_z - hfi_q05_z)

# 4.3 Build process-specific candidate formulas ----
candidate_formulas_60 <- stats::setNames(
  lapply(seq_len(nrow(model_covariate_metadata_60)),function(i) {
    process_name <- model_covariate_metadata_60$process[[i]]
    additional_terms <- c(
      model_covariate_metadata_60$hfi_term[[i]],
      model_covariate_structures_60[[model_covariate_metadata_60$model_structure[[i]]]])
    stats::update.formula(
      selected_backbone_formulas_60[[process_name]],
      stats::as.formula(paste(". ~ . +",paste(additional_terms,collapse = " + "))))
  }),
  model_covariate_metadata_60$model_id)

# 4.4 Fit all weighted binomial GLMMs ----
fit_candidate_model_60 <- function(model_id) {
  model_information <- model_covariate_metadata_60 %>%
    dplyr::filter(.data$model_id == .env$model_id)
  lme4::glmer(
    formula = candidate_formulas_60[[model_id]],
    data = backbone_process_data_60[[model_information$process[[1]]]],
    weights = backbone_weight,
    family = stats::binomial(link = "logit"),
    nAGQ = 1,
    control = lme4::glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5)))}

candidate_models_60 <- stats::setNames(
  lapply(model_covariate_metadata_60$model_id,fit_candidate_model_60),
  model_covariate_metadata_60$model_id)

# 4.5 Extract model fit and HFI effects ----
model_comparison_table_60 <- dplyr::bind_rows(
  lapply(names(candidate_models_60),function(model_id) {
    fitted_model <- candidate_models_60[[model_id]]
    model_information <- model_covariate_metadata_60 %>%
      dplyr::filter(.data$model_id == .env$model_id)
    process_name <- model_information$process[[1]]
    hfi_variable <- model_information$hfi_variable[[1]]
    hfi_term <- model_information$hfi_term[[1]]
    coefficient_table <- summary(fitted_model)$coefficients
    contrast <- hfi_contrast_scaling_60 %>%
      dplyr::filter(.data$process == .env$process_name,
                    .data$hfi_variable == .env$hfi_variable)
    model_loglik <- stats::logLik(fitted_model)
    convergence_message <- fitted_model@optinfo$conv$lme4$messages
    optimizer_code <- fitted_model@optinfo$conv$opt
    
    if(!hfi_term %in% rownames(coefficient_table)) {
      stop("HFI term ",hfi_term," was not found in model ",model_id,".")
    }
    if(nrow(contrast) != 1L) {
      stop("Exactly one HFI contrast was expected for model ",model_id,".")
    }
    
    hfi_coefficient <- coefficient_table[hfi_term,"Estimate"]
    hfi_standard_error <- coefficient_table[hfi_term,"Std. Error"]
    hfi_confidence_low <- hfi_coefficient - 1.96 * hfi_standard_error
    hfi_confidence_high <- hfi_coefficient + 1.96 * hfi_standard_error
    hfi_difference_z <- contrast$hfi_q95_q05_difference_z[[1]]
    
    tibble::tibble(
      model_id = model_id,
      process = process_name,
      hfi_variable = hfi_variable,
      model_structure = model_information$model_structure[[1]],
      n_observations = stats::nobs(fitted_model),
      model_df = attr(model_loglik,"df"),
      AIC = stats::AIC(fitted_model),
      hfi_coefficient = hfi_coefficient,
      hfi_standard_error = hfi_standard_error,
      hfi_confidence_low = hfi_confidence_low,
      hfi_confidence_high = hfi_confidence_high,
      hfi_p_value = coefficient_table[hfi_term,"Pr(>|z|)"],
      Q95_Q05_log_odds = hfi_coefficient * hfi_difference_z,
      Q95_Q05_log_odds_CI_low = hfi_confidence_low * hfi_difference_z,
      Q95_Q05_log_odds_CI_high = hfi_confidence_high * hfi_difference_z,
      converged = is.null(convergence_message) &&
        (is.null(optimizer_code) || all(optimizer_code == 0)),
      singular = lme4::isSingular(fitted_model,tol = 1e-4))
  })) %>%
  dplyr::group_by(process,hfi_variable) %>%
  dplyr::arrange(AIC,.by_group = TRUE) %>%
  dplyr::mutate(
    delta_AIC = AIC - min(AIC,na.rm = TRUE),
    akaike_weight = exp(-0.5 * delta_AIC) /
      sum(exp(-0.5 * delta_AIC)),
    supported = delta_AIC <= model_delta_aic_threshold_60
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(process) %>%
  dplyr::mutate(
    delta_AIC_process = AIC - min(AIC,na.rm = TRUE),
    akaike_weight_process = exp(-0.5 * delta_AIC_process) /
      sum(exp(-0.5 * delta_AIC_process))
  ) %>%
  dplyr::ungroup()
  
  
  

model_comparability_60 <- model_comparison_table_60 %>%
  dplyr::group_by(process) %>%
  dplyr::summarise(
    n_distinct_sample_sizes = dplyr::n_distinct(n_observations),
    sample_sizes = paste(sort(unique(n_observations)),collapse = ", "),
    .groups = "drop"
  )


# 4.6 Store separate process tables ----
gps_model_comparison_table_60 <- model_comparison_table_60 %>%
  dplyr::filter(process == "GPS_aerial_vs_terrestrial")

acc_aerial_model_comparison_table_60 <- model_comparison_table_60 %>%
  dplyr::filter(process == "ACC_aerial_vs_terrestrial")

acc_feeding_resting_model_comparison_table_60 <- model_comparison_table_60 %>%
  dplyr::filter(process == "ACC_feeding_vs_resting")


# 4.7 Select the four best models per process and human-pressure metric ----
human_pressure_metric_labels_60 <- c(
  population_density = "Population density",
  settlement_density = "Settlement density"
)

top_four_models_60 <- model_comparison_table_60 %>%
  dplyr::group_by(process,hfi_variable) %>%
  dplyr::slice_min(
    order_by = AIC,
    n = 4,
    with_ties = FALSE
  ) %>%
  dplyr::arrange(AIC,.by_group = TRUE) %>%
  dplyr::mutate(rank_within_metric = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(process) %>%
  dplyr::arrange(AIC,.by_group = TRUE) %>%
  dplyr::mutate(
    row_order = dplyr::row_number(),
    best_model = row_order == 1L
  ) %>%
  dplyr::ungroup()

# 4.8 Prepare one common landscape table ----
# 4.8 Prepare one readable common table ----
human_pressure_metric_labels_60 <- c(
  population_density = "Population\ndensity",
  settlement_density = "Settlement\ndensity"
)

model_structure_table_labels_60 <- c(
  elevation = "Elevation",
  ruggedness = "Ruggedness",
  elevation_ruggedness = "Elevation +\nruggedness",
  forest = "Elevation +\nforest",
  open_habitat = "Elevation +\nopen habitat",
  open_habitat_ruggedness = "Elevation + open habitat\n+ ruggedness",
  ridgeline = "Elevation + distance\nto ridgeline"
)

top_four_model_table_data_60 <- top_four_models_60 %>%
  dplyr::mutate(
    process_label = unname(process_labels_60[process]),
    metric_label = unname(human_pressure_metric_labels_60[hfi_variable]),
    model_structure_label = unname(model_structure_table_labels_60[model_structure]),
    rank_display = as.character(row_order),
    AIC_display = sprintf("%.1f",AIC),
    delta_AIC_display = sprintf("%.2f",delta_AIC_process),
    akaike_weight_display = sprintf("%.3f",akaike_weight_process),
    coefficient_display = sprintf(
      "%.3f\n[%.3f, %.3f]",
      hfi_coefficient,
      hfi_confidence_low,
      hfi_confidence_high
    ),
    contrast_display = sprintf(
      "%.3f\n[%.3f, %.3f]",
      Q95_Q05_log_odds,
      Q95_Q05_log_odds_CI_low,
      Q95_Q05_log_odds_CI_high
    ),
    convergence_display = dplyr::case_when(
      !converged ~ "No",
      singular ~ "Singular",
      TRUE ~ "Yes"
    ),
    df_display = sprintf("%.0f",model_df)
  ) %>%
  dplyr::select(
    process_label,
    row_order,
    best_model,
    Rank = rank_display,
    Metric = metric_label,
    `Environmental\nadjustment` = model_structure_label,
    AIC = AIC_display,
    `Delta\nAIC` = delta_AIC_display,
    `Akaike\nweight` = akaike_weight_display,
    `Coefficient\n[95% CI]` = coefficient_display,
    `Q05-Q95 log-odds\n[95% CI]` = contrast_display,
    Converged = convergence_display,
    df = df_display
  ) %>%
  tidyr::pivot_longer(
    cols = -c(process_label,row_order,best_model),
    names_to = "table_column",
    values_to = "table_value"
  ) %>%
  dplyr::mutate(
    table_column = factor(
      table_column,
      levels = c(
        "Rank",
        "Metric",
        "Environmental\nadjustment",
        "AIC",
        "Delta\nAIC",
        "Akaike\nweight",
        "Coefficient\n[95% CI]",
        "Q05-Q95 log-odds\n[95% CI]",
        "Converged",
        "df"
      )
    ),
    row_order = factor(
      row_order,
      levels = rev(seq_len(max(row_order)))
    )
  )

top_four_model_table_plot_60 <- ggplot2::ggplot(
  top_four_model_table_data_60,
  ggplot2::aes(
    x = table_column,
    y = row_order
  )
) +
  ggplot2::geom_tile(
    ggplot2::aes(fill = best_model),
    color = "grey72",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = table_value),
    size = 3.1,
    lineheight = 0.85
  ) +
  ggplot2::facet_wrap(
    ~ process_label,
    ncol = 1,
    scales = "free_y"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      `FALSE` = "white",
      `TRUE` = "grey82"
    ),
    guide = "none"
  ) +
  ggplot2::scale_x_discrete(
    position = "top",
    expand = c(0,0)
  ) +
  ggplot2::scale_y_discrete(
    expand = c(0,0)
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    title = "Comparison of population and settlement density models",
    subtitle = paste0(
      "Four lowest-AIC models are shown for each behavioural process. ",
      "Delta AIC and Akaike weights compare both human-pressure metrics. ",
      "Grey indicates the overall lowest-AIC model within each process."
    )
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = 8.5,
      face = "bold",
      lineheight = 0.85,
      hjust = 0.5,
      margin = ggplot2::margin(b = 8)
    ),
    axis.text.y = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      size = 11,
      face = "bold",
      margin = ggplot2::margin(t = 6,b = 6)
    ),
    strip.background = ggplot2::element_rect(
      fill = "grey92",
      color = "grey65"
    ),
    plot.title = ggplot2::element_text(
      size = 15,
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 9.5,
      lineheight = 1
    ),
    panel.spacing = grid::unit(0.8,"lines"),
    plot.margin = ggplot2::margin(15,15,15,15)
  )

print(top_four_model_table_plot_60)

#' we find that the gps and acc aerial and terrestrial model are performing better
#' using the low vegetation, elevation and ruggedness while the acc feeding versus
#' resting form comprise elevation and distance to ridgeline. 


# ------------------------------------------------------------------------------ STEP 5: sensitivity to environmental adjustment ----
#' **Steps:**
#' (i) for GPS and ACC aerial-versus-terrestrial models, fit the full
#'     elevation + open habitat + ruggedness structure and remove each
#'     environmental term separately;
#' (ii) for the conditional ACC feeding-versus-resting model, fit the full
#'      elevation + distance-to-ridgeline structure and remove each term;
#' (iii) quantify changes in AIC, HFI coefficients and common Q05-Q95
#'       log-odds contrasts relative to the corresponding full model;
#' (iv) save the three sensitivity tables in one A3 landscape PDF.
 
# Parameters ----
hfi_sensitivity_variable_60 <- c("settlement_density")
hfi_sensitivity_term_60 <- paste0(hfi_sensitivity_variable_60,"_z")
confidence_level_60 <- 0.95
confidence_critical_value_60 <- stats::qnorm(1 - (1 - confidence_level_60) / 2)

environmental_term_labels_60 <- c(
  elevation = "Elevation",
  open_habitat = "Open habitat",
  ruggedness = "Ruggedness")

sensitivity_environmental_terms_60 <- list(
  GPS_aerial_vs_terrestrial = c(
    elevation = "elevation_100m_z",
    open_habitat = "prop_low_vegetation_5cells_z",
    ruggedness = "ruggedness_100m_z"),
  ACC_aerial_vs_terrestrial = c(
    elevation = "elevation_100m_z",
    open_habitat = "prop_low_vegetation_5cells_z",
    ruggedness = "ruggedness_100m_z"),
  ACC_feeding_vs_resting = c(
    elevation = "elevation_100m_z",
    ruggedness = "ruggedness_100m_z"))

# 5.1 Define full and reduced structures for each process ----
build_sensitivity_structures_60 <- function(named_terms) {
  structures <- list(full = unname(named_terms))
  for(term_name in names(named_terms)) {
    structures[[paste0("without_",term_name)]] <- unname(named_terms[names(named_terms) != term_name])
  }
  structures
}

sensitivity_structures_by_process_60 <- stats::setNames(
  lapply(sensitivity_environmental_terms_60,build_sensitivity_structures_60),
  names(sensitivity_environmental_terms_60))

sensitivity_model_metadata_60 <- dplyr::bind_rows(
  lapply(names(sensitivity_structures_by_process_60),function(process_name) {
    structure_names <- names(sensitivity_structures_by_process_60[[process_name]])
    removed_keys <- sub("^without_","",structure_names)
    tibble::tibble(
      process = process_name,
      model_structure = structure_names,
      structure_order = seq_along(structure_names),
      removed_variable = dplyr::if_else(
        structure_names == "full",
        "None — full model",
        unname(environmental_term_labels_60[removed_keys])),
      model_id = paste(process_name,structure_names,sep = "__"))
  }))

# 5.2 Build the process-specific sensitivity formulas ----
sensitivity_formulas_60 <- stats::setNames(
  lapply(seq_len(nrow(sensitivity_model_metadata_60)),function(i) {
    process_name <- sensitivity_model_metadata_60$process[[i]]
    structure_name <- sensitivity_model_metadata_60$model_structure[[i]]
    environmental_terms <- sensitivity_structures_by_process_60[[process_name]][[structure_name]]
    additional_terms <- c(hfi_sensitivity_term_60,environmental_terms)
    stats::update.formula(
      selected_backbone_formulas_60[[process_name]],
      stats::as.formula(paste(". ~ . +",paste(additional_terms,collapse = " + "))))
  }),
  sensitivity_model_metadata_60$model_id)

# 5.3 Fit the weighted binomial GLMMs ----
fit_sensitivity_model_60 <- function(model_id) {
  model_information <- sensitivity_model_metadata_60 %>%
    dplyr::filter(.data$model_id == .env$model_id)
  process_name <- model_information$process[[1]]
  lme4::glmer(
    formula = sensitivity_formulas_60[[model_id]],
    data = backbone_process_data_60[[process_name]],
    weights = backbone_weight,
    family = stats::binomial(link = "logit"),
    nAGQ = 1,
    control = control_glmer_60)
}

sensitivity_models_60 <- stats::setNames(
  lapply(sensitivity_model_metadata_60$model_id,fit_sensitivity_model_60),
  sensitivity_model_metadata_60$model_id)

# 5.4 Extract AIC, HFI effects and common Q05-Q95 contrasts ----
sensitivity_model_results_60 <- dplyr::bind_rows(
  lapply(names(sensitivity_models_60),function(model_id) {
    fitted_model <- sensitivity_models_60[[model_id]]
    model_information <- sensitivity_model_metadata_60 %>%
      dplyr::filter(.data$model_id == .env$model_id)
    process_name <- model_information$process[[1]]
    coefficient_table <- summary(fitted_model)$coefficients
    model_loglik <- stats::logLik(fitted_model)
    convergence_message <- fitted_model@optinfo$conv$lme4$messages
    optimizer_code <- fitted_model@optinfo$conv$opt
    contrast <- hfi_contrast_scaling_60 %>%
      dplyr::filter(
        .data$process == .env$process_name,
        .data$hfi_variable == .env$hfi_sensitivity_variable_60)
    
    if(!hfi_sensitivity_term_60 %in% rownames(coefficient_table)) {
      stop("HFI term ",hfi_sensitivity_term_60," was not found in model ",model_id,".")
    }
    if(nrow(contrast) != 1L) {
      stop("Exactly one common HFI contrast was expected for model ",model_id,".")
    }
    
    hfi_coefficient <- coefficient_table[hfi_sensitivity_term_60,"Estimate"]
    hfi_standard_error <- coefficient_table[hfi_sensitivity_term_60,"Std. Error"]
    hfi_confidence_low <- hfi_coefficient - confidence_critical_value_60 * hfi_standard_error
    hfi_confidence_high <- hfi_coefficient + confidence_critical_value_60 * hfi_standard_error
    hfi_difference_z <- contrast$hfi_q95_q05_difference_z[[1]]
    
    tibble::tibble(
      model_id = model_id,
      process = process_name,
      model_structure = model_information$model_structure[[1]],
      structure_order = model_information$structure_order[[1]],
      removed_variable = model_information$removed_variable[[1]],
      n_observations = stats::nobs(fitted_model),
      model_df = attr(model_loglik,"df"),
      AIC = stats::AIC(fitted_model),
      hfi_coefficient = hfi_coefficient,
      hfi_standard_error = hfi_standard_error,
      hfi_confidence_low = hfi_confidence_low,
      hfi_confidence_high = hfi_confidence_high,
      hfi_p_value = coefficient_table[hfi_sensitivity_term_60,"Pr(>|z|)"],
      Q95_Q05_difference_z = hfi_difference_z,
      Q95_Q05_log_odds = hfi_coefficient * hfi_difference_z,
      Q95_Q05_log_odds_CI_low = hfi_confidence_low * hfi_difference_z,
      Q95_Q05_log_odds_CI_high = hfi_confidence_high * hfi_difference_z,
      converged = is.null(convergence_message) &&
        (is.null(optimizer_code) || all(optimizer_code == 0)),
      singular = lme4::isSingular(fitted_model,tol = 1e-4))
  }))

# 5.5 Calculate changes relative to each full model ----
hfi_coefficient_sensitivity_60 <- sensitivity_model_results_60 %>%
  dplyr::group_by(process) %>%
  dplyr::mutate(
    full_model_AIC = AIC[model_structure == "full"][[1]],
    full_hfi_coefficient = hfi_coefficient[model_structure == "full"][[1]],
    full_Q95_Q05_log_odds = Q95_Q05_log_odds[model_structure == "full"][[1]],
    delta_AIC_from_full = AIC - full_model_AIC,
    hfi_coefficient_change = hfi_coefficient - full_hfi_coefficient,
    Q95_Q05_log_odds_change = Q95_Q05_log_odds - full_Q95_Q05_log_odds,
    hfi_sign_changed = model_structure != "full" &
      sign(hfi_coefficient) != sign(full_hfi_coefficient)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(process_label = unname(process_labels_60[process])) %>%
  dplyr::arrange(process,structure_order)

# 5.6 Prepare compact process-specific tables ----
prepare_sensitivity_table_60 <- function(process_name) {
  hfi_coefficient_sensitivity_60 %>%
    dplyr::filter(process == process_name) %>%
    dplyr::arrange(structure_order) %>%
    dplyr::transmute(
      `Variable removed` = removed_variable,
      AIC = sprintf("%.1f",AIC),
      `Delta AIC\nvs full` = sprintf("%+.2f",delta_AIC_from_full),
      `HFI coefficient\n[95% CI]` = sprintf(
        "%.4f\n[%.4f, %.4f]",
        hfi_coefficient,hfi_confidence_low,hfi_confidence_high),
      `Change in HFI\ncoefficient` = sprintf("%+.4f",hfi_coefficient_change),
      `Q05-Q95 log-odds\n[95% CI]` = sprintf(
        "%.4f\n[%.4f, %.4f]",
        Q95_Q05_log_odds,Q95_Q05_log_odds_CI_low,Q95_Q05_log_odds_CI_high),
      `Change in Q05-Q95\nlog-odds` = sprintf("%+.4f",Q95_Q05_log_odds_change),
      `Sign changed` = dplyr::if_else(
        model_structure == "full","—",
        dplyr::if_else(hfi_sign_changed,"Yes","No")))
}

gps_sensitivity_table_60 <- prepare_sensitivity_table_60(
  "GPS_aerial_vs_terrestrial")

acc_aerial_sensitivity_table_60 <- prepare_sensitivity_table_60(
  "ACC_aerial_vs_terrestrial")

acc_feeding_resting_sensitivity_table_60 <- prepare_sensitivity_table_60(
  "ACC_feeding_vs_resting")

# 5.7 Convert one data frame into a table plot ----
make_sensitivity_table_plot_60 <- function(data,title,subtitle,text_size = 3.1) {
  table_plot_data <- data %>%
    dplyr::mutate(
      row_id = dplyr::row_number(),
      row_fill = dplyr::if_else(row_id %% 2L == 0L,"grey96","white")) %>%
    tidyr::pivot_longer(
      cols = -c(row_id,row_fill),
      names_to = "column",
      values_to = "value") %>%
    dplyr::mutate(
      column = factor(column,levels = names(data)),
      row_id = factor(row_id,levels = rev(seq_len(nrow(data)))))
  
  ggplot2::ggplot(table_plot_data,ggplot2::aes(x = column,y = row_id)) +
    ggplot2::geom_tile(ggplot2::aes(fill = row_fill),color = "grey80") +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(ggplot2::aes(label = value),size = text_size,lineheight = 0.9) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL,y = NULL,title = title,subtitle = subtitle) +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 8.5,face = "bold",lineheight = 0.9,
        angle = 20,hjust = 0,margin = ggplot2::margin(b = 8)),
      plot.title = ggplot2::element_text(face = "bold",size = 14),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(15,15,15,15))
}

# 5.8 Create the three sensitivity table plots ----
gps_sensitivity_table_plot_60 <- make_sensitivity_table_plot_60(
  gps_sensitivity_table_60,
  "GPS: sensitivity of the HFI effect",
  "Aerial versus terrestrial; full model = elevation + open habitat + ruggedness")

acc_aerial_sensitivity_table_plot_60 <- make_sensitivity_table_plot_60(
  acc_aerial_sensitivity_table_60,
  "ACC: sensitivity of the HFI effect",
  "Aerial versus terrestrial; full model = elevation + open habitat + ruggedness")

acc_feeding_resting_sensitivity_table_plot_60 <- make_sensitivity_table_plot_60(
  acc_feeding_resting_sensitivity_table_60,
  "ACC conditional model: sensitivity of the HFI effect",
  "Feeding versus resting; full model = elevation + ruggedness")

# 5.9 Save the three tables in an A3 landscape PDF ----
sensitivity_pdf_file_60 <- file.path(
  results_directory_60,
  "STEP5_HFI_environmental_sensitivity_A3_landscape.pdf")

grDevices::pdf(
  file = sensitivity_pdf_file_60,
  width = 16.54,
  height = 11.69,
  onefile = TRUE,
  family = "Helvetica")

print(gps_sensitivity_table_plot_60)
print(acc_aerial_sensitivity_table_plot_60)
print(acc_feeding_resting_sensitivity_table_plot_60)
grDevices::dev.off()

message("Sensitivity PDF saved to: ",sensitivity_pdf_file_60)



# ------------------------------------------------------------------------------ STEP 6: VIF of retained models ----
#' Check that all VIF values remain below 3.
hfi_vif_term_60 <- "settlement_density"

vif_model_predictors_60 <- list(
  GPS_aerial_vs_terrestrial = c(hfi_vif_term_60,"elevation_100m_z","ruggedness_100m_z"),
  ACC_aerial_vs_terrestrial = c(hfi_vif_term_60,"elevation_100m_z","ruggedness_100m_z"),
  ACC_feeding_vs_resting = c(hfi_vif_term_60,"elevation_100m_z", "ruggedness_100m_z"))

calculate_vif_60 <- function(data,predictors) {
  dplyr::bind_rows(lapply(predictors,function(variable) {
    auxiliary_model <- stats::lm(stats::reformulate(setdiff(predictors,variable),response = variable),data = data)
    tibble::tibble(variable = variable,vif = 1 / (1 - summary(auxiliary_model)$r.squared))
  }))
}

vif_selected_models_60 <- dplyr::bind_rows(lapply(names(vif_model_predictors_60),function(process_name) {
  calculate_vif_60(
    data = backbone_process_data_60[[process_name]],
    predictors = vif_model_predictors_60[[process_name]]
  ) %>%
    dplyr::mutate(process = process_name,.before = 1)
})) %>%
  dplyr::group_by(process) %>%
  dplyr::mutate(maximum_vif = max(vif),all_vif_below_3 = maximum_vif < 3) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(process,dplyr::desc(vif))

print(vif_selected_models_60,n = Inf,width = Inf)


# CONCLUSION ----
# Retained models:
# GPS aerial vs terrestrial: HFI + elevation + ruggedness
# ACC aerial vs terrestrial: HFI + elevation + ruggedness
# ACC feeding vs resting: HFI + elevation

# Parameters ----
final_hfi_variable_60 <- "settlement_density"
final_hfi_term_60 <- paste0(final_hfi_variable_60,"_z")
final_environmental_terms_60 <- list(
  GPS_aerial_vs_terrestrial = c("elevation_100m_z","ruggedness_100m_z"),
  ACC_aerial_vs_terrestrial = c("elevation_100m_z","ruggedness_100m_z"),
  ACC_feeding_vs_resting = c("elevation_100m_z", "ruggedness_100m_z"))
dir.create(results_directory_60,recursive = TRUE,showWarnings = FALSE)

# 1. Define and fit the retained models ----
final_model_formulas_60 <- stats::setNames(lapply(names(final_environmental_terms_60),function(process_name) {
  additional_terms <- c(final_hfi_term_60,final_environmental_terms_60[[process_name]])
  stats::update.formula(
    selected_backbone_formulas_60[[process_name]],
    stats::as.formula(paste(". ~ . +",paste(additional_terms,collapse = " + "))))
}),names(final_environmental_terms_60))

fit_final_model_60 <- function(process_name) {
  lme4::glmer(
    formula = final_model_formulas_60[[process_name]],
    data = backbone_process_data_60[[process_name]],
    weights = backbone_weight,
    family = stats::binomial(link = "logit"),
    nAGQ = 1,
    control = control_glmer_60)
}

final_models_60 <- stats::setNames(
  lapply(names(final_model_formulas_60),fit_final_model_60),
  names(final_model_formulas_60))

final_model_gps_60 <- final_models_60[["GPS_aerial_vs_terrestrial"]]
final_model_acc_aerial_60 <- final_models_60[["ACC_aerial_vs_terrestrial"]]
final_model_acc_feeding_resting_60 <- final_models_60[["ACC_feeding_vs_resting"]]

control_final_models_60 <- dplyr::bind_rows(lapply(names(final_models_60),function(process_name) {
  fitted_model <- final_models_60[[process_name]]
  convergence_message <- fitted_model@optinfo$conv$lme4$messages
  optimizer_code <- fitted_model@optinfo$conv$opt
  tibble::tibble(
    process = process_name,n_observations = stats::nobs(fitted_model),
    AIC = stats::AIC(fitted_model),
    converged = is.null(convergence_message) &&
      (is.null(optimizer_code) || all(optimizer_code == 0)),
    singular = lme4::isSingular(fitted_model,tol = 1e-4))
}))

# 2. Extract the common raw Q05-Q95 HFI contrast ----
final_common_hfi_thresholds_60 <- common_hfi_thresholds_60 %>%
  dplyr::filter(hfi_variable == final_hfi_variable_60)

if(nrow(final_common_hfi_thresholds_60) != 1L) {
  stop("Exactly one common HFI contrast was expected.")
}

get_final_hfi_levels_60 <- function(process_name) {
  levels <- hfi_contrast_scaling_60 %>%
    dplyr::filter(
      .data$process == .env$process_name,
      .data$hfi_variable == .env$final_hfi_variable_60)
  if(nrow(levels) != 1L) stop("One HFI contrast was expected for ",process_name,".")
  c(Q05 = levels$hfi_q05_z[[1]],Q95 = levels$hfi_q95_z[[1]])
}

predict_fixed_probability_60 <- function(model,data,hfi_value_z) {
  prediction_data <- data
  prediction_data[[final_hfi_term_60]] <- hfi_value_z
  as.numeric(stats::predict(
    model,newdata = prediction_data,type = "response",
    re.form = NA,allow.new.levels = TRUE))
}

# 3. Predict GPS aerial and terrestrial probabilities ----
gps_prediction_data_60 <- backbone_process_data_60[["GPS_aerial_vs_terrestrial"]]
gps_hfi_levels_60 <- get_final_hfi_levels_60("GPS_aerial_vs_terrestrial")
gps_p_aerial_q05_60 <- predict_fixed_probability_60(
  final_model_gps_60,gps_prediction_data_60,gps_hfi_levels_60[["Q05"]])
gps_p_aerial_q95_60 <- predict_fixed_probability_60(
  final_model_gps_60,gps_prediction_data_60,gps_hfi_levels_60[["Q95"]])

gps_probability_rows_60 <- dplyr::bind_rows(
  tibble::tibble(
    dataset = "GPS",individual_id = gps_prediction_data_60$individual_id,
    hfi_level = "Q05",aerial = gps_p_aerial_q05_60,
    terrestrial = 1 - gps_p_aerial_q05_60),
  tibble::tibble(
    dataset = "GPS",individual_id = gps_prediction_data_60$individual_id,
    hfi_level = "Q95",aerial = gps_p_aerial_q95_60,
    terrestrial = 1 - gps_p_aerial_q95_60)) %>%
  tidyr::pivot_longer(
    cols = c(aerial,terrestrial),
    names_to = "state",values_to = "probability")

# 4. Predict and combine the two conditional ACC models ----
acc_prediction_data_60 <- backbone_process_data_60[["ACC_aerial_vs_terrestrial"]]
acc_aerial_hfi_levels_60 <- get_final_hfi_levels_60("ACC_aerial_vs_terrestrial")
acc_conditional_hfi_levels_60 <- get_final_hfi_levels_60("ACC_feeding_vs_resting")

acc_p_aerial_q05_60 <- predict_fixed_probability_60(
  final_model_acc_aerial_60,acc_prediction_data_60,
  acc_aerial_hfi_levels_60[["Q05"]])
acc_p_aerial_q95_60 <- predict_fixed_probability_60(
  final_model_acc_aerial_60,acc_prediction_data_60,
  acc_aerial_hfi_levels_60[["Q95"]])

acc_p_feeding_given_terrestrial_q05_60 <- predict_fixed_probability_60(
  final_model_acc_feeding_resting_60,acc_prediction_data_60,
  acc_conditional_hfi_levels_60[["Q05"]])
acc_p_feeding_given_terrestrial_q95_60 <- predict_fixed_probability_60(
  final_model_acc_feeding_resting_60,acc_prediction_data_60,
  acc_conditional_hfi_levels_60[["Q95"]])

acc_probability_rows_60 <- dplyr::bind_rows(
  tibble::tibble(
    dataset = "ACC",individual_id = acc_prediction_data_60$individual_id,
    hfi_level = "Q05",
    aerial = acc_p_aerial_q05_60,
    feeding = (1 - acc_p_aerial_q05_60) *
      acc_p_feeding_given_terrestrial_q05_60,
    resting = (1 - acc_p_aerial_q05_60) *
      (1 - acc_p_feeding_given_terrestrial_q05_60)),
  tibble::tibble(
    dataset = "ACC",individual_id = acc_prediction_data_60$individual_id,
    hfi_level = "Q95",
    aerial = acc_p_aerial_q95_60,
    feeding = (1 - acc_p_aerial_q95_60) *
      acc_p_feeding_given_terrestrial_q95_60,
    resting = (1 - acc_p_aerial_q95_60) *
      (1 - acc_p_feeding_given_terrestrial_q95_60))) %>%
  tidyr::pivot_longer(
    cols = c(aerial,feeding,resting),
    names_to = "state",values_to = "probability")

# 5. Average predictions equally among individuals ----
average_transition_probabilities_60 <- function(data) {
  data %>%
    dplyr::group_by(dataset,individual_id,hfi_level,state) %>%
    dplyr::summarise(probability = mean(probability),.groups = "drop") %>%
    dplyr::group_by(dataset,hfi_level,state) %>%
    dplyr::summarise(probability = mean(probability),.groups = "drop")
}

final_transition_probabilities_60 <- dplyr::bind_rows(
  average_transition_probabilities_60(gps_probability_rows_60),
  average_transition_probabilities_60(acc_probability_rows_60)) %>%
  dplyr::mutate(
    hfi_level = factor(hfi_level,levels = c("Q05","Q95")),
    dataset_label = dplyr::recode(
      dataset,GPS = "GPS behavioural classification",
      ACC = "ACC behavioural classification"),
    response = dplyr::case_when(
      dataset == "GPS" & state == "aerial" ~ "Remain aerial",
      dataset == "GPS" & state == "terrestrial" ~ "Transition to terrestrial",
      dataset == "ACC" & state == "aerial" ~ "Remain aerial",
      dataset == "ACC" & state == "feeding" ~ "Transition to feeding",
      dataset == "ACC" & state == "resting" ~ "Transition to resting"),
    response_order = dplyr::case_when(
      response == "Remain aerial" ~ 1L,
      response == "Transition to terrestrial" ~ 2L,
      response == "Transition to resting" ~ 3L,
      response == "Transition to feeding" ~ 4L)) %>%
  dplyr::arrange(dataset,response_order,hfi_level)

# 6. Compile the final probability summary table ----
confidence_critical_value_60 <- stats::qnorm(0.975)

extract_final_hfi_effect_60 <- function(process_name,dataset,response) {
  fitted_model <- final_models_60[[process_name]]
  coefficient_table <- summary(fitted_model)$coefficients
  contrast <- hfi_contrast_scaling_60 %>%
    dplyr::filter(
      .data$process == .env$process_name,
      .data$hfi_variable == .env$final_hfi_variable_60
    )
  
  if(!final_hfi_term_60 %in% rownames(coefficient_table)) {
    stop("HFI coefficient absent from model ",process_name,".")
  }
  if(nrow(contrast) != 1L) {
    stop("Exactly one HFI contrast was expected for ",process_name,".")
  }
  
  estimate <- coefficient_table[final_hfi_term_60,"Estimate"]
  standard_error <- coefficient_table[final_hfi_term_60,"Std. Error"]
  confidence_low <- estimate - confidence_critical_value_60 * standard_error
  confidence_high <- estimate + confidence_critical_value_60 * standard_error
  hfi_difference_z <- contrast$hfi_q95_q05_difference_z[[1]]
  
  tibble::tibble(
    dataset = dataset,
    response = response,
    hfi_coefficient = estimate,
    hfi_confidence_low = confidence_low,
    hfi_confidence_high = confidence_high,
    Q95_Q05_log_odds = estimate * hfi_difference_z,
    Q95_Q05_log_odds_CI_low = confidence_low * hfi_difference_z,
    Q95_Q05_log_odds_CI_high = confidence_high * hfi_difference_z
  )
}

retained_hfi_summary_table_60 <- dplyr::bind_rows(
  extract_final_hfi_effect_60(
    process_name = "GPS_aerial_vs_terrestrial",
    dataset = "GPS model",
    response = "Remain aerial versus terrestrial"
  ),
  extract_final_hfi_effect_60(
    process_name = "ACC_aerial_vs_terrestrial",
    dataset = "ACC sequential model",
    response = "Remain aerial versus terrestrial"
  ),
  extract_final_hfi_effect_60(
    process_name = "ACC_feeding_vs_resting",
    dataset = "ACC sequential model",
    response = "Feeding versus resting | terrestrial"
  )
) %>%
  dplyr::mutate(
    hfi_coefficient_CI = sprintf(
      "[%.3f; %.3f]",
      hfi_confidence_low,
      hfi_confidence_high
    ),
    Q95_Q05_log_odds_CI = sprintf(
      "[%.3f; %.3f]",
      Q95_Q05_log_odds_CI_low,
      Q95_Q05_log_odds_CI_high
    )
  ) %>%
  dplyr::select(
    dataset,
    response,
    hfi_coefficient,
    hfi_coefficient_CI,
    Q95_Q05_log_odds,
    Q95_Q05_log_odds_CI
  )

# 7. Create the compact final results table ----
retained_hfi_summary_gt_60 <- retained_hfi_summary_table_60 %>%
  gt::gt(groupname_col = "dataset") %>%
  gt::tab_header(
    title = gt::md("**Settlement density influence on the decision to land in GPS and ACC models**"),
    subtitle = paste0(
      "Settlement density; common raw contrast: Q05 = ",
      sprintf("%.3f",final_common_hfi_thresholds_60$hfi_q05_raw),
      "; Q95 = ",
      sprintf("%.3f",final_common_hfi_thresholds_60$hfi_q95_raw)
    )
  ) %>%
  gt::cols_label(
    response = "Model response",
    hfi_coefficient = "Settlement coeff",
    hfi_coefficient_CI = "Settlement coefficient 95% CI",
    Q95_Q05_log_odds = "Q95−Q05 log-odds",
    Q95_Q05_log_odds_CI = "Log-odds 95% CI"
  ) %>%
  gt::fmt_number(
    columns = c(hfi_coefficient,Q95_Q05_log_odds),
    decimals = 3
  ) %>%
  gt::cols_align(
    align = "left",
    columns = c(response,hfi_coefficient_CI,Q95_Q05_log_odds_CI)
  ) %>%
  gt::cols_align(
    align = "center",
    columns = c(hfi_coefficient,Q95_Q05_log_odds)
  ) %>%
  gt::tab_style(
    style = list(
      gt::cell_fill(color = "#ECECF8"),
      gt::cell_text(weight = "bold")
    ),
    locations = gt::cells_column_labels()
  ) %>%
  gt::tab_style(
    style = list(
      gt::cell_fill(color = "#F4F4FA"),
      gt::cell_text(weight = "bold")
    ),
    locations = gt::cells_row_groups()
  ) %>%
  gt::tab_options(
    table.width = gt::pct(100),
    table.font.size = gt::px(14),
    heading.title.font.size = gt::px(19),
    heading.subtitle.font.size = gt::px(13),
    column_labels.font.weight = "bold",
    data_row.padding = gt::px(7)
  ) %>%
  gt::tab_source_note(
    source_note = gt::md(
      paste0(
        "For the first two models, a positive coefficient indicates greater ",
        "aerial persistence. For the conditional ACC model, a positive ",
        "coefficient indicates a shift toward feeding rather than resting ",
        "among terrestrial transitions."
      )
    )
  )

print(retained_hfi_summary_gt_60)



# 9. Export VISUALISATION n°1 and n°2 ----
ggplot2::ggsave(
  filename = file.path(
    results_directory_60,
    "VISUALISATION_1_HFI_absolute_class_support.png"),
  plot = hfi_class_support_plot_60,
  width = 16.54,height = 11.69,units = "in",dpi = 300)

ggplot2::ggsave(
  filename = file.path(
    results_directory_60,
    "VISUALISATION_2_HFI_quantile_individual_support.png"),
  plot = hfi_quantile_individual_support_plot_60,
  width = 16.54,height = 11.69,units = "in",dpi = 300)

message(
  "Final visualisations saved in: ",
  results_directory_60)



# 6. Reconstruct marginal ACC state probabilities ----
acc_prediction_data_60 <- backbone_process_data_60[["ACC_aerial_vs_terrestrial"]]

acc_aerial_hfi_levels_60 <- get_final_hfi_levels_60(
  "ACC_aerial_vs_terrestrial")

acc_conditional_hfi_levels_60 <- get_final_hfi_levels_60(
  "ACC_feeding_vs_resting")

predict_acc_state_probabilities_60 <- function(hfi_level) {
  p_aerial <- predict_fixed_probability_60(
    final_model_acc_aerial_60,
    acc_prediction_data_60,
    acc_aerial_hfi_levels_60[[hfi_level]])
  
  p_feeding_given_terrestrial <- predict_fixed_probability_60(
    final_model_acc_feeding_resting_60,
    acc_prediction_data_60,
    acc_conditional_hfi_levels_60[[hfi_level]])
  
  tibble::tibble(
    individual_id = acc_prediction_data_60$individual_id,
    hfi_level = hfi_level,
    aerial = p_aerial,
    terrestrial = 1 - p_aerial,
    feeding_given_terrestrial = p_feeding_given_terrestrial,
    feeding = (1 - p_aerial) * p_feeding_given_terrestrial,
    resting = (1 - p_aerial) * (1 - p_feeding_given_terrestrial))
}

acc_state_probability_rows_60 <- dplyr::bind_rows(
  predict_acc_state_probabilities_60("Q05"),
  predict_acc_state_probabilities_60("Q95"))

# Average first within individuals, then equally across individuals
acc_state_probabilities_60 <- acc_state_probability_rows_60 %>%
  dplyr::group_by(individual_id,hfi_level) %>%
  dplyr::summarise(
    aerial = mean(aerial),
    terrestrial = mean(terrestrial),
    feeding_given_terrestrial = mean(feeding_given_terrestrial),
    feeding = mean(feeding),
    resting = mean(resting),
    .groups = "drop") %>%
  dplyr::group_by(hfi_level) %>%
  dplyr::summarise(
    aerial = mean(aerial),
    terrestrial = mean(terrestrial),
    feeding_given_terrestrial = mean(feeding_given_terrestrial),
    feeding = mean(feeding),
    resting = mean(resting),
    .groups = "drop")
# 7. Compile Q05-Q95 marginal probability contrasts ----
acc_marginal_probability_summary_60 <- acc_state_probabilities_60 %>%
  dplyr::select(hfi_level,aerial,feeding,resting) %>%
  tidyr::pivot_longer(
    cols = c(aerial,feeding,resting),
    names_to = "state",
    values_to = "probability") %>%
  tidyr::pivot_wider(
    names_from = hfi_level,
    values_from = probability,
    names_prefix = "probability_") %>%
  dplyr::mutate(
    probability_difference_Q95_Q05 =
      probability_Q95 - probability_Q05,
    Q95_Q05_log_odds =
      stats::qlogis(probability_Q95) -
      stats::qlogis(probability_Q05))

print(acc_marginal_probability_summary_60,n = Inf,width = Inf)
# state   probability_Q05 probability_Q95 probability_difference_Q95_Q05 Q95_Q05_log_odds
#   1 aerial           0.371           0.384                         0.0135            0.0576
# 2 feeding          0.0556          0.0460                       -0.00960          -0.200 
# 3 resting          0.573           0.570                        -0.00393          -0.0161
