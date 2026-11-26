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
GE_acc_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/acc_20.rds")
GE_gps_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/gps_20.rds")

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

# for 60 minutes datasets
# process                   duration_form    AIC delta_AIC akaike_weight singular
# ACC_aerial_vs_terrestrial linear        12927.      2.75         0.202 FALSE   
# ACC_feeding_vs_resting    linear         3683.      0            0.665 FALSE   
# GPS_aerial_vs_terrestrial quadratic     14046.      0            0.978 FALSE 

# for 20 minutes datasets
# process                   duration_form    AIC delta_AIC akaike_weight singular
# ACC_aerial_vs_terrestrial quadratic     25369.     0             1.000 FALSE   
# ACC_feeding_vs_resting    linear         5185.     0.945         0.384 FALSE   
# GPS_aerial_vs_terrestrial quadratic     28559.     0             1.000 FALSE  

selected_backbone_formulas_60 <- list(
  GPS_aerial_vs_terrestrial =
    remain_aerial ~ cos_diel_c + sin_diel_c + duration_z + duration_z2 + (1 | individual_id),
  ACC_aerial_vs_terrestrial =
    remain_aerial ~ cos_diel_c + sin_diel_c + duration_z + duration_z2 + (1 | individual_id),
  ACC_feeding_vs_resting =
    feeding_vs_resting ~ cos_diel_c + sin_diel_c + duration_z + (1 | individual_id))


#' ----------------------------------------------------------------------------- STEP 2: inspect proxy of human presence ----
#' **Steps:**
#' (1) inspect Pearson correlations and list the pairs above the threshold;
#' (2) drop redundant environmental covariates on biological grounds;
#' (3) reshape the three human-pressure layers to long format; and compute each layer's empirical Q05 and Q95 within each process;
#' (4) VISUALISATION n°5: represent nbr of transition and individual in the distribution 
#' of the 3 standardized proxies of human presence
#' (5) TABLE n°1: summary of the distribution of values of human variables  

# Parameters ----
correlation_threshold_60 <- 0.70
environmental_covariates_60 <- c(
  "elevation_100m","ruggedness_100m","slope_100m","distance_to_ridgeline_100m",
  "prop_forest_5cells","prop_low_vegetation_5cells","prop_rocky_terrain_5cells",
  "population_density","settlement_density","settlement_density_raw")

process_labels_60 <- c(
  GPS_aerial_vs_terrestrial = "GPS: aerial vs terrestrial",
  ACC_aerial_vs_terrestrial = "ACC: aerial vs terrestrial",
  ACC_feeding_vs_resting    = "ACC: feeding vs resting")

human_pressure_variables_60 <- c("settlement_density","settlement_density_raw",
                                 "population_density")

human_pressure_labels_60 <- c(
  settlement_density     = "Built-area share (25 m)",
  settlement_density_raw = "Built-cell share (100 m)",
  population_density     = "Population density")

relative_class_labels_60 <- c("0-10%","10-20%","20-30%","30-40%","40-50%",
                              "50-60%","60-70%","70-80%","80-90%","90-100%")

#' TABLE 1 holds one dataset type per block. ACC_feeding_vs_resting is a subset
#' of ACC_aerial_vs_terrestrial and would duplicate the ACC block, so only the
#' two aerial-origin processes are tabulated.
table_datasets_60 <- c(ACC = "ACC_aerial_vs_terrestrial",GPS = "GPS_aerial_vs_terrestrial")
table_quantiles_60 <- c(0.05,0.90)
table_proxy_order_60 <- c("settlement_density_raw","settlement_density","population_density")

#' Both settlement layers are stored as percentages from 0 to 100, so their
#' values already read as "% of a km2 that is built". Population density is in
#' residents per km2 and is not a percentage, hence the unit in each label.
table_proxy_labels_60 <- c(
  settlement_density_raw = "Proportion of built 100 m cells per km2",
  settlement_density     = "Proportion of settlement (25 m built cells per 100m) per km2",
  population_density     = "Resident population density per km2")

results_directory_60 <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results"
dir.create(results_directory_60,recursive = TRUE,showWarnings = FALSE)

#' Geometry is never needed from here on and slows every operation below.
step2_process_data_60 <- lapply(backbone_process_data_60,function(data) {
  if(inherits(data,"sf")) sf::st_drop_geometry(data) else data})


# 2.1 Correlations among environmental covariates ----
correlation_matrix_list_60 <- stats::setNames(
  lapply(step2_process_data_60,function(data) {
    matrix_60 <- stats::cor(as.data.frame(data)[,environmental_covariates_60],
                            method = "pearson",use = "pairwise.complete.obs")
    dimnames(matrix_60) <- list(environmental_covariates_60,environmental_covariates_60)
    matrix_60}),names(step2_process_data_60))

correlation_summary_60 <- dplyr::bind_rows(
  lapply(names(correlation_matrix_list_60),function(process_name) {
    matrix_60 <- correlation_matrix_list_60[[process_name]]
    indices <- which(lower.tri(matrix_60),arr.ind = TRUE)
    tibble::tibble(
      process = process_name,
      variable_1 = rownames(matrix_60)[indices[,1]],
      variable_2 = colnames(matrix_60)[indices[,2]],
      pearson_r = matrix_60[indices])})) %>%
  dplyr::mutate(process_label = unname(process_labels_60[process])) %>%
  dplyr::arrange(process,dplyr::desc(abs(pearson_r)))

correlated_pairs_60 <- correlation_summary_60 %>%
  dplyr::filter(abs(pearson_r) >= correlation_threshold_60) %>%
  dplyr::transmute(process,variable_1,variable_2,pearson_r = round(pearson_r,3))

print(as.data.frame(correlated_pairs_60),row.names = FALSE)

invisible(lapply(names(correlation_matrix_list_60),function(process_name) {
  corrplot::corrplot(
    correlation_matrix_list_60[[process_name]],method = "color",type = "lower",
    order = "original",col = corrplot::COL2("RdBu",200),col.lim = c(-1,1),
    diag = FALSE,addCoef.col = "black",number.digits = 2,number.cex = 0.6,
    tl.col = "black",tl.cex = 0.6,tl.srt = 45,cl.cex = 0.8,
    addgrid.col = "white",
    title = paste("Pearson correlations -",process_labels_60[[process_name]]),
    mar = c(0,0,3,0))}))


# 2.2 Retain non-redundant environmental covariates ----
excluded_environmental_covariates_60 <- tibble::tibble(
  variable = c("slope_100m","prop_rocky_terrain_5cells"),
  reason = c("Redundant with ruggedness and less representative of topographic complexity",
             "Associated with elevation and of limited additional biological interpretability"))

selected_environmental_covariates_60 <- setdiff(
  environmental_covariates_60,excluded_environmental_covariates_60$variable)


# 2.3 Reshape the human-pressure layers and compute Q05 and Q95 ----
human_pressure_long_60 <- dplyr::bind_rows(
  lapply(names(step2_process_data_60),function(process_name) {
    step2_process_data_60[[process_name]] %>%
      dplyr::select(individual_id,dplyr::all_of(human_pressure_variables_60)) %>%
      dplyr::mutate(
        process = process_name,
        process_label = unname(process_labels_60[process_name])) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(human_pressure_variables_60),
        names_to = "human_variable",
        values_to = "human_value")})) %>%
  dplyr::filter(is.finite(human_value))

human_pressure_thresholds_60 <- human_pressure_long_60 %>%
  dplyr::group_by(process,process_label,human_variable) %>%
  dplyr::summarise(
    empirical_q05 = as.numeric(stats::quantile(human_value,0.05,names = FALSE)),
    empirical_q95 = as.numeric(stats::quantile(human_value,0.95,names = FALSE)),
    .groups = "drop")


# VISUALISATION n°1 ------------------------------------------------------------
human_pressure_class_support_60 <- human_pressure_long_60 %>%
  dplyr::left_join(human_pressure_thresholds_60,
                   by = c("process","process_label","human_variable")) %>%
  dplyr::filter(
    human_value >= empirical_q05,
    human_value <= empirical_q95,
    empirical_q95 > empirical_q05) %>%
  dplyr::mutate(
    relative_position = (human_value - empirical_q05) / (empirical_q95 - empirical_q05),
    relative_class_number = pmin(10L,floor(relative_position * 10) + 1L),
    relative_class = factor(relative_class_labels_60[relative_class_number],
                            levels = relative_class_labels_60)) %>%
  dplyr::group_by(process,process_label,human_variable,relative_class) %>%
  dplyr::summarise(
    n_transitions = dplyr::n(),
    n_individuals = dplyr::n_distinct(individual_id),
    .groups = "drop") %>%
  dplyr::group_by(process,process_label,human_variable) %>%
  tidyr::complete(
    relative_class = factor(relative_class_labels_60,levels = relative_class_labels_60),
    fill = list(n_transitions = 0L,n_individuals = 0L)) %>%
  dplyr::ungroup()

human_pressure_class_plot_data_60 <- human_pressure_class_support_60 %>%
  dplyr::select(process_label,human_variable,relative_class,n_transitions,n_individuals) %>%
  tidyr::pivot_longer(
    cols = c(n_transitions,n_individuals),
    names_to = "support_measure",
    values_to = "count") %>%
  dplyr::mutate(
    support_measure = dplyr::recode(support_measure,
                                    n_transitions = "Transitions",
                                    n_individuals = "Individuals"),
    support_measure = factor(support_measure,levels = c("Transitions","Individuals")),
    human_variable = factor(unname(human_pressure_labels_60[human_variable]),
                            levels = unname(human_pressure_labels_60[human_pressure_variables_60])),
    relative_class = factor(relative_class,levels = relative_class_labels_60)) %>%
  dplyr::group_by(process_label,human_variable,support_measure) %>%
  dplyr::mutate(
    relative_support = if(max(count,na.rm = TRUE) > 0) count / max(count,na.rm = TRUE) else 0) %>%
  dplyr::ungroup()

human_pressure_class_support_plot_60 <- ggplot2::ggplot(
  human_pressure_class_plot_data_60,
  ggplot2::aes(x = relative_class,y = human_variable,fill = relative_support)) +
  ggplot2::geom_tile(color = "white",linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = count),size = 2.8) +
  ggplot2::facet_grid(support_measure ~ process_label) +
  ggplot2::scale_fill_gradient(
    low = "white",high = "darkseagreen4",limits = c(0,1),
    breaks = c(0,0.25,0.50,0.75,1),
    labels = scales::percent_format(accuracy = 1),name = "Relative\nsupport") +
  ggplot2::labs(
    x = "Relative class between empirical Q05 and Q95",
    y = "Human-pressure variable",
    title = "Empirical support across relative human-pressure classes",
    subtitle = "Numbers show absolute counts; colour is relative within each variable, measure and process") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45,hjust = 1),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "right")

print(human_pressure_class_support_plot_60)


# TABLE 1: summary statistcs of the distribution of human variables ------------
human_pressure_table_data_60 <- human_pressure_long_60 %>%
  dplyr::filter(process %in% table_datasets_60) %>%
  dplyr::mutate(
    dataset_type = names(table_datasets_60)[match(process,table_datasets_60)],
    proxy = factor(unname(table_proxy_labels_60[human_variable]),
                   levels = unname(table_proxy_labels_60[table_proxy_order_60]))) %>%
  dplyr::group_by(dataset_type,proxy) %>%
  dplyr::summarise(
    q05 = as.numeric(stats::quantile(human_value,table_quantiles_60[1],names = FALSE)),
    q90 = as.numeric(stats::quantile(human_value,table_quantiles_60[2],names = FALSE)),
    median = stats::median(human_value),
    .groups = "drop") %>%
  dplyr::arrange(dataset_type,proxy)

human_pressure_table_60 <- human_pressure_table_data_60 %>%
  gt::gt(groupname_col = "dataset_type") %>%
  gt::cols_label(proxy = "Proxy of human presence on the ground",q05 = "Q05",q90 = "Q90",median = "Median") %>%
  gt::fmt_number(columns = c(q05,q90,median),rows = !grepl("Population",proxy),decimals = 2) %>%
  gt::fmt_number(columns = c(q05,q90,median),rows = grepl("Population",proxy),decimals = 1) %>%
  gt::cols_align(align = "left",columns = proxy) %>%
  gt::cols_align(align = "center",columns = c(q05,q90,median)) %>%
  gt::tab_style(
    style = list(gt::cell_fill(color = "#ECECF8"),gt::cell_text(weight = "bold")),
    locations = gt::cells_column_labels()) %>%
  gt::tab_style(
    style = list(gt::cell_fill(color = "#F4F4FA"),gt::cell_text(weight = "bold")),
    locations = gt::cells_row_groups()) %>%
  gt::tab_options(
    table.width = gt::pct(100),
    table.font.size = gt::px(14),
    heading.title.font.size = gt::px(17),
    heading.subtitle.font.size = gt::px(12),
    column_labels.font.weight = "bold",
    data_row.padding = gt::px(7))

print(human_pressure_table_60)
gt::gtsave(human_pressure_table_60,file.path(results_directory_60,"TABLE_1_human_pressure_distribution.html"),inline_css = TRUE)
writeLines(as.character(gt::as_latex(human_pressure_table_60)),file.path(results_directory_60,"TABLE_1_human_pressure_distribution.tex"))


#' ----------------------------------------------------------------------------- STEP 3: select model covariates ----
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
#' (v) TABLE n°2: export the model comparison as an annex in HTML and LaTeX


# Parameters ----
#' process_labels_60, human_pressure_variables_60 and results_directory_60 come
#' from STEP 2 and are not redefined here.
models_per_metric_60 <- 3

model_structure_labels_60 <- c(
  elevation = "Elevation",
  ruggedness = "Ruggedness",
  elevation_ruggedness = "Elevation + ruggedness",
  forest = "Elevation + forest",
  open_habitat = "Elevation + open habitat",
  open_habitat_ruggedness = "Elevation + open habitat + ruggedness",
  ridgeline = "Elevation + distance to ridgeline")

# 3.1 Define environmental adjustment structures ----
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
  hfi_variable = human_pressure_variables_60,
  model_structure = names(model_covariate_structures_60)) %>%
  dplyr::mutate(
    hfi_term = paste0(hfi_variable,"_z"),
    model_id = paste(process,hfi_variable,model_structure,sep = "__"))

# 3.2 Common raw Q05-Q95 contrast ----
#' The lowest Q95 and the highest Q05 across processes are retained so that all
#' models express their effect over one single raw interval.
hfi_thresholds_60 <- dplyr::bind_rows(
  lapply(names(backbone_process_data_60),function(process_name) {
    dplyr::bind_rows(lapply(human_pressure_variables_60,function(hfi_variable) {
      values <- backbone_process_data_60[[process_name]][[hfi_variable]]
      values <- values[is.finite(values)]
      tibble::tibble(
        process = process_name,
        hfi_variable = hfi_variable,
        hfi_q05 = as.numeric(stats::quantile(values,0.05,names = FALSE)),
        hfi_q95 = as.numeric(stats::quantile(values,0.95,names = FALSE)))}))}))

common_hfi_thresholds_60 <- hfi_thresholds_60 %>%
  dplyr::group_by(hfi_variable) %>%
  dplyr::summarise(
    hfi_q05_raw = max(hfi_q05,na.rm = TRUE),
    hfi_q95_raw = min(hfi_q95,na.rm = TRUE),
    .groups = "drop")

#' Recovers the centre and scale applied in Transition_Matrix_Preparation.R.
#' Exporting standardization_parameters from that script would remove the need
#' for this inversion.
derive_hfi_scaling_60 <- function(data,process_name,hfi_variable) {
  raw_value <- data[[hfi_variable]]
  standardized_value <- data[[paste0(hfi_variable,"_z")]]
  scale_value <- stats::cov(raw_value,standardized_value) / stats::var(standardized_value)
  tibble::tibble(process = process_name,hfi_variable = hfi_variable,
                 center = mean(raw_value) - scale_value * mean(standardized_value),
                 scale = scale_value)}

hfi_process_scaling_60 <- dplyr::bind_rows(
  lapply(names(backbone_process_data_60),function(process_name) {
    dplyr::bind_rows(lapply(human_pressure_variables_60,function(hfi_variable) {
      derive_hfi_scaling_60(backbone_process_data_60[[process_name]],process_name,hfi_variable)}))}))

hfi_contrast_scaling_60 <- tidyr::expand_grid(
  process = names(backbone_process_data_60),
  hfi_variable = human_pressure_variables_60) %>%
  dplyr::left_join(common_hfi_thresholds_60,by = "hfi_variable") %>%
  dplyr::left_join(hfi_process_scaling_60,by = c("process","hfi_variable")) %>%
  dplyr::mutate(
    hfi_q05_z = (hfi_q05_raw - center) / scale,
    hfi_q95_z = (hfi_q95_raw - center) / scale,
    hfi_q95_q05_difference_z = hfi_q95_z - hfi_q05_z)

# 3.3 Build process-specific candidate formulas ----
candidate_formulas_60 <- stats::setNames(
  lapply(seq_len(nrow(model_covariate_metadata_60)),function(i) {
    additional_terms <- c(
      model_covariate_metadata_60$hfi_term[[i]],
      model_covariate_structures_60[[model_covariate_metadata_60$model_structure[[i]]]])
    stats::update.formula(
      selected_backbone_formulas_60[[model_covariate_metadata_60$process[[i]]]],
      stats::as.formula(paste(". ~ . +",paste(additional_terms,collapse = " + "))))}),
  model_covariate_metadata_60$model_id)

# 3.4 Fit all weighted binomial GLMMs ----
fit_candidate_model_60 <- function(model_id) {
  model_information <- model_covariate_metadata_60 %>%
    dplyr::filter(.data$model_id == .env$model_id)
  lme4::glmer(
    formula = candidate_formulas_60[[model_id]],
    data = backbone_process_data_60[[model_information$process[[1]]]],
    weights = backbone_weight,
    family = stats::binomial(link = "logit"),
    nAGQ = 1,
    control = control_glmer_60)}

candidate_models_60 <- stats::setNames(
  lapply(model_covariate_metadata_60$model_id,fit_candidate_model_60),
  model_covariate_metadata_60$model_id)

# 3.5 Extract model fit and HFI effects ----
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
      dplyr::filter(.data$process == .env$process_name,.data$hfi_variable == .env$hfi_variable)
    model_loglik <- stats::logLik(fitted_model)
    convergence_message <- fitted_model@optinfo$conv$lme4$messages
    optimizer_code <- fitted_model@optinfo$conv$opt
    hfi_coefficient <- coefficient_table[hfi_term,"Estimate"]
    hfi_standard_error <- coefficient_table[hfi_term,"Std. Error"]
    hfi_confidence_low <- hfi_coefficient - stats::qnorm(0.975) * hfi_standard_error
    hfi_confidence_high <- hfi_coefficient + stats::qnorm(0.975) * hfi_standard_error
    hfi_difference_z <- contrast$hfi_q95_q05_difference_z[[1]]
    tibble::tibble(
      process = process_name,
      hfi_variable = hfi_variable,
      model_structure = model_information$model_structure[[1]],
      n_observations = stats::nobs(fitted_model),
      model_df = attr(model_loglik,"df"),
      AIC = stats::AIC(fitted_model),
      hfi_coefficient = hfi_coefficient,
      hfi_confidence_low = hfi_confidence_low,
      hfi_confidence_high = hfi_confidence_high,
      Q95_Q05_log_odds = hfi_coefficient * hfi_difference_z,
      Q95_Q05_log_odds_CI_low = hfi_confidence_low * hfi_difference_z,
      Q95_Q05_log_odds_CI_high = hfi_confidence_high * hfi_difference_z,
      converged = is.null(convergence_message) && (is.null(optimizer_code) || all(optimizer_code == 0)),
      singular = lme4::isSingular(fitted_model,tol = 1e-4))})) %>%
  dplyr::group_by(process) %>%
  dplyr::mutate(delta_AIC = AIC - min(AIC,na.rm = TRUE)) %>%
  dplyr::ungroup()

#' Convergence is no longer a table column, so the failures are listed here
#' instead. An AIC comparison is only valid if every model in a process was fitted
#' on the same number of observations.
cat("\n== Models that failed to converge or are singular ==\n")
model_convergence_issues_60 <- model_comparison_table_60 %>%
  dplyr::filter(!converged | singular) %>%
  dplyr::transmute(process,hfi_variable,model_structure,converged,singular)
if(nrow(model_convergence_issues_60) == 0L) cat("None.\n") else print(as.data.frame(model_convergence_issues_60),row.names = FALSE)

cat("\n== Sample sizes within each process ==\n")
print(as.data.frame(model_comparison_table_60 %>%
                      dplyr::group_by(process) %>%
                      dplyr::summarise(n_distinct_sample_sizes = dplyr::n_distinct(n_observations),
                                       sample_sizes = paste(sort(unique(n_observations)),collapse = ", "),
                                       .groups = "drop")),row.names = FALSE)


# TABLE 2: model comparison, annex ---------------------------------------------
top_models_60 <- model_comparison_table_60 %>%
  dplyr::group_by(process,hfi_variable) %>%
  dplyr::slice_min(order_by = AIC,n = models_per_metric_60,with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(process) %>%
  dplyr::arrange(AIC,.by_group = TRUE) %>%
  dplyr::mutate(rank = dplyr::row_number(),best_model = rank == 1L) %>%
  dplyr::ungroup()

model_comparison_table_data_60 <- top_models_60 %>%
  dplyr::transmute(
    process_label = unname(process_labels_60[process]),
    best_model,
    rank,
    metric = unname(human_pressure_labels_60[hfi_variable]),
    adjustment = unname(model_structure_labels_60[model_structure]),
    df = model_df,
    AIC = AIC,
    delta_AIC = delta_AIC,
    coefficient = sprintf("%.3f [%.3f, %.3f]",hfi_coefficient,hfi_confidence_low,hfi_confidence_high),
    log_odds = sprintf("%.3f [%.3f, %.3f]",Q95_Q05_log_odds,Q95_Q05_log_odds_CI_low,Q95_Q05_log_odds_CI_high))

model_comparison_gt_60 <- model_comparison_table_data_60 %>%
  dplyr::select(-best_model) %>%
  gt::gt(groupname_col = "process_label") %>%
  gt::cols_label(rank = "Rank",metric = "Proxy of human presence",adjustment = "Environmental covariates",
                 df = "Degrees of freedom",AIC = "AIC",delta_AIC = gt::md("&Delta;AIC"),
                 coefficient = "Coefficient [95% CI]",log_odds = gt::md("Q05&ndash;Q95 log-odds [95% CI]")) %>%
  gt::fmt_number(columns = AIC,decimals = 1) %>%
  gt::fmt_number(columns = delta_AIC,decimals = 2) %>%
  gt::cols_align(align = "left",columns = c(metric,adjustment,coefficient,log_odds)) %>%
  gt::cols_align(align = "center",columns = c(rank,df,AIC,delta_AIC)) %>%
  gt::tab_style(style = list(gt::cell_fill(color = "#ECECF8"),gt::cell_text(weight = "bold")),
                locations = gt::cells_column_labels()) %>%
  gt::tab_style(style = list(gt::cell_fill(color = "#F4F4FA"),gt::cell_text(weight = "bold")),
                locations = gt::cells_row_groups()) %>%
  gt::tab_style(style = gt::cell_text(weight = "bold"),
                locations = gt::cells_body(rows = model_comparison_table_data_60$best_model)) %>%
  gt::tab_options(table.font.size = gt::px(13),heading.title.font.size = gt::px(17),
                  heading.subtitle.font.size = gt::px(12),data_row.padding = gt::px(6),
                  row_group.padding = gt::px(8))

print(model_comparison_gt_60)
gt::gtsave(model_comparison_gt_60,file.path(results_directory_60,"TABLE_2_model_comparison.html"),inline_css = TRUE)
writeLines(as.character(gt::as_latex(model_comparison_gt_60)),file.path(results_directory_60,"TABLE_2_model_comparison.tex"))

#' RESULT OF STEP 3:
#' The proxy of human presence selected is "built cell share 100m"
#' The model formula selected for both ACC dataset (aerial versus terrestrial and
#' feeding versus resting) and the GPS dataset (aerial versus terrestrial) include 
#' elevation, open habitat and ruggedness.   


# ----------------------------------------------------------------------------- STEP 4: sensitivity to environmental adjustment ----
#' **Steps:**
#' (i) for both binary state GPS and ACC datasets and for the ACC feeding versus resting
#'  fit the full elevation + open habitat + ruggedness structure and remove each
#'     environmental term separately;
#' (ii) quantify changes in AIC, HFI coefficients and common Q05-Q95
#'       log-odds contrasts relative to the corresponding full model;
#' (iii) TABLE n°3: export one sensitivity table as an annex in HTML and LaTeX


# Parameters ----
hfi_sensitivity_variable_60 <- "settlement_density_raw"
hfi_sensitivity_term_60 <- paste0(hfi_sensitivity_variable_60,"_z")

sensitivity_environmental_terms_60 <- c(
  elevation = "elevation_100m_z",
  open_habitat = "prop_low_vegetation_5cells_z",
  ruggedness = "ruggedness_100m_z")

environmental_term_labels_60 <- c(
  elevation = "Elevation",
  open_habitat = "Open habitat",
  ruggedness = "Ruggedness")

# 4.1 Define the full and reduced structures ----
sensitivity_structures_60 <- c(
  list(full = unname(sensitivity_environmental_terms_60)),
  stats::setNames(
    lapply(names(sensitivity_environmental_terms_60),function(term_name) {
      unname(sensitivity_environmental_terms_60[names(sensitivity_environmental_terms_60) != term_name])}),
    paste0("without_",names(sensitivity_environmental_terms_60))))

sensitivity_model_metadata_60 <- tidyr::expand_grid(
  process = names(backbone_process_data_60),
  model_structure = names(sensitivity_structures_60)) %>%
  dplyr::mutate(
    structure_order = match(model_structure,names(sensitivity_structures_60)),
    removed_variable = dplyr::if_else(
      model_structure == "full","None - full model",
      unname(environmental_term_labels_60[sub("^without_","",model_structure)])),
    model_id = paste(process,model_structure,sep = "__"))

# 4.2 Build the sensitivity formulas ----
sensitivity_formulas_60 <- stats::setNames(
  lapply(seq_len(nrow(sensitivity_model_metadata_60)),function(i) {
    additional_terms <- c(hfi_sensitivity_term_60,
                          sensitivity_structures_60[[sensitivity_model_metadata_60$model_structure[[i]]]])
    stats::update.formula(
      selected_backbone_formulas_60[[sensitivity_model_metadata_60$process[[i]]]],
      stats::as.formula(paste(". ~ . +",paste(additional_terms,collapse = " + "))))}),
  sensitivity_model_metadata_60$model_id)

# 4.3 Fit the weighted binomial GLMMs ----
sensitivity_models_60 <- stats::setNames(
  lapply(seq_len(nrow(sensitivity_model_metadata_60)),function(i) {
    lme4::glmer(
      formula = sensitivity_formulas_60[[i]],
      data = backbone_process_data_60[[sensitivity_model_metadata_60$process[[i]]]],
      weights = backbone_weight,
      family = stats::binomial(link = "logit"),
      nAGQ = 1,
      control = control_glmer_60)}),
  sensitivity_model_metadata_60$model_id)

# 4.4 Extract AIC, HFI effects and common Q05-Q95 contrasts ----
sensitivity_model_results_60 <- dplyr::bind_rows(
  lapply(seq_len(nrow(sensitivity_model_metadata_60)),function(i) {
    fitted_model <- sensitivity_models_60[[i]]
    process_name <- sensitivity_model_metadata_60$process[[i]]
    coefficient_table <- summary(fitted_model)$coefficients
    convergence_message <- fitted_model@optinfo$conv$lme4$messages
    optimizer_code <- fitted_model@optinfo$conv$opt
    hfi_difference_z <- hfi_contrast_scaling_60$hfi_q95_q05_difference_z[
      hfi_contrast_scaling_60$process == process_name &
        hfi_contrast_scaling_60$hfi_variable == hfi_sensitivity_variable_60][[1]]
    hfi_coefficient <- coefficient_table[hfi_sensitivity_term_60,"Estimate"]
    hfi_standard_error <- coefficient_table[hfi_sensitivity_term_60,"Std. Error"]
    hfi_confidence_low <- hfi_coefficient - stats::qnorm(0.975) * hfi_standard_error
    hfi_confidence_high <- hfi_coefficient + stats::qnorm(0.975) * hfi_standard_error
    tibble::tibble(
      process = process_name,
      model_structure = sensitivity_model_metadata_60$model_structure[[i]],
      structure_order = sensitivity_model_metadata_60$structure_order[[i]],
      removed_variable = sensitivity_model_metadata_60$removed_variable[[i]],
      AIC = stats::AIC(fitted_model),
      hfi_coefficient = hfi_coefficient,
      hfi_confidence_low = hfi_confidence_low,
      hfi_confidence_high = hfi_confidence_high,
      Q95_Q05_log_odds = hfi_coefficient * hfi_difference_z,
      Q95_Q05_log_odds_CI_low = hfi_confidence_low * hfi_difference_z,
      Q95_Q05_log_odds_CI_high = hfi_confidence_high * hfi_difference_z,
      converged = is.null(convergence_message) && (is.null(optimizer_code) || all(optimizer_code == 0)),
      singular = lme4::isSingular(fitted_model,tol = 1e-4))}))

# 4.5 Calculate changes relative to each full model ----
hfi_coefficient_sensitivity_60 <- sensitivity_model_results_60 %>%
  dplyr::group_by(process) %>%
  dplyr::mutate(
    delta_AIC_from_full = AIC - AIC[model_structure == "full"][[1]],
    hfi_coefficient_change = hfi_coefficient - hfi_coefficient[model_structure == "full"][[1]],
    Q95_Q05_log_odds_change = Q95_Q05_log_odds - Q95_Q05_log_odds[model_structure == "full"][[1]],
    hfi_sign_changed = model_structure != "full" &
      sign(hfi_coefficient) != sign(hfi_coefficient[model_structure == "full"][[1]])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(process_label = unname(process_labels_60[process])) %>%
  dplyr::arrange(process_label,structure_order)

#' Convergence is not a table column, so the failures are listed here instead.
cat("\n== Sensitivity models that failed to converge or are singular ==\n")
sensitivity_convergence_issues_60 <- hfi_coefficient_sensitivity_60 %>%
  dplyr::filter(!converged | singular) %>%
  dplyr::transmute(process,model_structure,converged,singular)
if(nrow(sensitivity_convergence_issues_60) == 0L) cat("None.\n") else print(as.data.frame(sensitivity_convergence_issues_60),row.names = FALSE)


# TABLE 3: environmental sensitivity, annex ------------------------------------
sensitivity_table_data_60 <- hfi_coefficient_sensitivity_60 %>%
  dplyr::transmute(
    process_label,
    full_model = model_structure == "full",
    removed_variable,
    AIC,
    delta_AIC_from_full,
    coefficient = sprintf("%.4f [%.4f, %.4f]",hfi_coefficient,hfi_confidence_low,hfi_confidence_high),
    coefficient_change = sprintf("%+.4f",hfi_coefficient_change),
    log_odds = sprintf("%.4f [%.4f, %.4f]",Q95_Q05_log_odds,Q95_Q05_log_odds_CI_low,Q95_Q05_log_odds_CI_high),
    log_odds_change = sprintf("%+.4f",Q95_Q05_log_odds_change))

sensitivity_gt_60 <- sensitivity_table_data_60 %>%
  dplyr::select(-full_model) %>%
  gt::gt(groupname_col = "process_label") %>%
  gt::cols_label(removed_variable = "Variable removed",AIC = "AIC",
                 delta_AIC_from_full = gt::md("&Delta;AIC vs full"),
                 coefficient = "Coefficient [95% CI]",coefficient_change = "Change in coefficient",
                 log_odds = gt::md("Q05&ndash;Q95 log-odds [95% CI]"),log_odds_change = "Change in log-odds") %>%
  gt::fmt_number(columns = AIC,decimals = 1) %>%
  gt::fmt_number(columns = delta_AIC_from_full,decimals = 2,force_sign = TRUE) %>%
  gt::cols_align(align = "left",columns = c(removed_variable,coefficient,log_odds)) %>%
  gt::cols_align(align = "center",columns = c(AIC,delta_AIC_from_full,coefficient_change,log_odds_change)) %>%
  gt::tab_style(style = list(gt::cell_fill(color = "#ECECF8"),gt::cell_text(weight = "bold")),
                locations = gt::cells_column_labels()) %>%
  gt::tab_style(style = list(gt::cell_fill(color = "#F4F4FA"),gt::cell_text(weight = "bold")),
                locations = gt::cells_row_groups()) %>%
  gt::tab_style(style = gt::cell_text(weight = "bold"),
                locations = gt::cells_body(rows = sensitivity_table_data_60$full_model)) %>%
  gt::tab_options(table.font.size = gt::px(13),heading.title.font.size = gt::px(17),
                  heading.subtitle.font.size = gt::px(12),data_row.padding = gt::px(6),
                  row_group.padding = gt::px(8))

print(sensitivity_gt_60)
gt::gtsave(sensitivity_gt_60,file.path(results_directory_60,"TABLE_3_environmental_sensitivity.html"),inline_css = TRUE)
writeLines(as.character(gt::as_latex(sensitivity_gt_60)),file.path(results_directory_60,"TABLE_3_environmental_sensitivity.tex"))

# RESULTS OF STEP 4 ----
#' We select the variables that influence the coefficient of settlement density:
#' GPS aerial vs terrestrial: HFI + elevation + ruggedness
#' ACC aerial vs terrestrial: HFI + elevation + ruggedness
#' ACC feeding vs resting: HFI + elevation


# ----------------------------------------------------------------------------- STEP 5: compile the probabilities of transitions ----
#' **Note:** we keep all covariate at their mean (0 because these have been 
#' standardised) and the random intercept at 0. 
#'
#' The ACC classification has three states but only two models, so its
#' probabilities are recomposed by the chain rule:
#'   P(feeding) = (1 - P(aerial)) * P(feeding | terrestrial)
#'   P(resting) = (1 - P(aerial)) * (1 - P(feeding | terrestrial))
#'
#' **Steps:**
#' (i) fit the models retained at the end of STEP 4;
#' (ii) build the reference eagle and its prediction matrices at Q05 and Q95;
#' (iii) compute each transition probability and recompose the ACC states;
#' (iv) convert to Q95-Q05 log-odds with a parametric bootstrap interval;
#' (v) TABLE n°4: export the transition probabilities as HTML and LaTeX


# Parameters ----
final_hfi_variable_60 <- "settlement_density_raw"
final_hfi_term_60 <- paste0(final_hfi_variable_60,"_z")

final_environmental_terms_60 <- list(
  GPS_aerial_vs_terrestrial = c("elevation_100m_z","ruggedness_100m_z"),
  ACC_aerial_vs_terrestrial = c("elevation_100m_z","ruggedness_100m_z"),
  ACC_feeding_vs_resting = "elevation_100m_z")

transition_labels_60 <- c(
  GPS_aerial = "P(A→A)",
  GPS_terrestrial = "P(A→T)",
  ACC_aerial = "P(A→A)",
  ACC_feeding = "P(A→F)",
  ACC_resting = "P(A→R)")

transition_classification_60 <- c(GPS_aerial = "GPS",GPS_terrestrial = "GPS",
                                  ACC_aerial = "ACC",ACC_feeding = "ACC",ACC_resting = "ACC")

n_draws_60 <- 2000

# 5.1 Define and fit the retained models ----
final_model_formulas_60 <- stats::setNames(
  lapply(names(final_environmental_terms_60),function(process_name) {
    stats::update.formula(
      selected_backbone_formulas_60[[process_name]],
      stats::as.formula(paste(". ~ . +",paste(c(final_hfi_term_60,final_environmental_terms_60[[process_name]]),collapse = " + "))))}),
  names(final_environmental_terms_60))

final_models_60 <- stats::setNames(
  lapply(names(final_model_formulas_60),function(process_name) {
    lme4::glmer(
      formula = final_model_formulas_60[[process_name]],
      data = backbone_process_data_60[[process_name]],
      weights = backbone_weight,
      family = stats::binomial(link = "logit"),
      nAGQ = 1,
      control = control_glmer_60)}),
  names(final_model_formulas_60))

cat("\n== Retained models that failed to converge or are singular ==\n")
final_convergence_issues_60 <- dplyr::bind_rows(lapply(names(final_models_60),function(process_name) {
  fitted_model <- final_models_60[[process_name]]
  convergence_message <- fitted_model@optinfo$conv$lme4$messages
  optimizer_code <- fitted_model@optinfo$conv$opt
  tibble::tibble(
    process = process_name,
    n_observations = stats::nobs(fitted_model),
    converged = is.null(convergence_message) && (is.null(optimizer_code) || all(optimizer_code == 0)),
    singular = lme4::isSingular(fitted_model,tol = 1e-4))})) %>%
  dplyr::filter(!converged | singular)
if(nrow(final_convergence_issues_60) == 0L) cat("None.\n") else print(as.data.frame(final_convergence_issues_60),row.names = FALSE)

# 5.2 Reference eagle and prediction matrices at the common raw Q05 and Q95 ----
#' Both ACC models are evaluated on the same reference eagle, built from the
#' aerial-origin ACC data, so that the three states are recomposed consistently.
reference_data_60 <- list(
  GPS_aerial_vs_terrestrial = backbone_process_data_60[["GPS_aerial_vs_terrestrial"]],
  ACC_aerial_vs_terrestrial = backbone_process_data_60[["ACC_aerial_vs_terrestrial"]],
  ACC_feeding_vs_resting = backbone_process_data_60[["ACC_aerial_vs_terrestrial"]])

hfi_levels_60 <- stats::setNames(
  lapply(names(final_models_60),function(process_name) {
    contrast <- hfi_contrast_scaling_60[hfi_contrast_scaling_60$process == process_name &
                                          hfi_contrast_scaling_60$hfi_variable == final_hfi_variable_60,]
    c(Q05 = contrast$hfi_q05_z[[1]],Q95 = contrast$hfi_q95_z[[1]])}),
  names(final_models_60))

build_reference_matrix_60 <- function(model,data,hfi_value_z) {
  fixed_terms <- stats::delete.response(stats::terms(lme4::nobars(stats::formula(model))))
  profile <- as.data.frame(data)[1,all.vars(fixed_terms),drop = FALSE]
  for(variable in all.vars(fixed_terms)) {
    if(is.numeric(data[[variable]])) profile[[variable]] <- mean(data[[variable]],na.rm = TRUE) else {
      values <- factor(data[[variable]])
      profile[[variable]] <- factor(names(which.max(table(values))),levels = levels(values))}}
  profile[[final_hfi_term_60]] <- hfi_value_z
  stats::model.matrix(fixed_terms,data = profile)}

reference_matrices_60 <- stats::setNames(
  lapply(names(final_models_60),function(process_name) {
    stats::setNames(lapply(c("Q05","Q95"),function(level) {
      build_reference_matrix_60(final_models_60[[process_name]],reference_data_60[[process_name]],
                                hfi_levels_60[[process_name]][[level]])}),c("Q05","Q95"))}),
  names(final_models_60))

# can be deleted if nothing happens: fails if model.matrix reorders the columns
stopifnot(all(vapply(names(final_models_60),function(process_name) {
  identical(colnames(reference_matrices_60[[process_name]]$Q05),names(lme4::fixef(final_models_60[[process_name]])))},logical(1))))

# 5.3 Transition probabilities of the reference eagle ----
transition_probabilities_60 <- function(coefficients,level) {
  p_gps <- stats::plogis(as.vector(reference_matrices_60$GPS_aerial_vs_terrestrial[[level]] %*% coefficients$GPS_aerial_vs_terrestrial))
  p_acc <- stats::plogis(as.vector(reference_matrices_60$ACC_aerial_vs_terrestrial[[level]] %*% coefficients$ACC_aerial_vs_terrestrial))
  p_feeding_given <- stats::plogis(as.vector(reference_matrices_60$ACC_feeding_vs_resting[[level]] %*% coefficients$ACC_feeding_vs_resting))
  c(GPS_aerial = p_gps,
    GPS_terrestrial = 1 - p_gps,
    ACC_aerial = p_acc,
    ACC_feeding = (1 - p_acc) * p_feeding_given,
    ACC_resting = (1 - p_acc) * (1 - p_feeding_given))}

final_coefficients_60 <- lapply(final_models_60,lme4::fixef)
probability_q05_60 <- transition_probabilities_60(final_coefficients_60,"Q05")
probability_q95_60 <- transition_probabilities_60(final_coefficients_60,"Q95")
log_odds_estimate_60 <- stats::qlogis(probability_q95_60) - stats::qlogis(probability_q05_60)

# 5.4 Parametric bootstrap interval on the Q95-Q05 log-odds ----
#' The composed ACC probabilities have no closed-form standard error, so the
#' fixed effects are redrawn from their sampling distribution and the
#' probabilities of the reference eagle are recomputed. For the three
#' single-model transitions this interval matches the Wald interval of the
#' coefficient times the Q95-Q05 contrast. Drawing the two ACC models independently is
#' legitimate: the three-state likelihood factorises exactly into the aerial
#' versus terrestrial binomial and the feeding versus resting binomial, so their
#' fixed-effect blocks are asymptotically independent. The draws do assume that
#' the two individual random intercepts are independent, which is already
#' assumed by fitting the two models separately.
coefficient_cholesky_60 <- lapply(final_models_60,function(model) t(chol(as.matrix(stats::vcov(model)))))

set.seed(20261125)
bootstrap_log_odds_60 <- vapply(seq_len(n_draws_60),function(draw) {
  coefficients <- stats::setNames(lapply(names(final_models_60),function(process_name) {
    as.vector(final_coefficients_60[[process_name]] +
                coefficient_cholesky_60[[process_name]] %*% stats::rnorm(length(final_coefficients_60[[process_name]])))}),
    names(final_models_60))
  stats::qlogis(transition_probabilities_60(coefficients,"Q95")) -
    stats::qlogis(transition_probabilities_60(coefficients,"Q05"))},numeric(5))

log_odds_confidence_60 <- apply(bootstrap_log_odds_60,1,stats::quantile,probs = c(0.025,0.975))


# TABLE 4: transition probabilities --------------------------------------------
final_common_hfi_thresholds_60 <- common_hfi_thresholds_60[common_hfi_thresholds_60$hfi_variable == final_hfi_variable_60,]

transition_probability_table_data_60 <- tibble::tibble(
  classification = unname(transition_classification_60[names(probability_q05_60)]),
  transition = unname(transition_labels_60[names(probability_q05_60)]),
  probability_q05 = unname(probability_q05_60),
  probability_q95 = unname(probability_q95_60),
  log_odds = sprintf("%+.3f (%.3f – %.3f)",log_odds_estimate_60,
                     log_odds_confidence_60[1,],log_odds_confidence_60[2,]))

transition_probability_gt_60 <- transition_probability_table_data_60 %>%
  gt::gt(groupname_col = "classification") %>%
  gt::tab_header(
    title = gt::md(paste0("**Conditional transition probabilities at both ends of the ",tolower(human_pressure_labels_60[[final_hfi_variable_60]])," gradient**")),
    subtitle = sprintf("Common raw contrast: Q05 = %.3f, Q95 = %.3f. A = aerial, T = terrestrial, F = feeding, R = resting.",
                       final_common_hfi_thresholds_60$hfi_q05_raw,final_common_hfi_thresholds_60$hfi_q95_raw)) %>%
  gt::cols_label(transition = "Type of transition",
                 probability_q05 = gt::md("5th quantile<br>(Q05)"),
                 probability_q95 = gt::md("95th quantile<br>(Q95)"),
                 log_odds = gt::md("Q95 &minus; Q05<br>log-odds")) %>%
  gt::fmt_number(columns = c(probability_q05,probability_q95),decimals = 3) %>%
  gt::cols_align(align = "left",columns = transition) %>%
  gt::cols_align(align = "center",columns = c(probability_q05,probability_q95,log_odds)) %>%
  gt::tab_style(style = list(gt::cell_fill(color = "#ECECF8"),gt::cell_text(weight = "bold")),
                locations = gt::cells_column_labels()) %>%
  gt::tab_style(style = list(gt::cell_fill(color = "#F4F4FA"),gt::cell_text(weight = "bold")),
                locations = gt::cells_row_groups()) %>%
  gt::tab_options(table.font.size = gt::px(14),heading.title.font.size = gt::px(17),
                  heading.subtitle.font.size = gt::px(12),data_row.padding = gt::px(7),
                  row_group.padding = gt::px(8)) %>%
  gt::tab_source_note(source_note = gt::md(paste0(
    "Probabilities are predicted for a reference eagle (random intercept at 0, numeric covariates at their mean, ",
    "factors at their most frequent level); only the human-footprint proxy changes between the two quantiles. ",
    "Log-odds are the difference between the two quantiles, with a 2.5-97.5 % parametric bootstrap interval over ",
    n_draws_60," draws of the fixed effects. A positive value means the transition becomes more likely as human pressure increases, all else being equal. ",
    "The two ACC states conditional on landing are obtained from a sequential decomposition whose likelihood factorises, ",
    "so their two models are drawn independently in the bootstrap.")))

print(transition_probability_gt_60)
gt::gtsave(transition_probability_gt_60,file.path(results_directory_60,"TABLE_4_transition_probabilities.html"),inline_css = TRUE)
writeLines(as.character(gt::as_latex(transition_probability_gt_60)),file.path(results_directory_60,"TABLE_4_transition_probabilities.tex"))
