#' -----------------------------------------------------------------------------
# Title: Model Validation on GPS-bassed classification of behaviors ----
#' Authors : Louise Faure
#' Date : 24.07.26
#' 
#' Info : this script follow the Covariates&HFI_Selection.R script where covariates
#' and one model representation is selected.  
#' 
# Main steps:
#' (1) fit two open habitat models Open Habitat models 
#'     (i) the first model assume individuals have the same response to HFI 
#'     (ii) the second model assume varying slope per individuals to HFI
#'     (iii) resume results : (a) a table for model comparison (stability of HFI, 
#'     AIC and EDF) and (b) retain the names of the individuals with a negative 
#'     response to HFI.
#'
#' (2) Model controls on autocorrelation within pseudo residuals
#'    (i) generate DHARMa pseudo-residuals for the common-slope and individual-
#'    random-slope models;
#'    (ii) test temporal autocorrelation within each individual and burst at lags 
#'    of 60, 120 and 180 minutes;
#'    (iii) summarise temporal residual correlations across individuals using 
#'    bootstrap 95% confidence intervals;
#'    (iv) aggregate residuals within 1-km spatial cells and calculate Moran’s I 
#'    separately for each individual;
#'    (v) correct spatial test p-values for multiple comparisons using the 
#'    Benjamini–Hochberg FDR procedure;
#'
#' (3) Confirmation of the results: bootstraping on individuals to confirm the 
#' confidence interval for HFI



# library
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(ggplot2)
library(mgcv)


# data
GE_standardize_weighted <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/hfi_model_data_60_standardized_weighted.rds")

# parameters
gam_selection_method <- "ML"
gam_final_method <- "REML"
hfi_variable_60 <- "hfi_mean_1000m_z"
confidence_level_60 <- 0.95
confidence_critical_value_60 <- stats::qnorm(1 - (1 - confidence_level_60) / 2)

# output directory
results_directory_60 <- paste0(
  "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/",
  "THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results")

dir.create(results_directory_60,recursive = TRUE,showWarnings = FALSE)


#------------------------------------------------------------------------------- STEP 1: fit common- and individual-HFI-slope open-habitat models ----
#' **Steps:**
#' (i) verify the modelling dataset and individual weights;
#' (ii) decompose HFI into within- and between-individual components;
#' (iii) fit a common-slope and an individual-random-slope model;
#' (iv) compare HFI stability, confidence intervals, AIC and EDF;
#' (v) extract individual HFI slopes and individual-specific Q05-Q95 contrasts;
#' (vi) retain individuals with negative HFI-response estimates.

# 1.1 Prepare the common modelling dataset ----
# hfi_mean_1000m_z = hfi standardized.
# hfi_between_z = mean globally standardized HFI used by each individual.
# hfi_within_z = deviation of each observation from the individual's mean HFI.
open_habitat_model_data_60 <-
  GE_standardize_weighted %>%
  dplyr::mutate(
    individual_id = droplevels(factor(individual_id))) %>%
  dplyr::group_by(individual_id) %>%
  dplyr::mutate(hfi_between_z =mean(.data[[hfi_variable_60]]),
    hfi_within_z =.data[[hfi_variable_60]] - hfi_between_z) %>%
  dplyr::ungroup()

# 1.2 Define the common-slope model ----
# hfi_within_z estimates the response to HFI variation encountered within
# individuals.
# hfi_between_z controls differences in average exposure among individuals.
common_hfi_slope_formula_60 <-
  remain_aerial ~
  cos_diel +
  sin_time +
  age_z +
  duration_z +
  duration_z2 +
  hfi_within_z +
  hfi_between_z +
  dem_elevation_z +
  prop_low_vegetation_5cells_z +
  s(
    individual_id,
    bs = "re")


# 1.3 Define the individual-random-slope model ----
# The additional random-effect term estimates individual deviations from the
# population within-individual HFI slope. The model answer the following question :
# "Do individuals differ in their response to a increase of HFI measured on a common scale?
individual_hfi_slope_formula_60 <-
  remain_aerial ~
  cos_diel +
  sin_time +
  age_z +
  duration_z +
  duration_z2 +
  hfi_within_z +
  hfi_between_z +
  dem_elevation_z +
  prop_low_vegetation_5cells_z +
  s(
    individual_id,
    bs = "re"
  ) +
  s(
    hfi_within_z,
    individual_id,
    bs = "re"
  )


# 1.4 Fit both models using ML ----
# ML is retained because the models are compared using AIC.
# Both models use exactly the same observations and individual weights.
common_hfi_slope_model_ml_60 <-
  mgcv::gam(
    formula =
      common_hfi_slope_formula_60,
    
    data =
      open_habitat_model_data_60,
    
    weights =
      individual_weight,
    
    family =
      stats::binomial(
        link = "logit"
      ),
    
    method =
      gam_selection_method
  )

individual_hfi_slope_model_ml_60 <-
  mgcv::gam(
    formula =
      individual_hfi_slope_formula_60,
    
    data =
      open_habitat_model_data_60,
    
    weights =
      individual_weight,
    
    family =
      stats::binomial(
        link = "logit"
      ),
    
    method =
      gam_selection_method
  )

# 1.5 Compact comparison of common- and random-slope models ----
models_to_compare_60 <- list(
  `Common HFI slope` =
    common_hfi_slope_model_ml_60,
  
  `Individual random HFI slopes` =
    individual_hfi_slope_model_ml_60)


# Extract the requested statistics from one model
extract_model_statistics_60 <- function(model) {
  
  model_summary <- summary(model,re.test = TRUE)
  
  coefficient_table <- model_summary$p.table
  
  hfi_estimate <- coefficient_table["hfi_within_z","Estimate"]
  
  hfi_se <-coefficient_table["hfi_within_z","Std. Error"]
  
  hfi_ci_low <-hfi_estimate - 1.96 * hfi_se
  
  hfi_ci_high <- hfi_estimate +1.96 * hfi_se
  
  # Identify the individual random-slope term
  random_slope_term <- grep("hfi_within_z.*individual_id|individual_id.*hfi_within_z",
      rownames(model_summary$s.table),value = TRUE)
  
  random_slope_edf <- if (length(random_slope_term) == 1L) {
      model_summary$s.table[random_slope_term,"edf"]} else {NA_real_}
  
  random_slope_p <-if (length(random_slope_term) == 1L) {
      model_summary$s.table[random_slope_term,ncol(model_summary$s.table)]} else {NA_real_}
  
  tibble::tibble(Statistic = c("Number of observations","Number of individuals","AIC","Total EDF","Within-HFI coefficient",
      "Within-HFI SE","Within-HFI 95% CI","Random-slope EDF","Random-slope p-value","Deviance explained (%)","Converged"),
    
    Value = c(nrow(stats::model.frame(model)),
      
      dplyr::n_distinct(stats::model.frame(model)$individual_id),
      
      sprintf("%.1f",stats::AIC(model)),
      
      sprintf("%.2f",sum(model$edf)),
      
      sprintf("%.3f",hfi_estimate),
      
      sprintf("%.3f",hfi_se),
      
      sprintf("[%.3f, %.3f]",hfi_ci_low,hfi_ci_high),
      
      ifelse(is.na(random_slope_edf),"--",sprintf("%.2f", random_slope_edf)),
      
      ifelse(is.na(random_slope_p),"--",format.pval(random_slope_p,digits = 3,eps = 0.001)),
      
      sprintf("%.2f",100 * model_summary$dev.expl),
      
      as.character(model$converged)))}


# Extract the statistics and place models in columns
open_habitat_model_comparison_table_60 <-
  dplyr::bind_rows(
    lapply(models_to_compare_60,extract_model_statistics_60),
    .id = "Model") %>%
  tidyr::pivot_wider(
    names_from = Model,
    values_from = Value)

print(open_habitat_model_comparison_table_60,n = Inf)
# Statistic              `Common HFI slope` `Individual random HFI slopes`
# Number of observations 10872              10872                         
# Number of individuals  62                 62                            
# AIC                    15707.0            15693.5                       
# Total EDF              49.59              70.70                         
# Within-HFI coefficient 0.053              0.061                         
# Within-HFI SE          0.028              0.033                         
# Within-HFI 95% CI      [-0.002, 0.107]    [-0.003, 0.125]               
# Random-slope EDF       --                 20.98                         
# Random-slope p-value   --                 0.00366                       
# Deviance explained (%) 3.73               4.11                          
# Converged              TRUE               TRUE   

# 1.6 Extract and rank individual HFI slopes ----
# HFI_slope = population HFI slope + individual random-slope deviation.
individual_slope_template_60 <- stats::model.frame(individual_hfi_slope_model_ml_60) %>%
  tibble::as_tibble() %>%
  dplyr::group_by(individual_id) %>%
  dplyr::mutate(n_transitions =dplyr::n()) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# Model matrices at hfi_within_z = 0 and hfi_within_z = 1
slope_matrix_0_60 <- stats::predict(individual_hfi_slope_model_ml_60,
    newdata = individual_slope_template_60 %>%
      dplyr::mutate(hfi_within_z = 0),
    type = "lpmatrix")

slope_matrix_1_60 <- stats::predict(individual_hfi_slope_model_ml_60,
    newdata =individual_slope_template_60 %>%
      dplyr::mutate(hfi_within_z = 1),
    type = "lpmatrix")

# Individual HFI slopes, ranked from most negative to most positive
individual_hfi_slopes_60 <-
  individual_slope_template_60 %>%
  dplyr::transmute(
    individual_id,
    n_transitions,
    HFI_slope = as.numeric((slope_matrix_1_60 - slope_matrix_0_60) %*%
          stats::coef(individual_hfi_slope_model_ml_60))) %>%
  dplyr::arrange(HFI_slope)

print(individual_hfi_slopes_60,n = Inf)





#------------------------------------------------------------------------------- STEP 2: residual temporal and spatial autocorrelation diagnostics ----
#' **Steps:**
#' (i) generate DHARMa pseudo-residuals for the common-slope and individual-random-slope models;
#'  (ii) test temporal autocorrelation within each individual and burst at lags of 60, 120 and 180 minutes;
#'  (iii) summarise temporal residual correlations across individuals using bootstrap 95% confidence intervals;
#'  (iv) aggregate residuals within 1-km spatial cells and calculate Moran’s I separately for each individual;
#'  (v) correct spatial test p-values for multiple comparisons using the Benjamini–Hochberg FDR procedure;

library(DHARMa)
library(spdep)

# Parameters
n_dharma_simulations_60 <- 1000L
temporal_lags_60 <- 1:3
nominal_interval_min_60 <- 60
temporal_tolerance_min_60 <- 15
n_bootstrap_60 <- 999L

spatial_grid_size_m_60 <- 1000
spatial_k_neighbours_60 <- 4L
n_spatial_permutations_60 <- 999L


# Models to diagnose
diagnostic_models_60 <- list(
  `Common HFI slope` =
    common_hfi_slope_model_ml_60,
  
  `Individual HFI slopes` =
    individual_hfi_slope_model_ml_60)

# 2.1 Generate DHARMa pseudo-residuals ----
set.seed(20260724)

# Simulate conditional DHARMa residuals for weighted binary GAMs ----
simulate_weighted_binary_DHARMa_60 <- function(model,n_simulations = 1000L,seed = 1L) {
  
  model_frame <-
    stats::model.frame(model)
  
  observed_response <-
    stats::model.response(model_frame)
  
  fitted_probability <-
    as.numeric(
      stats::fitted(model)
    )
  
  # CONTROL: response and fitted probabilities
  stopifnot(
    length(observed_response) ==
      length(fitted_probability),
    
    all(
      observed_response %in% c(0, 1)
    ),
    
    all(
      is.finite(fitted_probability)
    ),
    
    all(
      fitted_probability >= 0 &
        fitted_probability <= 1
    )
  )
  
  set.seed(seed)
  
  simulated_response <-
    vapply(
      seq_len(n_simulations),
      function(simulation_number) {
        stats::rbinom(
          n =length(fitted_probability),
          
          size =1L,
          
          prob =fitted_probability
        )
      },
      numeric(
        length(fitted_probability)
      )
    )
  
  DHARMa::createDHARMa(
    simulatedResponse =
      simulated_response,
    
    observedResponse =
      observed_response,
    
    fittedPredictedResponse =
      fitted_probability,
    
    integerResponse =
      TRUE)}

dharma_residuals_60 <-
  stats::setNames(
    lapply(
      seq_along(
        diagnostic_models_60
      ),
      function(i) {
        simulate_weighted_binary_DHARMa_60(
          model =
            diagnostic_models_60[[i]],
          
          n_simulations =
            n_dharma_simulations_60,
          
          seed =
            20260724L + i)}),
    names(diagnostic_models_60))

# Attach residuals to the model data
residual_data_60 <- dplyr::bind_rows(
  lapply(
    names(dharma_residuals_60),
    function(model_name) {
      
      residual_object <-
        dharma_residuals_60[[model_name]]
      
      stopifnot(
        length(residual_object$scaledResiduals) ==
          nrow(open_habitat_model_data_60)
      )
      
      open_habitat_model_data_60 %>%
        dplyr::transmute(
          model = model_name,
          individual_id,
          burst_id,
          
          timestamp =
            as.POSIXct(timestamp, tz = "UTC"),
          
          x_3035,
          y_3035,
          
          # Normal-score transformation of DHARMa residuals
          residual_normal =
            stats::qnorm(
              pmin(
                pmax(residual_object$scaledResiduals,1e-6),1 - 1e-6)))
    }
  )
)


# 2.2 Temporal autocorrelation ----
temporal_pairs_60 <- dplyr::bind_rows(
  lapply(temporal_lags_60, function(lag_value) {
      
      residual_data_60 %>%
        dplyr::arrange(model,individual_id,burst_id,timestamp) %>%
        dplyr::group_by(model,
          individual_id,
          burst_id
        ) %>%
        dplyr::mutate(
          residual_previous =
            dplyr::lag(
              residual_normal,
              lag_value
            ),
          
          timestamp_previous =
            dplyr::lag(
              timestamp,
              lag_value
            ),
          
          time_lag_min =
            as.numeric(
              difftime(
                timestamp,
                timestamp_previous,
                units = "mins"
              )
            )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::filter(
          !is.na(residual_previous),
          
          abs(
            time_lag_min -
              lag_value *
              nominal_interval_min_60
          ) <=
            temporal_tolerance_min_60 *
            sqrt(lag_value)
        ) %>%
        dplyr::mutate(
          lag = lag_value
        )
    }
  )
)


# Correlation calculated separately for each individual
temporal_correlations_by_individual_60 <-
  temporal_pairs_60 %>%
  dplyr::group_by(
    model,
    lag,
    individual_id
  ) %>%
  dplyr::summarise(
    n_pairs =
      dplyr::n(),
    
    residual_correlation =
      if (
        dplyr::n() >= 5L &&
        stats::sd(residual_previous) > 0 &&
        stats::sd(residual_normal) > 0
      ) {
        stats::cor(
          residual_previous,
          residual_normal
        )
      } else {
        NA_real_
      },
    
    .groups = "drop"
  ) %>%
  dplyr::filter(
    is.finite(residual_correlation)
  )


# Equal-individual bootstrap using Fisher-transformed correlations
bootstrap_temporal_correlation_60 <- function(
    correlations,
    n_bootstrap,
    seed
) {
  
  correlations <- pmin(
    pmax(correlations, -0.999999),
    0.999999
  )
  
  fisher_z <- atanh(correlations)
  
  set.seed(seed)
  
  bootstrap_values <- replicate(
    n_bootstrap,
    tanh(
      mean(
        sample(
          fisher_z,
          replace = TRUE
        )
      )
    )
  )
  
  c(
    estimate =
      tanh(mean(fisher_z)),
    
    confidence_low =
      stats::quantile(
        bootstrap_values,
        0.025,
        names = FALSE
      ),
    
    confidence_high =
      stats::quantile(
        bootstrap_values,
        0.975,
        names = FALSE
      )
  )
}


temporal_autocorrelation_diagnostics_60 <-
  temporal_correlations_by_individual_60 %>%
  dplyr::group_by(
    model,
    lag
  ) %>%
  dplyr::group_modify(
    ~ {
      result <- bootstrap_temporal_correlation_60(
        correlations =
          .x$residual_correlation,
        
        n_bootstrap =
          n_bootstrap_60,
        
        seed =
          20260724 +
          as.integer(.y$lag)
      )
      
      tibble::tibble(
        time_lag_min =
          as.integer(.y$lag) *
          nominal_interval_min_60,
        
        residual_correlation =
          result[["estimate"]],
        
        confidence_low =
          result[["confidence_low"]],
        
        confidence_high =
          result[["confidence_high"]],
        
        n_individuals =
          nrow(.x),
        
        temporal_autocorrelation_detected =
          confidence_low > 0 |
          confidence_high < 0
      )
    }
  ) %>%
  dplyr::ungroup()

print(temporal_autocorrelation_diagnostics_60,n = Inf,width = Inf)
# model               lag time_lag_min residual_correlation confidence_low confidence_high n_individuals temporal_autocorrelation_detected
# Common HFI slope          1           60               0.0190       -0.0181           0.0540            62 FALSE                            
# Common HFI slope          2          120               0.0376       -0.0328           0.102             60 FALSE                            
# Common HFI slope          3          180               0.129        -0.00491          0.260             48 FALSE                            
# Individual HFI slopes     1           60               0.0171       -0.0203           0.0529            62 FALSE                            
# Individual HFI slopes     2          120               0.0405       -0.0287           0.105             60 FALSE                            
# Individual HFI slopes     3          180               0.126        -0.00385          0.253             48 FALSE 

# 2.3 Spatial autocorrelation within individuals ----
# Residuals are first averaged within individual × 1-km cell.
spatial_cell_residuals_60 <-
  residual_data_60 %>%
  dplyr::mutate(
    cell_x =
      floor(
        x_3035 /
          spatial_grid_size_m_60
      ),
    
    cell_y =
      floor(
        y_3035 /
          spatial_grid_size_m_60
      )
  ) %>%
  dplyr::group_by(
    model,
    individual_id,
    cell_x,
    cell_y
  ) %>%
  dplyr::summarise(
    x =
      mean(x_3035),
    
    y =
      mean(y_3035),
    
    residual =
      mean(residual_normal),
    
    n_locations =
      dplyr::n(),
    
    .groups = "drop"
  )


calculate_individual_moran_60 <- function(data) {
  
  if (
    nrow(data) < 6L ||
    stats::sd(data$residual) == 0
  ) {
    return(
      tibble::tibble(
        n_spatial_cells = nrow(data),
        moran_I = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  k_value <- min(
    spatial_k_neighbours_60,
    nrow(data) - 1L
  )
  
  neighbours <-
    spdep::knearneigh(
      cbind(
        data$x,
        data$y
      ),
      k = k_value
    ) %>%
    spdep::knn2nb()
  
  spatial_weights <-
    spdep::nb2listw(
      neighbours,
      style = "W",
      zero.policy = TRUE
    )
  
  moran_test <-
    spdep::moran.mc(
      data$residual,
      spatial_weights,
      nsim =
        n_spatial_permutations_60,
      alternative = "greater",
      zero.policy = TRUE
    )
  
  tibble::tibble(
    n_spatial_cells =
      nrow(data),
    
    moran_I =
      as.numeric(
        moran_test$statistic
      ),
    
    p_value =
      moran_test$p.value
  )
}


spatial_autocorrelation_by_individual_60 <-
  spatial_cell_residuals_60 %>%
  dplyr::group_by(
    model,
    individual_id
  ) %>%
  dplyr::group_modify(
    ~ calculate_individual_moran_60(.x)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(
    model
  ) %>%
  dplyr::mutate(
    adjusted_p_value =
      stats::p.adjust(
        p_value,
        method = "BH"
      ),
    
    positive_spatial_autocorrelation =
      moran_I > 0 &
      adjusted_p_value < 0.05
  ) %>%
  dplyr::ungroup()


spatial_autocorrelation_diagnostics_60 <-
  spatial_autocorrelation_by_individual_60 %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(n_individuals_tested = sum(is.finite(moran_I)),
    
    median_moran_I = stats::median(moran_I,na.rm = TRUE),
    
    n_positive_after_FDR = sum(positive_spatial_autocorrelation,na.rm = TRUE),
    
    proportion_positive_after_FDR =mean(positive_spatial_autocorrelation,na.rm = TRUE),
    
    spatial_autocorrelation_detected = n_positive_after_FDR > 0,
    
    .groups = "drop")

print(spatial_autocorrelation_diagnostics_60,n = Inf,width = Inf)

# model                 n_individuals_tested median_moran_I n_positive_after_FDR proportion_positive_after_FDR spatial_autocorrelation_detected
# Common HFI slope                        62       -0.00532                    0                             0 FALSE                           
# Individual HFI slopes                   62       -0.00670                    0                             0 FALSE

# 2.4 Final diagnostic summary ----
autocorrelation_summary_60 <-
  temporal_autocorrelation_diagnostics_60 %>%
  dplyr::filter(
    lag == 1L
  ) %>%
  dplyr::select(
    model,
    
    temporal_lag1_correlation =
      residual_correlation,
    
    temporal_confidence_low =
      confidence_low,
    
    temporal_confidence_high =
      confidence_high,
    
    temporal_autocorrelation_detected
  ) %>%
  dplyr::left_join(
    spatial_autocorrelation_diagnostics_60,
    by = "model"
  )

print(autocorrelation_summary_60,n = Inf,width = Inf)
# no spatio temporal autocorrelation within the residual of the models have been detected



#------------------------------------------------------------------------------- STEP 3: individual-cluster bootstrap of the common HFI-slope model ----
#' **Steps:**
#' (i) sample individuals with replacement;
#' (ii) retain all transitions from every sampled individual;
#' (iii) assign a new ID and equal total weight to every sampled copy;
#' (iv) refit the common-slope model;
#' (v) estimate 95% bootstrap confidence intervals for the HFI coefficient
#'     and the Q05-Q95 log-odds contrast.


# 3.1 Parameters ----
n_bootstrap_60 <- 1000L       # Use 300L for an exploratory run
bootstrap_seed_60 <- 20260724L
n_bootstrap_cores_60 <- max(1L,min(4L,parallel::detectCores() - 1L))

# 3.2 Fix the original Q05-Q95 HFI contrast ----
# hfi_within_z remains on the global HFI-standard-deviation scale.
hfi_q05_within_60 <- as.numeric(stats::quantile(open_habitat_model_data_60$hfi_within_z,0.05))
hfi_q95_within_60 <- as.numeric(stats::quantile(open_habitat_model_data_60$hfi_within_z,0.95))
hfi_q05_q95_range_60 <- hfi_q95_within_60 - hfi_q05_within_60

# 3.3 Observed estimates ----
observed_hfi_coefficient_60 <- unname(stats::coef(common_hfi_slope_model_ml_60)[["hfi_within_z"]])
observed_log_odds_contrast_60 <- observed_hfi_coefficient_60 * hfi_q05_q95_range_60

# 3.4 Prepare bootstrap sampling ----
individual_data_60 <- split(open_habitat_model_data_60, open_habitat_model_data_60$individual_id, drop = TRUE)
n_individuals_60 <- length(individual_data_60)
set.seed(bootstrap_seed_60)
bootstrap_draws_60 <- replicate(n_bootstrap_60,sample.int(n_individuals_60,n_individuals_60,replace = TRUE),simplify = FALSE)

# 3.5 Bootstrap model formula ----
bootstrap_common_formula_60 <- remain_aerial ~
  cos_diel +
  sin_time +
  age_z +
  duration_z +
  duration_z2 +
  hfi_within_z +
  hfi_between_z +
  dem_elevation_z +
  prop_low_vegetation_5cells_z +
  s(
    bootstrap_individual_id,
    bs = "re")


# 3.6 Run one bootstrap replication ----
run_individual_bootstrap_60 <- function(iteration) {
  
  sampled_individuals <-bootstrap_draws_60[[iteration]]
  
  bootstrap_data <- dplyr::bind_rows(lapply(seq_along(sampled_individuals),
        function(copy_id) {individual_data_60[[sampled_individuals[copy_id]]] %>%
            dplyr::mutate(bootstrap_individual_id =copy_id)})) %>%
    dplyr::mutate(bootstrap_individual_id =factor(bootstrap_individual_id)) %>%
    dplyr::group_by(bootstrap_individual_id) %>%
    dplyr::mutate(bootstrap_weight_raw =1 / dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bootstrap_weight =bootstrap_weight_raw /
        mean(bootstrap_weight_raw))
  
  tryCatch(
    {bootstrap_model <- withCallingHandlers(
        mgcv::gam(
          formula =bootstrap_common_formula_60,
          data =bootstrap_data,
          weights =bootstrap_weight,
          family =stats::binomial(link = "logit"),
          method =gam_selection_method),
        
        # Suppress only the expected warning caused by fractional weights
        warning = function(w) {if (
            grepl("non-integer #successes",conditionMessage(w),
              fixed = TRUE)) {
            invokeRestart("muffleWarning")}})
      
      hfi_coefficient <-unname(stats::coef(bootstrap_model)[["hfi_within_z"]])
      
      tibble::tibble(
        iteration =iteration,
        converged =bootstrap_model$converged,
        hfi_coefficient =hfi_coefficient,
        log_odds_contrast =hfi_coefficient *hfi_q05_q95_range_60)},
    
    error = function(e) {
      tibble::tibble(iteration =iteration,
        converged =FALSE,
        hfi_coefficient =NA_real_,
        log_odds_contrast =NA_real_)})}


# 3.7 Run bootstrap ----
bootstrap_results_60 <- dplyr::bind_rows(
  parallel::mclapply(X =seq_len(n_bootstrap_60),
    FUN =run_individual_bootstrap_60,
    mc.cores =n_bootstrap_cores_60,
    mc.preschedule =FALSE))

successful_bootstrap_results_60 <-bootstrap_results_60 %>%
  dplyr::filter(converged,is.finite(hfi_coefficient),
    is.finite(log_odds_contrast))

# 3.8 Final bootstrap report ----
final_bootstrap_report_60 <- tibble::tibble(
  measure = c("HFI coefficient","Q95-Q05 log-odds contrast"),
  estimate = c(observed_hfi_coefficient_60,observed_log_odds_contrast_60),
  confidence_low = c(stats::quantile( successful_bootstrap_results_60$hfi_coefficient,0.025,names = FALSE),
    stats::quantile(successful_bootstrap_results_60$log_odds_contrast,0.025,names = FALSE)),
  confidence_high = c(stats::quantile(successful_bootstrap_results_60$hfi_coefficient,0.975,names = FALSE),
    stats::quantile(successful_bootstrap_results_60$log_odds_contrast, 0.975,names = FALSE)))

print(final_bootstrap_report_60,n = Inf)
# measure                   estimate confidence_low confidence_high
# HFI coefficient             0.0525       0.000202           0.104
# Q95-Q05 log-odds contrast   0.153        0.000587           0.302
