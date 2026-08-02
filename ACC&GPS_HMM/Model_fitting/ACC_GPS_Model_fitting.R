#' ----------------------------------------------------------------------------- 
# Title: ACC and GPS-classified behaviours dataset comparison and model validation ----
#' Author: Louise Faure
#' Date: 22.07.26
#'
#' **Info:** This script follows "Transition_Matrix_Preparation.R".
#'
#' **Main steps:** 
#' (1) fit model with a common slope response to HFI for gps data and a 
#' "sequential" two stage model for the gps data
#' (2) control the spatio temporal autocorrelation in the pseudo residual of the
#' model. 
#'    (i) generate DHARMa pseudo-residuals for gps and acc models;
#'    (ii) test temporal autocorrelation within each individual and burst at lags 
#'    of 60, 120 and 180 minutes;
#'    (iii) summarise temporal residual correlations across individuals using 
#'    bootstrap 95% confidence intervals;
#'    (iv) aggregate residuals within 1-km spatial cells and calculate Moran’s I 
#'    separately for each individual;
#'    (v) correct spatial test p-values for multiple comparisons using the 
#'    Benjamini–Hochberg FDR procedure;
#'
#' (3) Confirmation of the results of both gps and acc: bootstraping on 
#' individuals to confirm the confidence interval for HFI
#' -----------------------------------------------------------------------------


# Libraries ----
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(lme4)
library(corrplot)

# Data ----
GE_acc_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/acc_weighted.rds")

# Parameters ----




#-------------------------------------------------------------------------------
# STEP 3: individual within-individual HFI slopes for ACC states ----
#'
#' Aerial is the reference state.
#' The model estimates separate HFI slopes for:
#'   (i) resting versus aerial;
#'   (ii) feeding versus aerial.
#'
#' The random intercept and random HFI slope are fitted as independent
#' random-effect components, equivalent to || in lme4.

# Parameters ----
hfi_variable_60 <- "hfi_mean_1000m_z"
acc_backbone_terms_60 <- c("cos_diel_c","sin_diel_c","age_z","duration_z","duration_z2")
acc_environmental_terms_60 <- c("elevation_100m_z","ruggedness_100m_z")
acc_numeric_terms_60 <- c(acc_backbone_terms_60,acc_environmental_terms_60,hfi_variable_60)

# 3.1 Prepare the ACC model dataset and decompose HFI ----
acc_hfi_random_slope_data_60 <- GE_acc_hfi %>%
  dplyr::select(
    transition_destination,individual_id,
    dplyr::all_of(acc_numeric_terms_60)
  ) %>%
  tidyr::drop_na() %>%
  dplyr::filter(
    dplyr::if_all(dplyr::all_of(acc_numeric_terms_60),is.finite)
  ) %>%
  dplyr::mutate(
    transition_destination = stats::relevel(
      factor(transition_destination,levels = c("aerial","resting","feeding")),
      ref = "aerial"
    ),
    individual_id = droplevels(factor(individual_id)),
    hfi_raw_z = .data[[hfi_variable_60]]
  ) %>%
  dplyr::group_by(individual_id) %>%
  dplyr::mutate(
    hfi_between_z = mean(hfi_raw_z),
    hfi_within_z = hfi_raw_z - hfi_between_z,
    n_transitions_individual = dplyr::n(),
    individual_weight_raw = 1 / n_transitions_individual
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    individual_weight = individual_weight_raw / mean(individual_weight_raw),
    individual_intercept_id = factor(individual_id),
    individual_slope_id = factor(individual_id)
  )

# 3.2 Define the common fixed-effect formula ----
acc_hfi_formula_60 <- stats::reformulate(
  termlabels = c(
    acc_backbone_terms_60,
    acc_environmental_terms_60,
    "hfi_within_z",
    "hfi_between_z"
  ),
  response = "transition_destination"
)

# 3.3 Fit a common-slope model and an individual-slope model ----
acc_common_hfi_slope_model_ml_60 <- mclogit::mblogit(
  formula = acc_hfi_formula_60,
  data = acc_hfi_random_slope_data_60,
  random = ~ 1 | individual_intercept_id,
  catCov = "diagonal",
  weights = individual_weight,
  method = "PQL",
  estimator = "ML",
  dispersion = FALSE,
  aggregate = FALSE,
  na.action = stats::na.fail
)

acc_individual_hfi_slope_model_ml_60 <- mclogit::mblogit(
  formula = acc_hfi_formula_60,
  data = acc_hfi_random_slope_data_60,
  random = list(
    ~ 1 | individual_intercept_id,
    ~ 0 + hfi_within_z | individual_slope_id
  ),
  catCov = "diagonal",
  weights = individual_weight,
  method = "PQL",
  estimator = "ML",
  dispersion = FALSE,
  aggregate = FALSE,
  na.action = stats::na.fail
)

# 3.4 Compile the model comparison ----
acc_hfi_slope_model_comparison_60 <- tibble::tibble(
  model = c("Common within-individual HFI slopes",
            "Individual within-individual HFI slopes"),
  AIC = c(
    stats::AIC(acc_common_hfi_slope_model_ml_60),
    stats::AIC(acc_individual_hfi_slope_model_ml_60)
  ),
  model_df = c(
    attr(stats::logLik(acc_common_hfi_slope_model_ml_60),"df"),
    attr(stats::logLik(acc_individual_hfi_slope_model_ml_60),"df")
  ),
  converged = c(
    isTRUE(acc_common_hfi_slope_model_ml_60$converged),
    isTRUE(acc_individual_hfi_slope_model_ml_60$converged)
  )
) %>%
  dplyr::mutate(delta_AIC = AIC - min(AIC)) %>%
  dplyr::arrange(AIC)

# 3.5 Extract individual random effects ----
# 3.5 Identify the category-specific random-effect blocks ----
random_effect_blocks_60 <- as.list(
  acc_individual_hfi_slope_model_ml_60$random.effects
)

random_effect_labels_60 <- vapply(
  acc_individual_hfi_slope_model_ml_60$VarCov,
  function(variance_matrix) {
    block_names <- rownames(variance_matrix)
    if(length(block_names) != 1L) {
      stop("Each random-effect block was expected to contain one term.")
    }
    block_names[[1]]
  },
  character(1)
)

names(random_effect_blocks_60) <- random_effect_labels_60

control_random_effect_blocks_60 <- tibble::tibble(
  block = seq_along(random_effect_blocks_60),
  block_label = random_effect_labels_60,
  n_individuals = vapply(random_effect_blocks_60,length,integer(1))
)

expected_slope_blocks_60 <- c(
  "aerial~hfi_within_z",
  "resting~hfi_within_z",
  "feeding~hfi_within_z"
)

if(!all(expected_slope_blocks_60 %in% random_effect_labels_60)) {
  stop(
    "The expected category-specific HFI slope blocks were not found. Available blocks: ",
    paste(random_effect_labels_60,collapse = ", ")
  )
}

# 3.6 Extract one random-slope block and retain individual names ----
extract_random_slope_block_60 <- function(block_label,slope_name) {
  random_effect_block <- random_effect_blocks_60[[block_label]]
  
  individual_names <- rownames(random_effect_block)
  if(is.null(individual_names)) {
    individual_names <- names(random_effect_block)
  }
  if(is.null(individual_names)) {
    individual_names <- levels(
      acc_hfi_random_slope_data_60$individual_id
    )
  }
  
  if(length(individual_names) != length(random_effect_block)) {
    stop("Individual names and random effects have different lengths for ",block_label,".")
  }
  
  tibble::tibble(
    individual_id = as.character(individual_names),
    !!slope_name := as.numeric(random_effect_block)
  )
}

aerial_random_slopes_60 <- extract_random_slope_block_60(
  block_label = "aerial~hfi_within_z",
  slope_name = "aerial_random_slope"
)

resting_random_slopes_60 <- extract_random_slope_block_60(
  block_label = "resting~hfi_within_z",
  slope_name = "resting_random_slope"
)

feeding_random_slopes_60 <- extract_random_slope_block_60(
  block_label = "feeding~hfi_within_z",
  slope_name = "feeding_random_slope"
)

# 3.7 Extract the population-level HFI slopes ----
fixed_effect_matrix_60 <- acc_individual_hfi_slope_model_ml_60$coefmat

if(!"hfi_within_z" %in% colnames(fixed_effect_matrix_60)) {
  stop("The fixed hfi_within_z coefficient was not found.")
}

resting_fixed_row_60 <- which(
  tolower(rownames(fixed_effect_matrix_60)) == "resting"
)

feeding_fixed_row_60 <- which(
  tolower(rownames(fixed_effect_matrix_60)) == "feeding"
)

if(length(resting_fixed_row_60) != 1L ||
   length(feeding_fixed_row_60) != 1L) {
  stop(
    "The resting and feeding fixed-effect equations could not be identified. ",
    "Available equations: ",
    paste(rownames(fixed_effect_matrix_60),collapse = ", ")
  )
}

population_resting_HFI_slope_60 <- fixed_effect_matrix_60[
  resting_fixed_row_60,
  "hfi_within_z"
]

population_feeding_HFI_slope_60 <- fixed_effect_matrix_60[
  feeding_fixed_row_60,
  "hfi_within_z"
]

# 3.8 Calculate individual resting-versus-aerial and feeding-versus-aerial slopes ----
individual_acc_hfi_slopes_60 <- aerial_random_slopes_60 %>%
  dplyr::inner_join(
    resting_random_slopes_60,
    by = "individual_id"
  ) %>%
  dplyr::inner_join(
    feeding_random_slopes_60,
    by = "individual_id"
  ) %>%
  dplyr::mutate(
    resting_HFI_slope =
      population_resting_HFI_slope_60 +
      resting_random_slope -
      aerial_random_slope,
    
    feeding_HFI_slope =
      population_feeding_HFI_slope_60 +
      feeding_random_slope -
      aerial_random_slope
  )

# 3.9 Classify and rank individuals by their HFI responses ----
individual_acc_hfi_slopes_60 <- individual_acc_hfi_slopes_60 %>%
  dplyr::mutate(
    response_group = dplyr::case_when(
      feeding_HFI_slope > 0 & resting_HFI_slope > 0 ~
        "1. Positive feeding and resting",
      
      feeding_HFI_slope > 0 & resting_HFI_slope <= 0 ~
        "2. Positive feeding only",
      
      feeding_HFI_slope <= 0 & resting_HFI_slope > 0 ~
        "3. Positive resting only",
      
      TRUE ~
        "4. All other individuals"
    ),
    
    response_group_order = dplyr::case_when(
      feeding_HFI_slope > 0 & resting_HFI_slope > 0 ~ 1L,
      feeding_HFI_slope > 0 & resting_HFI_slope <= 0 ~ 2L,
      feeding_HFI_slope <= 0 & resting_HFI_slope > 0 ~ 3L,
      TRUE ~ 4L
    ),
    
    ranking_score = dplyr::case_when(
      response_group_order == 1L ~
        (feeding_HFI_slope + resting_HFI_slope) / 2,
      
      response_group_order == 2L ~
        feeding_HFI_slope,
      
      response_group_order == 3L ~
        resting_HFI_slope,
      
      TRUE ~
        (feeding_HFI_slope + resting_HFI_slope) / 2
    )
  ) %>%
  dplyr::arrange(
    response_group_order,
    dplyr::desc(ranking_score),
    individual_id
  ) %>%
  dplyr::select(
    individual_id,
    response_group,
    feeding_HFI_slope,
    resting_HFI_slope
  )

base::print(
  as.data.frame(individual_acc_hfi_slopes_60),
  row.names = FALSE
)
# individual_id                  response_group feeding_HFI_slope resting_HFI_slope
# Trimmis20 (eobs 7041) 1. Positive feeding and resting      0.1053125398      0.1378792039
# Flüela1 21 (eobs 6995) 1. Positive feeding and resting      0.0701610343      0.1027323892
# Mals1_23 (eobs 11919) 1. Positive feeding and resting      0.0663891217      0.0989738685
# Almen19 (eobs 7001) 1. Positive feeding and resting      0.0522503872      0.0848121467
# Ettenberg22 (eobs 10539) 1. Positive feeding and resting      0.0488477783      0.0814053615
# Siat20 (eobs 7037) 1. Positive feeding and resting      0.0292900375      0.0618579766
# Punteglias20 (eobs 6483) 1. Positive feeding and resting      0.0232082398      0.0557396070
# Fahrntal19 (eobs 7014) 1. Positive feeding and resting      0.0114228044      0.0439741059
# Flüela20 (eobs 7040) 1. Positive feeding and resting      0.0080087521      0.0405862136
# Dischma1 19 (eobs 7006) 1. Positive feeding and resting      0.0009252196      0.0334809968
# Flüela19 (eobs 7007)        3. Positive resting only     -0.0029403722      0.0296080244
# Mals2_23 (eobs 11920)        3. Positive resting only     -0.0042003906      0.0283519696
# Dischma2 19 (eobs 7009)        3. Positive resting only     -0.0081705512      0.0243820336
# Grabernock23 (eobs 11916)        3. Positive resting only     -0.0091762003      0.0233820757
# Vernuga22 (eobs 10538)        3. Positive resting only     -0.0229385033      0.0096198329
# Umbrail18 (eobs 5859)        3. Positive resting only     -0.0321759301      0.0003773899
# Krn20 (eobs 7549)        4. All other individuals     -0.0392847684     -0.0067308647
# Gaming24 (eobs 11914)        4. All other individuals     -0.0429579741     -0.0103988527
# Matsch19 (eobs 7035)        4. All other individuals     -0.0438700119     -0.0113209823
# Schlappin1 18 (eobs 5858)        4. All other individuals     -0.0455205095     -0.0129502346
# ValSozzine21 (eobs 7500)        4. All other individuals     -0.0464473955     -0.0138930816
# Sinestra1 19 (eobs 7003)        4. All other individuals     -0.0481524272     -0.0155958337
# Sinestra2 19 (eobs 7005)        4. All other individuals     -0.0494163899     -0.0168560497
# Laas1_23 (eobs 11917)        4. All other individuals     -0.0508440229     -0.0182868590
# Tabland22 (eobs 10534)        4. All other individuals     -0.0524415199     -0.0198941938
# Trenzeira19 (eobs 5858)        4. All other individuals     -0.0604222920     -0.0278715814
# Johnsbach2_23 (eobs 11913)        4. All other individuals     -0.0628426641     -0.0302913396
# Grabernock21 (eobs 7506)        4. All other individuals     -0.0637383870     -0.0311780915
# Nalps18 (eobs 5860)        4. All other individuals     -0.0658667301     -0.0333057518
# Vrata20 (eobs 7551)        4. All other individuals     -0.0709867993     -0.0384381790
# Laas2_23 (eobs 11918)        4. All other individuals     -0.0710863599     -0.0385420015
# Kastelbell19 (eobs 7034)        4. All other individuals     -0.0722691318     -0.0397172417
# Stürfis20 (eobs 7049)        4. All other individuals     -0.0723984150     -0.0398256265
# Reschen21 (eobs 7503)        4. All other individuals     -0.0727957428     -0.0402628092
# Adamello20 (eobs 7548)        4. All other individuals     -0.0758050162     -0.0432495926
# Sampuoir2 19 (eobs 6462)        4. All other individuals     -0.0763130119     -0.0437697051
# Sampuoir1 19 (eobs 5943)        4. All other individuals     -0.0790097092     -0.0464719535
# Lischana22 (eobs 5941)        4. All other individuals     -0.0797099684     -0.0471620138
# Tuors2 19 (eobs 7011)        4. All other individuals     -0.0805596917     -0.0480088636
# Grosio 19 (eobs 7000)        4. All other individuals     -0.0869214433     -0.0543740690
# ValGrande19 (eobs 7033)        4. All other individuals     -0.0876250617     -0.0550829065
# Johnsbach1_23 (eobs 11645)        4. All other individuals     -0.0916035892     -0.0590542500
# Burgum21 (eobs 7502)        4. All other individuals     -0.0920845494     -0.0595407054
# Cornasc20 (eobs 7039)        4. All other individuals     -0.0922545960     -0.0597061797
# Lassingbach24 (eobs 11915)        4. All other individuals     -0.0938277080     -0.0612775737
# Art San Romerio18 (eobs 5941)        4. All other individuals     -0.0972576729     -0.0647061479
# Tasna18 (eobs 5940)        4. All other individuals     -0.1027035460     -0.0701502534
# Flüela2 21 (eobs 7043)        4. All other individuals     -0.1047659800     -0.0722180207
# Reschen20 (eobs 7556)        4. All other individuals     -0.1087913542     -0.0762436675
# Valdidentro_Braulio20 (eobs 7581)        4. All other individuals     -0.1101782732     -0.0776206663
# Schlappin22 (eobs 5944)        4. All other individuals     -0.1159976377     -0.0834658588
# Torta19 (eobs 7002)        4. All other individuals     -0.1211010879     -0.0885511996
# Güstizia18 (eobs 5942)        4. All other individuals     -0.1270445629     -0.0945019755
# Tuors1 19 (eobs 7010)        4. All other individuals     -0.1291542956     -0.0966297412
# Seta19 (eobs 5796)        4. All other individuals     -0.1453008623     -0.1127686350
# Sils20 (eobs 7038)        4. All other individuals     -0.1732960288     -0.1407582592
# Avers20 (eobs 7101)        4. All other individuals     -0.1991352223     -0.1666066502
#Nalps19 (eobs 5861)        4. All other individuals     -0.2356461207     -0.2031114214

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
# model                 temporal_lag1_correlation temporal_confidence_low temporal_confidence_high temporal_autocorrelation_detected
# Common HFI slope                         0.0190                 -0.0181                   0.0540 FALSE                            
# Individual HFI slopes                    0.0171                 -0.0203                   0.0529 FALSE                            
# n_individuals_tested median_moran_I n_positive_after_FDR proportion_positive_after_FDR spatial_autocorrelation_detected
#                  62       -0.00532                    0                             0 FALSE                           
#                  62       -0.00670                    0                             0 FALSE  


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