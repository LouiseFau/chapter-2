#'-------------------------------------------------------------------------------
#' Title: iSSF fitted per individuals ----
#' Authors : Louise Faure
#' Date : 29.09.26
#' **Info:** this script follow the Data_processing_annotation.R script where I 
#' generated random step based on a gamma and uniform distribution for step lenght
#' and turning angle, extracted the environmental values below each data point, and 
#' created several similarity indexes. This models used the covariates identified in 
#' "iSSF_covariate_selection.R", by adding in a two step approach the influence of
#' natal havitat on the individual selection of landing sites. 
#' **Expected results:**
#' (1) Individuals differ in their avoidance of settlements
#' (2) Individuals from more built natal territories avoid settlements less
#' (3) The natal descriptor reduces the variance between individuals
#' **Purpose:**
#' (1) dataset preparation (variable standardization)
#' (2) fit one iSSF per individual and extract the individual settlement density slope
#' (3) regress the individual slopes on the natal descriptor (random-effects
#' meta-regression fitted by maximum likelihood in base R, territory as random effect)
#' (4) display only the statistical indexes that answer the expected results
#' ------------------------------------------------------------------------------

# Libraries ----
library(tidyverse)
library(glmmTMB)

# Parameters ----
minimum_landings_60 <- 30L
fixed_stratum_sd_60 <- 1e3
model_control_60 <- glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000))
natal_descriptor_60 <- "natal_q90_NT"

# Colors ----
clr <- oce::oceColorsPalette(100)[9] #was [2] before
clr_light <- oce::oceColorsPalette(100)[10]
clr2 <- oce::oceColorsPalette(100)[80]

# Golden eagle dataset ----
annotated_data <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/issf_generated_observed_location_annotated(2).rds")

#------------------------------------------------------------------------------- STEP 1: dataset preparation ----
#' **Steps:**
#' (i) retain complete choice sets and individuals with at least 30 landings;
#' (ii) standardize environmental covariates and prepare movement terms.

# 1.1 Prepare complete choice sets and retain supported individuals ----
environmental_covariates_60 <- c("settlement_density","elevation_100m","ruggedness_100m","prop_forest_5cells","prop_low_vegetation_5cells","distance_to_ridgeline_100m")
required_columns_60 <- c("used","stratum","individual.local.identifier","step_length_km","turning_angle_rad","territory_id",natal_descriptor_60,environmental_covariates_60)
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
#' Population standardization: individual slopes are expressed on the same scale.
standardization_parameters_60 <- tibble::tibble(
  variable = environmental_covariates_60,
  center = vapply(environmental_covariates_60,\(x) mean(data_model_60[[x]]),numeric(1)),
  scale = vapply(environmental_covariates_60,\(x) stats::sd(data_model_60[[x]]),numeric(1))
)

data_model_60 <- data_model_60 %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(environmental_covariates_60),~ as.numeric(scale(.x)),.names = "{.col}_z"),
                log_step_length_km = log(step_length_km),
                cos_turning_angle = cos(turning_angle_rad))


#------------------------------------------------------------------------------- STEP 2: fit one iSSF per individual ----
#' **Steps:**
#' (i) fit the selected iSSF separately for each individual (no individual random slope);
#' (ii) extract the individual settlement density slope, its standard error and the natal descriptor.

# 2.1 Individual model ----
formula_individual_60 <- used ~ -1 + settlement_density_z + elevation_100m_z + prop_low_vegetation_5cells_z + ruggedness_100m_z + step_length_km + log_step_length_km + cos_turning_angle + (1 | stratum_ID)

fit_individual_issf_60 <- function(data){
  tryCatch(glmmTMB::glmmTMB(formula = formula_individual_60,family = poisson(link = "log"),data = data,map = list(theta = factor(NA)),start = list(theta = log(fixed_stratum_sd_60)),control = model_control_60),error = \(e) NULL)}

# 2.2 Fit one model per individual ----
individual_data_60 <- lapply(split(data_model_60,data_model_60$animal_ID,drop = TRUE),droplevels)
individual_models_60 <- lapply(individual_data_60,fit_individual_issf_60)

# 2.3 Extract individual slopes ----
individual_slopes_60 <- dplyr::bind_rows(lapply(names(individual_models_60),\(ind){
  m <- individual_models_60[[ind]]; d <- individual_data_60[[ind]]
  ct <- if(is.null(m)) NULL else summary(m)$coefficients$cond
  tibble::tibble(
    individual.local.identifier = as.character(d$individual.local.identifier[1]),
    territory_id = d$territory_id[1],
    n_landings = sum(d$used == 1L),
    converged = !is.null(m) && isTRUE(m$fit$convergence == 0) && isTRUE(m$sdr$pdHess),
    slope = if(is.null(ct)) NA_real_ else unname(ct["settlement_density_z","Estimate"]),
    slope_se = if(is.null(ct)) NA_real_ else unname(ct["settlement_density_z","Std. Error"]),
    natal_value = d[[natal_descriptor_60]][1])
}))

controle_individual_fits_60 <- individual_slopes_60 %>%
  dplyr::summarise(n_individuals = dplyr::n(),n_converged = sum(converged),n_missing_natal = sum(is.na(natal_value)),
                   median_slope_se = stats::median(slope_se,na.rm = TRUE),max_slope_se = max(slope_se,na.rm = TRUE))
controle_unstable_individuals_60 <- individual_slopes_60 %>%
  dplyr::filter(!converged | slope_se > 5 * stats::median(slope_se,na.rm = TRUE)) %>% dplyr::arrange(dplyr::desc(slope_se))


#------------------------------------------------------------------------------- STEP 3: regress individual slopes on the natal descriptor ----
#' **Model:** slope_i = X_i beta + u_territory + u_individual + e_i
#' e_i ~ N(0, slope_se_i^2): sampling error of step 2, known;
#' u_territory ~ N(0, tau2_territory): individuals from the same natal territory share a polygon;
#' u_individual ~ N(0, tau2_individual): remaining differences between individuals.
#' Each slope is therefore weighted by its precision. Fitted by maximum likelihood,
#' equivalent to a random-effects meta-regression (metafor::rma.mv, method = "ML").
#' Model A: no predictor (expected result 1). Model B: natal descriptor (expected results 2 and 3).

# 3.1 Prepare the individual table ----
two_step_data_60 <- individual_slopes_60 %>%
  dplyr::filter(converged,is.finite(slope),is.finite(slope_se),!is.na(natal_value),!is.na(territory_id)) %>%
  dplyr::mutate(natal_value_z = as.numeric(scale(natal_value)),territory_id = factor(territory_id))

# 3.2 Maximum likelihood fit of the meta-regression ----
fit_meta_regression_60 <- function(moderator_formula,data){
  X <- stats::model.matrix(moderator_formula,data); y <- data$slope
  Z <- stats::model.matrix(~ 0 + territory_id,data)
  sigma_matrix <- \(log_tau2) diag(data$slope_se^2 + exp(log_tau2[2]),nrow(data)) + exp(log_tau2[1]) * tcrossprod(Z)
  gls <- \(S){XtSi <- crossprod(X,solve(S)); cov_beta <- solve(XtSi %*% X)
  list(beta = stats::setNames(drop(cov_beta %*% XtSi %*% y),colnames(X)),cov_beta = cov_beta)}
  negative_loglik <- \(log_tau2){S <- sigma_matrix(log_tau2); r <- y - drop(X %*% gls(S)$beta)
  0.5 * (as.numeric(determinant(S,logarithm = TRUE)$modulus) + drop(crossprod(r,solve(S,r))) + length(y) * log(2 * pi))}
  opt <- stats::optim(c(log(0.1),log(0.1)),negative_loglik,method = "L-BFGS-B",lower = c(-20,-20),upper = c(5,5))
  g <- gls(sigma_matrix(opt$par))
  list(beta = g$beta,se = stats::setNames(sqrt(diag(g$cov_beta)),colnames(X)),
       tau2 = stats::setNames(exp(opt$par),c("territory","individual")),
       logLik = -opt$value,AIC = 2 * opt$value + 2 * (ncol(X) + 2),convergence = opt$convergence)
}

# 3.3 Fit the two models ----
two_step_models_60 <- list(
  A_no_predictor = fit_meta_regression_60(~ 1,two_step_data_60),
  B_natal = fit_meta_regression_60(~ natal_value_z,two_step_data_60))

# 3.4 Heterogeneity test (Cochran's Q): do slopes vary more than their sampling error? ----
fixed_effect_A_60 <- stats::lm(slope ~ 1,data = two_step_data_60,weights = 1 / slope_se^2)
heterogeneity_Q_60 <- sum(stats::weights(fixed_effect_A_60) * stats::residuals(fixed_effect_A_60)^2)


#------------------------------------------------------------------------------- STEP 4: results ----
#' Expected result 1: heterogeneity_p < 0.05 (individuals differ more than their imprecision predicts)
#' Expected result 2: natal_coefficient > 0 with a CI excluding 0 (less avoidance when the NT is more built)
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
  prop_individuals_avoiding = stats::pnorm(0,mean = mean_slope,sd = between_individual_sd),
  heterogeneity_Q = heterogeneity_Q_60,
  heterogeneity_p = stats::pchisq(heterogeneity_Q_60,df = nrow(two_step_data_60) - 1,lower.tail = FALSE),
  # expected result 2
  natal_coefficient = unname(model_B_60$beta["natal_value_z"]),
  natal_CI_95 = sprintf("[%.3f; %.3f]",natal_coefficient - 1.96 * model_B_60$se["natal_value_z"],natal_coefficient + 1.96 * model_B_60$se["natal_value_z"]),
  natal_p = 2 * stats::pnorm(-abs(natal_coefficient / model_B_60$se["natal_value_z"])),
  # expected result 3
  variance_explained_percent = 100 * (1 - sum(model_B_60$tau2) / sum(model_A_60$tau2)),
  delta_AIC_B_minus_A = model_B_60$AIC - model_A_60$AIC,
  LRT_p = stats::pchisq(2 * (model_B_60$logLik - model_A_60$logLik),df = 1,lower.tail = FALSE))

print(two_step_results_60,width = Inf)
