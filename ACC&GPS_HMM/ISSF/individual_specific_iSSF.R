#'-------------------------------------------------------------------------------
#' Title: iSSF and ilmm fitted per individuals ----
#' Authors : Louise Faure
#' Date : 29.09.26
#' **Info:** this script follow the Data_processing_annotation.R script where I 
#' generated random step based on a gamma and uniform distribution for step lenght
#' and turning angle, extracted the environmental values below each data point, and 
#' created several similarity indexes. This models used the covariates identified in 
#' "iSSF_covariate_selection.R", by adding in a two step approach the influence of
#' natal havitat on the individual selection of landing sites. 
#' **Research questions:**
#' (1) Do individuals differ in their response to settlement densities (persistence
#' in flight and avoidance of human-dominated landscape to land)
#' (2) Are these differences resulting from differences in anthropogenic exposition
#' at the natal territory?
#' (3) If not, is the absence of detectable influence of natal exposition an absence
#' of effect or the result of noise ? 
#' **Purpose:**
#' (1) dataset preparation, one model per individual, compilation of heterogeneity of the
#' individual settlement density slopes (Cochran's Q, I2, between-individual SD).
#' (2) natal descriptor and meta-regression of the individual slopes on it
#' (territory and individual random effects, slopes weighted by their precision).
#' (3): natal effect (coefficient, CI 95 %, variance explained, LRT).
#' (4): statistical power of the negative result: observed effect, smallest
#' detectable effect, interesting effect (a quarter less avoidance between the least
#' and the most built natal territories), equivalence test.
#' ------------------------------------------------------------------------------

# Libraries ----
library(tidyverse)
library(sf)
library(terra)
library(survival)

# Paths ----
annotated_path_60 <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/issf_generated_observed_location_annotated(2).rds"
persistence_path_60 <- "/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/gps_20_weighted.rds"
id_lookup_path_60 <- "/Users/louisefaure/Desktop/dossier sans titre/donnees aigles gps burst/gps_bursts_raw_move2.rds"
natal_polygon_path_60 <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/natal_polygons/natal_polygons.gpkg"
settlement_density_60 <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/settlement_density_1km2_100m.tif")

# Parameters ----
minimum_landings_60 <- 30L
natal_descriptor_60 <- "natal_q90_NT"
confidence_critical_value_60 <- stats::qnorm(0.975)
power_level_60 <- 0.80
interesting_reduction_60 <- 0.25                    # a quarter less avoidance between the least and the most built natal territories
natal_contrast_probabilities_60 <- c(0.10,0.90)     # least and most built natal territories (quantiles over territories)


#------------------------------------------------------------------------------- STEP 1: do individuals differ in their response to settlement density? ----
#' **Steps:**
#' (i) prepare natal territories by creating a territory_id (that can be used later as a random effect) and natal_q90_NT (which is the 90th percentil of
#' settlement density in natal territory);
#' (ii) iSSF dataset: filtration of individuals with at least 30 landings, standardization of variables;
#' (iii) persistence dataset: GPS transitions from an aerial location (already standardized);
#' (iv) one model per individual and its settlement density slope;
#' (v) heterogeneity of the slopes: Q, I2, between-individual SD.

# 1.1 Natal territories ----
id_lookup_60 <- move2::mt_track_data(readRDS(id_lookup_path_60)) %>%
  dplyr::transmute(id = as.character(individual_id),key = stringr::str_squish(as.character(individual_local_identifier)))

NT_poly_60 <- sf::st_read(natal_polygon_path_60,quiet = TRUE) %>% dplyr::transmute(id = as.character(id)) %>%
  sf::st_transform(3035) %>% sf::st_make_valid() %>% sf::st_collection_extract("POLYGON") %>%
  dplyr::inner_join(id_lookup_60,by = "id") %>%
  dplyr::mutate(territory_id = sapply(sf::st_equals_exact(.,.,par = 1),min))

weighted_quantile_60 <- function(v,w,p) {o <- order(v); v[o][which(cumsum(w[o]) / sum(w) >= p)[1]]}

natal_territories_60 <- terra::extract(settlement_density_60,terra::vect(NT_poly_60),exact = TRUE) %>%
  stats::setNames(c("row","value","fraction")) %>% dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(key = NT_poly_60$key[row],territory_id = NT_poly_60$territory_id[row],value = round(value,6)) %>%
  dplyr::group_by(key,territory_id) %>%
  dplyr::summarise(natal_q90_NT = weighted_quantile_60(value,fraction,0.90),.groups = "drop")

# 1.2 iSSF dataset: usable rows, complete choice sets, supported individuals, standardization ----
environmental_covariates_60 <- c("settlement_density","elevation_100m","ruggedness_100m","prop_low_vegetation_5cells")

issf_data_60 <- readRDS(annotated_path_60) %>%
  dplyr::mutate(used = as.integer(used),key = stringr::str_squish(as.character(individual.local.identifier))) %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(environmental_covariates_60),~ !is.na(.x) & is.finite(.x))) %>%
  dplyr::group_by(stratum) %>%
  dplyr::filter(sum(used == 1L) == 1L,sum(used == 0L) >= 1L) %>%
  dplyr::group_by(key) %>%
  dplyr::filter(sum(used == 1L) >= minimum_landings_60) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(environmental_covariates_60),~ as.numeric(scale(.x)),.names = "{.col}_z"),
                log_step_length_km = log(step_length_km),
                cos_turning_angle = cos(turning_angle_rad))

# 1.3 Persistence dataset: transitions from an aerial location (GPS, 20 min, variables standardized upstream) ----
backbone_terms_60 <- c("cos_diel_c","sin_diel_c","duration_z","duration_z2")
environmental_terms_60 <- c("elevation_100m_z","ruggedness_100m_z")
hfi_variable_60 <- "settlement_density_z"
persistence_terms_60 <- c(backbone_terms_60,environmental_terms_60,hfi_variable_60)

persistence_data_60 <- readRDS(persistence_path_60) %>%
  dplyr::mutate(key = stringr::str_squish(as.character(individual_id)),remain_aerial = as.integer(remain_aerial)) %>%
  dplyr::select(key,remain_aerial,dplyr::all_of(persistence_terms_60)) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(key) %>%
  dplyr::filter(sum(remain_aerial == 0L) >= minimum_landings_60) %>%
  dplyr::ungroup()

# 1.4 One model per individual and its settlement density slope ----
#' stratum is already unique across individuals (landing_id is a row number over the whole dataset).
formula_issf_60 <- used ~ settlement_density_z + elevation_100m_z + prop_low_vegetation_5cells_z + ruggedness_100m_z + step_length_km + log_step_length_km + cos_turning_angle + survival::strata(stratum)
formula_persistence_60 <- stats::as.formula(paste("remain_aerial ~",paste(persistence_terms_60,collapse = " + ")))

fit_issf_60 <- function(d) tryCatch(survival::clogit(formula_issf_60,data = droplevels(d)),error = \(e) NULL)
fit_persistence_60 <- function(d) tryCatch(stats::glm(formula_persistence_60,family = stats::binomial(link = "logit"),data = d),error = \(e) NULL)

extract_slope_60 <- function(models,data,n_landings,converged,coefficients,model_name) dplyr::bind_rows(lapply(names(models),\(ind){
  m <- models[[ind]]; ct <- if(is.null(m)) NULL else coefficients(m)
  tibble::tibble(model = model_name,key = ind,n_landings = n_landings(data[[ind]]),
                 converged = !is.null(m) && converged(m) && "settlement_density_z" %in% rownames(ct),
                 slope = if(is.null(ct)) NA_real_ else unname(ct["settlement_density_z","Estimate"]),
                 slope_se = if(is.null(ct)) NA_real_ else unname(ct["settlement_density_z","Std. Error"]))}))

issf_individual_data_60 <- split(issf_data_60,issf_data_60$key,drop = TRUE)
issf_models_60 <- lapply(issf_individual_data_60,fit_issf_60)
persistence_individual_data_60 <- split(persistence_data_60,persistence_data_60$key,drop = TRUE)
persistence_models_60 <- lapply(persistence_individual_data_60,fit_persistence_60)

#' clogit stores coefficients and their covariance separately: they are put in the same layout as
#' the glm coefficient table so that both models are extracted by the same function.
coefficients_clogit_60 <- function(m) cbind(Estimate = stats::coef(m),`Std. Error` = sqrt(diag(stats::vcov(m))))

individual_slopes_60 <- dplyr::bind_rows(
  extract_slope_60(issf_models_60,issf_individual_data_60,\(d) sum(d$used == 1L),\(m) all(is.finite(stats::coef(m))) && all(is.finite(sqrt(diag(stats::vcov(m))))),coefficients_clogit_60,"iSSF"),
  extract_slope_60(persistence_models_60,persistence_individual_data_60,\(d) sum(d$remain_aerial == 0L),\(m) isTRUE(m$converged),\(m) summary(m)$coefficients,"persistence")) %>%
  dplyr::left_join(natal_territories_60,by = "key") %>%
  dplyr::filter(converged,is.finite(slope),is.finite(slope_se),!is.na(territory_id))

slopes_by_model_60 <- lapply(split(individual_slopes_60,individual_slopes_60$model),\(d) dplyr::mutate(d,territory_id = factor(territory_id)))

# 1.5 Do individuals differ? ----
#' Each eagle has one slope and one margin of error on that slope, both from step 1.4. Do the
#' eagles really behave differently, or do their slopes differ only because none of them is
#' measured exactly? Nothing else is asked here: no natal territory, no sibling structure.
#' The model treats an eagle's slope as the average of all eagles, plus what is specific to that
#' eagle, plus its own measurement error. That last part is not guessed by the model: it is given,
#' since step 1.4 computed it. Whatever spread is left once it is removed is therefore real
#' difference between eagles (between_individual_sd), and an eagle with many landings (small
#' margin of error) weighs more than a poorly documented one.
#' Cochran's Q measures how far the slopes fall from their average, each counted the more it is
#' precise. If all eagles were identical, Q would be about the number of eagles minus one; a
#' clearly larger Q means real differences. I2 turns this into a percentage: the share of the
#' observed differences between eagles that is real rather than measurement error.
#' cluster is used only in STEP 2, where siblings sharing a natal polygon must not count twice.
fit_meta_regression_60 <- function(moderator_formula,data,cluster = NULL){
  X <- stats::model.matrix(moderator_formula,data); y <- data$slope
  Z <- if(is.null(cluster)) NULL else stats::model.matrix(~ 0 + factor(cluster))
  sigma_matrix <- \(log_tau2) diag(data$slope_se^2 + exp(log_tau2[1]),nrow(data)) + if(is.null(Z)) 0 else exp(log_tau2[2]) * tcrossprod(Z)
  gls <- \(S){XtSi <- crossprod(X,solve(S)); cov_beta <- solve(XtSi %*% X)
  list(beta = stats::setNames(drop(cov_beta %*% XtSi %*% y),colnames(X)),cov_beta = cov_beta)}
  negative_loglik <- \(log_tau2){S <- sigma_matrix(log_tau2); r <- y - drop(X %*% gls(S)$beta)
  0.5 * (as.numeric(determinant(S,logarithm = TRUE)$modulus) + drop(crossprod(r,solve(S,r))) + length(y) * log(2 * pi))}
  n_tau <- if(is.null(Z)) 1L else 2L
  opt <- stats::optim(rep(log(0.1),n_tau),negative_loglik,method = "L-BFGS-B",lower = rep(-20,n_tau),upper = rep(5,n_tau))
  g <- gls(sigma_matrix(opt$par))
  list(beta = g$beta,se = stats::setNames(sqrt(diag(g$cov_beta)),colnames(X)),
       tau2 = stats::setNames(exp(opt$par),c("individual","territory")[seq_len(n_tau)]),
       logLik = -opt$value,AIC = 2 * opt$value + 2 * (ncol(X) + n_tau),convergence = opt$convergence)
}

individual_variation_models_60 <- lapply(slopes_by_model_60,\(d) fit_meta_regression_60(~ 1,d))

heterogeneity_60 <- dplyr::bind_rows(lapply(names(slopes_by_model_60),\(m){
  d <- slopes_by_model_60[[m]]; A <- individual_variation_models_60[[m]]
  w <- 1 / d$slope_se^2; Q <- sum(w * (d$slope - sum(w * d$slope) / sum(w))^2); df <- nrow(d) - 1L
  tibble::tibble(model = m,n_individuals = nrow(d),
                 mean_slope = unname(A$beta["(Intercept)"]),
                 median_slope_se = stats::median(d$slope_se),
                 between_individual_sd = sqrt(A$tau2[["individual"]]),
                 signal_to_noise = between_individual_sd / median_slope_se,
                 heterogeneity_Q = Q,df = df,heterogeneity_p = stats::pchisq(Q,df,lower.tail = FALSE),
                 I2_percent = 100 * max(0,(Q - df) / Q),
                 prop_individuals_mean_sign = stats::pnorm(0,mean_slope,between_individual_sd,lower.tail = mean_slope < 0))}))

print(heterogeneity_60,width = Inf)
#' Expected: persistence heterogeneity_p > 0.05 and I2 close to 0 (no detectable differences between
#' individuals; signal_to_noise says whether the slopes could have revealed them);
#' iSSF heterogeneity_p < 0.05 and I2 clearly above 0 (real differences). 


#------------------------------------------------------------------------------- STEP 2: iSSF: natal descriptor and meta-regression ----
#' Model B adds natal_q90_NT (standardized among individuals) as the single predictor of the
#' individual slopes. Expected under the hypothesis: positive coefficient (less avoidance for
#' individuals from more built natal territories). Both models are now clustered by natal
#' territory, because siblings share a polygon and therefore the same natal value: without it the
#' natal coefficient would be credited with more independent information than there is.
#' Model A is refitted here with that clustering so that the two models stay comparable.

two_step_data_60 <- slopes_by_model_60$iSSF %>% dplyr::mutate(natal_value_z = as.numeric(scale(.data[[natal_descriptor_60]])))

model_A_60 <- fit_meta_regression_60(~ 1,two_step_data_60,cluster = two_step_data_60$territory_id)
model_B_60 <- fit_meta_regression_60(~ natal_value_z,two_step_data_60,cluster = two_step_data_60$territory_id)


#------------------------------------------------------------------------------- STEP 3: iSSF: natal effect ----
natal_coefficient_60 <- unname(model_B_60$beta["natal_value_z"])
natal_se_60 <- unname(model_B_60$se["natal_value_z"])
mean_slope_60 <- unname(model_A_60$beta["(Intercept)"])

natal_effect_60 <- tibble::tibble(
  natal_descriptor = natal_descriptor_60,
  n_territories = dplyr::n_distinct(two_step_data_60$territory_id),
  optimizer_ok = model_A_60$convergence == 0 & model_B_60$convergence == 0,
  mean_slope = mean_slope_60,
  natal_coefficient = natal_coefficient_60,
  natal_CI_95 = sprintf("[%.3f; %.3f]",natal_coefficient_60 - confidence_critical_value_60 * natal_se_60,natal_coefficient_60 + confidence_critical_value_60 * natal_se_60),
  natal_coefficient_percent_of_mean_slope = 100 * natal_coefficient_60 / abs(mean_slope_60),
  variance_explained_percent = 100 * (1 - sum(model_B_60$tau2) / sum(model_A_60$tau2)),
  delta_AIC_B_minus_A = model_B_60$AIC - model_A_60$AIC)


print(natal_effect_60,width = Inf)


#------------------------------------------------------------------------------- STEP 4: iSSF: size of the natal effect ----
#' observed_effect_percent = change in avoidance between an eagle from a lightly built natal territory (q10)
#' and one from a heavily built territory (q90), divided by the population mean slope and multiplied by 100.
#' observed_effect_CI_95 = the 95 % confidence interval on that percentage scale.
#' smallest_detectable_effect = the natal coefficient that, given the observed precision, would
#' have been declared significant in 80 % of repeated studies.
#' smallest_detectable_effect_percent = the same quantity on the percentage scale above. 
natal_contrast_60 <- stats::quantile(dplyr::distinct(two_step_data_60,territory_id,natal_value_z)$natal_value_z,natal_contrast_probabilities_60,names = FALSE)
contrast_width_60 <- diff(natal_contrast_60)
to_percent_60 <- \(coefficient) 100 * coefficient * contrast_width_60 / abs(mean_slope_60)
confidence_bounds_60 <- natal_coefficient_60 + c(-1,1) * confidence_critical_value_60 * natal_se_60
smallest_detectable_60 <- (confidence_critical_value_60 + stats::qnorm(power_level_60)) * natal_se_60

natal_effect_size_60 <- tibble::tibble(
  natal_coefficient = natal_coefficient_60,
  natal_se = natal_se_60,
  natal_CI_95 = sprintf("[%.3f; %.3f]",confidence_bounds_60[1],confidence_bounds_60[2]),
  observed_effect_percent = to_percent_60(natal_coefficient_60),
  observed_effect_CI_95 = sprintf("[%.1f; %.1f]",to_percent_60(confidence_bounds_60[1]),to_percent_60(confidence_bounds_60[2])),
  smallest_detectable_effect = smallest_detectable_60,
  smallest_detectable_effect_percent = to_percent_60(smallest_detectable_60))

print(natal_effect_size_60,width = Inf)

