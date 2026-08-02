#'-------------------------------------------------------------------------------
#' Title: iSSF fitting ----
#' Authors : Louise Faure
#' Date : 28.07.26
#' Info : this script follow the Data_processing_annotation.R script where I 
#' generated random step based on a gamma and uniform distribution for step lenght
#' and turning angle, and extracted the values below each data point.  
#' Purpose : 
#' (1) fit an iSSF for all individuals 
#' (2) validate the model through RMSE
#' (3) plot coefficient estimates and individuals coefficient
#' (4) predict interaction form 
#' ------------------------------------------------------------------------------

library(tidyverse)
library(glmmTMB)
library(patchwork) #patching up interaction plots
library(oce) #color palette for interaction plots
library(corrr)
library(terra)
library(sf)
library(stars) #use for st_rasterize
library(move2)
library(mapview)
library(ggnewscale)
library(performance)

wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")

#colors
clr <- oce::oceColorsPalette(100)[9] #was [2] before
clr_light <- oce::oceColorsPalette(100)[10]
clr2 <- oce::oceColorsPalette(100)[80]

#open data. data was prepared in 03_energy_landscape_modeling_method1.R. 
data <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/issf_generated_observed_location_annotated(2).rds") %>%
  mutate(animal_ID = as.numeric(as.factor(individual.local.identifier)), #animal ID and stratum ID should be numeric
         stratum_ID = as.numeric(as.factor(stratum)))

#look at correlation
data %>% 
  dplyr::select(c("ruggedness_100m", "step_length_km", "weeks_since_emig", "elevation_100m")) %>% 
  correlate() 

# term             ruggedness_100m step_length_km weeks_since_emig elevation_100m
# ruggedness_100m         NA           -0.0558          -0.00270          0.178 
# step_length_m           -0.0558      NA                0.000641        -0.0545
# weeks_since_emig        -0.00270      0.000641        NA                0.145 
# elevation_100m           0.178       -0.0545           0.145           NA  

### calculate total hours of flight for the manuscript
landing_sampling_summary_60 <- data %>%
  dplyr::filter(used == 1L) %>%
  dplyr::summarise(
    n_landing_steps = dplyr::n_distinct(stratum),
    total_sampled_hours = n_landing_steps
  )

print(landing_sampling_summary_60) # 5777 hours


### standardize environmental covariates (hfi_mean_1000m, elevation_100m, ruggedness_100m, weeks_since_emig)
variables_to_standardize_60 <- c(
  "hfi_mean_1000m",
  "elevation_100m",
  "ruggedness_100m",
  "weeks_since_emig",
  "step_length_km"
)

required_model_variables_60 <- c(
  "used",
  "animal_ID",
  "stratum_ID",
  variables_to_standardize_60
)

missing_model_variables_60 <- setdiff(
  required_model_variables_60,
  names(data)
)

if (length(missing_model_variables_60) > 0L) {
  stop(
    "Missing variables: ",
    paste(missing_model_variables_60, collapse = ", ")
  )
}


# Retain complete and finite observations
data_model_60 <- data %>%
  dplyr::mutate(
    used = as.integer(used),
    animal_ID = factor(animal_ID),
    stratum_ID = factor(stratum_ID)
  ) %>%
  dplyr::filter(
    !is.na(used),
    !is.na(animal_ID),
    !is.na(stratum_ID),
    dplyr::if_all(
      dplyr::all_of(variables_to_standardize_60),
      ~ !is.na(.x) & is.finite(.x)
    )
  ) %>%
  dplyr::group_by(stratum_ID) %>%
  dplyr::filter(
    sum(used == 1L) == 1L,
    sum(used == 0L) >= 1L
  ) %>%
  dplyr::ungroup()


# Calculate and retain standardization parameters
standardization_parameters_60 <- tibble::tibble(
  variable = variables_to_standardize_60,
  
  center = vapply(
    variables_to_standardize_60,
    function(variable) {
      mean(data_model_60[[variable]])
    },
    numeric(1)
  ),
  
  scale = vapply(
    variables_to_standardize_60,
    function(variable) {
      stats::sd(data_model_60[[variable]])
    },
    numeric(1)
  )
)

if (
  any(
    !is.finite(standardization_parameters_60$scale) |
    standardization_parameters_60$scale <= 0
  )
) {
  stop(
    "At least one variable has a non-finite or zero standard deviation."
  )
}

print(standardization_parameters_60)


# Create standardized variables
data_model_60 <- data_model_60 %>%
  dplyr::mutate(
    hfi_mean_1000m_z =
      (
        hfi_mean_1000m -
          standardization_parameters_60$center[
            standardization_parameters_60$variable ==
              "hfi_mean_1000m"
          ]
      ) /
      standardization_parameters_60$scale[
        standardization_parameters_60$variable ==
          "hfi_mean_1000m"
      ],
    
    elevation_100m_z =
      (
        elevation_100m -
          standardization_parameters_60$center[
            standardization_parameters_60$variable ==
              "elevation_100m"
          ]
      ) /
      standardization_parameters_60$scale[
        standardization_parameters_60$variable ==
          "elevation_100m"
      ],
    
    ruggedness_100m_z =
      (
        ruggedness_100m -
          standardization_parameters_60$center[
            standardization_parameters_60$variable ==
              "ruggedness_100m"
          ]
      ) /
      standardization_parameters_60$scale[
        standardization_parameters_60$variable ==
          "ruggedness_100m"
      ],
    
    weeks_since_emig_z =
      (
        weeks_since_emig -
          standardization_parameters_60$center[
            standardization_parameters_60$variable ==
              "weeks_since_emig"
          ]
      ) /
      standardization_parameters_60$scale[
        standardization_parameters_60$variable ==
          "weeks_since_emig"
      ],
    
    step_length_z =
      (
        step_length_km -
          standardization_parameters_60$center[
            standardization_parameters_60$variable ==
              "step_length_km"
          ]
      ) /
      standardization_parameters_60$scale[
        standardization_parameters_60$variable ==
          "step_length_km"
      ]
  )

# test 1 
# STEP 1: fit the iSSF -----------------------------------------------------------

TMB_struc <- glmmTMB::glmmTMB(
  used ~ -1 +
    
    # HFI selection and its variation with step length and age
    hfi_mean_1000m_z *
    step_length_z *
    weeks_since_emig_z +
    
    # Elevation selection and its variation with step length and age
    elevation_100m_z *
    step_length_z *
    weeks_since_emig_z +
    
    # Ruggedness selection and its variation with step length and age
    ruggedness_100m_z *
    step_length_z *
    weeks_since_emig_z +
    
    # Choice-set intercept
    (1 | stratum_ID) +
    
    # Individual variation only in the HFI response
    (0 + hfi_mean_1000m_z | animal_ID),
  
  family = poisson(link = "log"),
  data = data_model_60,
  doFit = FALSE,
  
  # Fix the stratum-intercept SD and estimate the HFI random-slope SD
  map = list(
    theta = factor(c(NA, 1L))
  ),
  
  # theta[1]: fixed stratum-intercept SD
  # theta[2]: initial HFI random-slope SD
  start = list(
    theta = c(log(1e3), 0)
  )
)


# Check that the model contains exactly two random-effect parameters:
# 1. stratum intercept SD;
# 2. individual HFI-slope SD.
stopifnot(
  length(TMB_struc$parameters$theta) == 2L
)


# Fit the model
TMB_M <- glmmTMB::fitTMB(TMB_struc)

summary(TMB_M)

saveRDS(
  TMB_M,
  file = "TMB_model_HFI_elevation_ruggedness.rds"
)




# # STEP 1: run the model ------------------------------------------------------------------ 
# #this is based on Muff et al:
# #https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40&isAllowed=y#glmmtmb-1
# 
# TMB_struc <- glmmTMB(used ~ -1 + TRI_100_z * step_length_z * weeks_since_emig_z + 
#                        ridge_100_z * step_length_z * weeks_since_emig_z + (1|stratum_ID) + 
#                        (0 + ridge_100_z | animal_ID) + 
#                        (0 + TRI_100_z | animal_ID), 
#                      family = poisson, data = data, doFit = FALSE,
#                      #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
#                      map = list(theta = factor(c(NA, 1:2))), #2 is the n of random slopes
#                      #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
#                      start = list(theta = c(log(1e3), 0, 0))) #add a 0 for each random slope. in this case, 2
# 
# 
# TMB_M <- glmmTMB:::fitTMB(TMB_struc)
# summary(TMB_M)
# 
# saveRDS(TMB_M, file = "TMB_model.rds")
# 
# #extract coefficient estimates and confidence intervals
# confint(TMB_M)
# 
# #extract individual-specific random effects:
# ranef(TMB_M)[[1]]$animal_ID

# STEP 2: model validation ------------------------------------------------------------------ 

#calculate the RMSE
performance_rmse(TMB_M)

# STEP 3: PLOT coefficient estimates ------------------------------------------------------------------ 

graph <- confint(TMB_M) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Factor") %>% 
  filter(!(Factor %in% c("weeks_since_emig_z", "Std.Dev.hfi_mean_1000m_z|animal_ID")))

colnames(graph)[c(2,3)] <- c("Lower", "Upper") 

labels <- rev(c("hfi_mean_1000m_z", "Step length", "Distance to ridge", "TRI: Step length", "TRI: Week",
                "Step length: Week", "Step length: Distance to ridge", "Distance to ridge: Week",
                "TRI: Step length: Week", "Distance to ridge: Step length: Week"))

VarOrder <- rev(unique(graph$Factor))
graph$Factor <- factor(graph$Factor, levels = VarOrder)

#plot
X11(width = 7, height = 2) 

coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", linewidth = 0.5) +
  geom_point(color = clr, size = 1.7)  +
  geom_linerange(aes(xmin = Lower, xmax = Upper),color = clr, linewidth = .7) +
  labs(x = "Estimate", y = "") +
  scale_y_discrete(labels = labels) +
  xlim(-.71, .25) +
  theme_minimal() +
  theme(text = element_text(size = 8), #font size should be between 6-8
        axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
        axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6)))

ggsave(coefs, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/coeffs.pdf", 
       width = 7, height = 2, dpi = 300)



# test 2 
# Prepare movement variables ----------------------------------------------------

data_model_60 <- data_model_60 %>%
  dplyr::mutate(
    cos_turning_angle = cos(turning_angle_rad)
  )

stopifnot(
  all(is.finite(data_model_60$step_length_z)),
  all(is.finite(data_model_60$cos_turning_angle)),
  all(is.finite(data_model_60$weeks_since_emig_z))
)


# STEP 1: fit the iSSF -----------------------------------------------------------

TMB_struc <- glmmTMB::glmmTMB(
  used ~ -1 +
    
    # Mean HFI selection at mean age
    hfi_mean_1000m_z +
    
    # Ontogenetic change in HFI selection
    hfi_mean_1000m_z:weeks_since_emig_z +
    
    # Environmental control variables
    elevation_100m_z +
    ruggedness_100m_z +
    
    # Movement terms without interactions with age
    step_length_z +
    cos_turning_angle +
    
    # Choice-set intercept
    (1 | stratum_ID) +
    
    # Individual variation only in the HFI response
    (0 + hfi_mean_1000m_z | animal_ID),
  
  family = poisson(link = "log"),
  data = data_model_60,
  doFit = FALSE,
  
  # theta[1] = stratum-intercept SD, fixed at 1000
  # theta[2] = individual HFI random-slope SD, freely estimated
  map = list(
    theta = factor(c(NA, 1L))
  ),
  
  start = list(
    theta = c(log(1e3), 0)
  )
)


# Check the expected random-effect structure:
# 1. stratum-intercept SD;
# 2. individual HFI-slope SD.
stopifnot(
  length(TMB_struc$parameters$theta) == 2L
)


# Fit the model
TMB_M <- glmmTMB::fitTMB(TMB_struc)

summary(TMB_M)


# Model controls
stopifnot(
  TMB_M$fit$convergence == 0L,
  isTRUE(TMB_M$sdr$pdHess)
)

glmmTMB::diagnose(TMB_M)


# Save model
saveRDS(
  TMB_M,
  file = "TMB_model_HFI_age_elevation_ruggedness.rds"
)


# STEP 3: extract fixed-effect coefficients -------------------------------------

fixed_effect_estimates_60 <- glmmTMB::fixef(TMB_M)$cond

fixed_effect_confint_60 <- stats::confint(
  TMB_M,
  parm = "beta_",
  method = "wald"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term")

names(fixed_effect_confint_60)[2:3] <- c(
  "confidence_low",
  "confidence_high"
)

fixed_effect_table_60 <- tibble::tibble(
  term = names(fixed_effect_estimates_60),
  estimate = unname(fixed_effect_estimates_60)
) %>%
  dplyr::left_join(
    fixed_effect_confint_60,
    by = "term"
  )

print(fixed_effect_table_60)




# fin test 2

# STEP 4: PLOT individual-specific coefficients ------------------------------------------------------------------ 

rnd <- ranef(TMB_M) %>% 
  as.data.frame() %>% 
  filter(grpvar == "animal_ID") %>% 
  mutate(ID = rep(levels(as.factor(data$ind1)),2),
         variable = c(rep("Distance_to_ridge", 55), rep("TRI", 55)),
         Lower = condval - condsd,
         Upper = condval + condsd)  #assign individuals' names

order <- rnd %>% 
  filter(variable == "Distance_to_ridge") %>% 
  arrange(desc(condval)) %>% 
  pull(ID)

rnd$ID <- factor(rnd$ID, levels = order)

#cols <- c(TRI = "lightcoral",
#          Distance_to_ridge = "cornflowerblue")

cols <- c(TRI = clr2,
          Distance_to_ridge = clr)

X11(width = 7, height = 9)
coefs_inds <- ggplot(rnd, aes(x = condval, y = ID, color = variable)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", linewidth = 0.5) +
  geom_point(size = 2, position = position_dodge(width = .7))  +
  geom_linerange(aes(xmin = Lower, xmax = Upper), size = 0.8, position = position_dodge(width = .7)) +
  scale_color_manual(values = cols) + 
  scale_y_discrete(labels = order) +
  labs(x = "Difference from fixed effect estimate", y = "") +
  theme_minimal() +
  theme(text = element_text(size = 8), #font size should be between 6-8
        axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
        axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6)))

ggsave(plot = coefs_inds, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/ind_coefs.pdf", 
       width = 7, height = 9, dpi = 300)


# STEP 5: predictions for interaction terms + PLOT ------------------------------------------------------------------ 

#new dataset created in 03_energy_landscape_modeling_method1.r ... this file was updated to have stratum_ID and animal_ID as numeric variables
#accoring to the predict.glmmTMB help file: "To compute population-level predictions for a given grouping variable 
#(i.e., setting all random effects for that grouping variable to zero), set the grouping variable values to NA."
new_data <- readRDS("new_data_only_ssf_preds.rds") %>% 
  mutate(stratum_ID = NA)
#       animal_ID = NA) #set these to NA for population-level predictions


#predict using the model
#preds <- predict(TMB_M, newdata = new_data, type = "link", se.fit = T)
preds <- predict(TMB_M, newdata = new_data, type = "link")

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>% 
  mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob

#prepare for plotting
y_axis_var <- c("TRI_100", "ridge_100")
x_axis_var <- "weeks_since_emig"
labels <- data.frame(vars = y_axis_var, 
                     label = c("Topographic Ruggedness Index", "Distance to ridge line (km)"))

for (i in y_axis_var){
  
  label <- labels %>% 
    filter(vars == i) %>% 
    pull(label)
  
  #interaction to be plotted
  interaction_term <- paste0("wk_", i)
  pred_r <- preds_pr %>% 
    filter(interaction == interaction_term) %>%  #only keep the rows that contain data for this interaction term
    dplyr::select(c(which(names(.) %in% c(x_axis_var, i)), "probs")) %>% 
    terra::rast(type = "xyz") %>%
    #focal(w = 3, fun = mean, na.policy = "all", na.rm = T) %>%
    as.data.frame(xy = T) #%>% 
  #rename(probs = focal_mean)
  
  #X11(width = 6.9, height = 3.5)
  if(i == "ridge_100"){
    pred_p <- pred_r %>% 
      ggplot() +
      geom_tile(aes(x = x, y = y, fill = probs)) +
      scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                           na.value = "white", name = "Flyability")+
      guides(fill = guide_colourbar(title.vjust = .95)) + #the legend title needs to move up a bit
      labs(x = "Week since dispersal", y = label) + #add label for GRC plot
      theme_minimal() +
      theme(plot.margin = margin(0, 15, 0, 0, "pt"),
            legend.direction="horizontal",
            legend.position = "bottom",
            legend.key.width=unit(.7,"cm"),
            legend.key.height=unit(.25,"cm"),
            text = element_text(size = 8), #font size should be between 6-8
            axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
            axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6))) 
  } else{
    pred_p <- pred_r %>% 
      ggplot() +
      geom_tile(aes(x = x, y = y, fill = probs)) +
      scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                           na.value = "white", name = "Flyability")+
      guides(fill = guide_colourbar(title.vjust = .95)) + #the legend title needs to move up a bit
      labs(x = "", y = label) + #add label for GRC plot
      theme_minimal() +
      theme(plot.margin = margin(0, 15, 0, 0, "pt"),
            legend.direction="horizontal",
            legend.position = "bottom",
            legend.key.width=unit(.7,"cm"),
            legend.key.height=unit(.25,"cm"),
            text = element_text(size = 8), #font size should be between 6-8
            axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
            axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6))) 
  }
  assign(paste0(i, "_p"), pred_p)
  
}

#plot all interaction plots together
X11(width = 4.5, height = 4.8)
combined <- TRI_100_p + ridge_100_p & theme(legend.position = "bottom")
p <- combined + plot_layout( guides = "collect", nrow = 2)

ggsave(p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/interactions_flyability.pdf", 
       width = 4.5, height = 4.8, dpi = 400)

# STEP 6: predictions for the Alps + PLOTS ------------------------------------------------------------------ 

#read in the ssf model
TMB_M <- readRDS( "TMB_model.rds")

#Alpine df prepared in 03_04_clogit_workflow.R
topo_df <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/topo_df_100_LF.rds")

#prepare dist to ridge to be used as the base layer for plotting
ridge_100 <- rast("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/ridge_100_LF.tif") %>% 
  #aggregate(fact = 5) %>% #200*200 m
  as.data.frame(xy = T) %>% 
  drop_na(distance_to_ridge_line_mask)
ridge_100[ridge_100$distance_to_ridge_line_mask > 5000, "distance_to_ridge_line_mask"] <- NA #do this for the plot to look nicer... NA values will be white

#saveRDS(ridge_100, file = "ridge_100_LF_df.rds")

#ridge_100 <- readRDS("ridge_100_LF_df.rds")

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/alpine_preds_7/")

wks_ls <- split(data, data$weeks_since_emig) #density maps were only made for four timestamps

gc()

(start <- Sys.time())
for(x in wks_ls){
  
  #extract week number
  one_week <- x %>% 
    distinct(weeks_since_emig) %>% 
    pull(weeks_since_emig)
  
  #create week ID to be used for naming files and plots
  week_i <- one_week %>% 
    str_pad(3,"left","0")
  
  #generate a new dataset
  topo_df <- topo_df %>%
    mutate(step_length = mean(x$step_length),
           step_length_z = 0,
           weeks_since_emig = one_week,
           weeks_since_emig_z = (one_week - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale'),
           stratum_ID = NA,
           animal_ID = sample(x$animal_ID, nrow(topo_df), replace = T)) #set the grouping variables to NA as per glmm.TMB recommendation
  
  #predict using the model
  #(b <- Sys.time())
  topo_df$preds <- predict(TMB_M, newdata = topo_df, type = "link")
  #Sys.time() - b #24 min for one week
  
  topo_df$probs <- gtools::inv.logit(topo_df$preds)
  
  gc()
  
  
  #calculate suitable areas
  #area_.7 <- topo_df %>% 
  #  filter(probs >= .7) %>% 
  #  summarize(pixels = n()) %>% #count the 
  #  mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
  #         area_km2 = round(area_m2/1e6,3),
  #         week_since_dispersal = week_i)
  
  #save area as a file
  #save(area_.7, file = paste0("area_alps_wk_", week_i, ".r"))
  
  
  #plot the raw map for a few weeks
  #  t <- topo_df %>% 
  #          ggplot() +
  #          geom_tile(aes(x = location.long, y = location.lat, fill = probs)) +
  #          scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
  #                               na.value = "white", name = "Intensity of use") +
  #          labs(x = "", y = "", title = paste0("Week ", one_week, " since dispersal")) +
  #          theme_void()
  
  #  ggsave(plot = t, filename = paste0("alps_wk_", week_i, ".tiff"), device = "tiff", width = 7, height = 4.5, dpi = 400)
  
  #density plot
  p <- ggplot() +  
    geom_tile(data = ridge_100, aes(x = x, y = y, fill = scale(distance_to_ridge_line_mask))) +
    scale_fill_gradientn(colors = grey.colors(100), guide = "none", na.value = "white") +
    new_scale_fill() +
    stat_density_2d(data = topo_df %>% filter(probs >= .7) %>% dplyr::select("location.lat", "location.long", "probs"), 
                    aes(x = location.long, y = location.lat, fill = after_stat(level)), geom = "polygon") +
    scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(100)[51:100], alpha = .2)) +
    labs(title = paste0("Week ", one_week, " since dispersal"), x = "", y = "") +
    theme_void()
  
  #dev.off()
  ggsave(plot = p, filename = paste0(week_i, "_alpine_pred.tiff"), device = "tiff", width = 7, height = 4.5, dpi = 400) #3min
  rm(p)
  gc()
  
  print(paste0("week ", week_i, " of 156 done!"))
  
}

Sys.time() - start #30 min per week
