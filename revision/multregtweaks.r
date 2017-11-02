# Try with different predator types.

pred <- c('Invertebrate','VertFishScav')
herb <- c('PlantSeed', 'FruiNect')
omni <- c('Omnivore')

multregdat <- sister_data %>% mutate(spatial_temp = nontropical_spatial_cv_temp - tropical_spatial_cv_temp,
                                     interann_temp = nontropical_interannual_var_temp - tropical_interannual_var_temp,
                                     spatial_precip = nontropical_spatial_cv_precip - tropical_spatial_cv_precip,
                                     interann_precip = nontropical_interannual_var_precip - tropical_interannual_var_precip,
                                     rangesize = log10(nontropical_range_size) - log10(tropical_range_size),
                                     total_richness = nontropical_total_richness - tropical_total_richness,
                                     congener_richness = nontropical_congener_richness - tropical_congener_richness,
                                     n_pop = nontropical_nsubpop - tropical_nsubpop,
                                     elev_cv = nontropical_elevation_cv - tropical_elevation_cv,
                                     seasonal_temp = nontropical_seasonal_var_temp - tropical_seasonal_var_temp,
                                     seasonal_precip = nontropical_seasonal_var_precip, tropical_seasonal_var_precip,
                                     area = log10(nontropical_area + 1) - log10(tropical_area + 1),
                                     meanlogbodymass = apply(cbind(log10(nontropical_BodyMass.Value), log10(tropical_BodyMass.Value)), 1, mean),
                                     predator = nontropical_Diet.5Cat %in% pred & tropical_Diet.5Cat %in% pred,
                                     herbivore = nontropical_Diet.5Cat %in% herb & tropical_Diet.5Cat %in% herb,
                                     migrant = nontropical_migrant_status %in% c('partial','obligate') | tropical_migrant_status %in% c('partial','obligate'),
                                     dist_phylo = dist_phy,
                                     dcv = nontropical_cv_logmass - tropical_cv_logmass) %>%
  dplyr::select(dcv, dlat, dist_phylo, dist_slc, spatial_temp, interann_temp, spatial_precip, interann_precip, rangesize, total_richness, congener_richness, n_pop, elev_cv, seasonal_temp, seasonal_precip, area, meanlogbodymass, predator, herbivore, migrant) 

multreg_complete <- filter(multregdat, complete.cases(multregdat)) %>% as.data.frame

multreg_complete <- mutate(multreg_complete, omnivore = !predator & !herbivore) %>% dplyr::select(-herbivore)

# Forward stepwise model selection using AICc
# Must load a modified version of a function from the MASS package in order to perform automated forward selection using AICc.

require(MASS)
require(MuMIn)

min_model <- lm(dcv ~ 1, data = multreg_complete)
full_model <- lm(dcv ~ ., data = multreg_complete)
full_model_formula <- formula(full_model)

source('stepaicc.r')
aicc_step_result <- stepAICc(min_model, direction = 'forward', scope = full_model_formula)

# Standardized coefficients
best_model_data <- multreg_complete %>%
  transmute(dcv = dcv,
            migrant = migrant,
            herbivore = herbivore,
            spatial_temp = spatial_temp/sd(spatial_temp),
            interann_temp = interann_temp/sd(interann_temp),
            seasonal_precip = seasonal_precip/sd(seasonal_precip))

best_model_data <- multreg_complete %>%
  transmute(dcv = dcv,
            predator = predator,
            spatial_temp = spatial_temp/sd(spatial_temp),
            interann_temp = interann_temp/sd(interann_temp),
            seasonal_precip = seasonal_precip/sd(seasonal_precip))

best_model_std <- lm(dcv ~ ., data = best_model_data)

# Correlations
cor(multreg_complete[,c('seasonal_precip','interann_temp','spatial_temp')])

# Polychoric
library(polycor)
hetcor(multreg_complete[,c('seasonal_precip','interann_temp','spatial_temp','migrant','herbivore')])
