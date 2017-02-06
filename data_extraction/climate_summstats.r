# Get the climate records and locations of where each species' body mass records come from in order to look at the range of where the specimens are from and the range of climatic conditions in there.

bird_clim <- read.csv(file.path(fp, 'bird_clim.csv'))

vnsis <- cbind(vnsis, bird_clim)

# Calculate climate stats for specimen locations, and the space covered by the specimen locations

birdrange_summstats <- vnsis %>% group_by(binomial) %>% 
  summarize(spatial_cv_temp = sd(bio1 + 273, na.rm = T)/mean(bio1 + 273, na.rm = T),
            spatial_cv_precip = sd(bio12, na.rm = T)/mean(bio12, na.rm = T),
            seasonal_var_temp = mean(bio4, na.rm = T),
            seasonal_var_precip = mean(bio15, na.rm = T),
            interannual_var_temp = mean(temp_interannual, na.rm = T),
            interannual_var_precip = mean(prec_interannual, na.rm = T),
            spatial_spread = mean(dist(cbind(decimallatitude, decimallongitude)), na.rm = T))

write.csv(birdrange_summstats, file.path(fp, 'bird_clim_summstats.csv'), row.names = FALSE)
