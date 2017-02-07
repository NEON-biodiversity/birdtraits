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

# Get a more accurate measure of spatial spread by converting the coordinates into a global CRS.
library(sp)
grid_crs <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
bird_clim_sp <- spTransform(SpatialPoints(coords=vnsis[!is.na(vnsis$decimallatitude),c('decimallongitude','decimallatitude')], proj4string = CRS('+proj=longlat +datum=WGS84')), CRSobj = grid_crs)

birdrange_summstats <- vnsis %>%
  filter(!is.na(decimallatitude)) %>%
  cbind(data.frame(latm = bird_clim_sp@coords[,2], longm = bird_clim_sp@coords[,1])) %>%
  group_by(binomial) %>%
  summarize(spatial_cv_temp = sd(bio1 + 273, na.rm = T)/mean(bio1 + 273, na.rm = T),
            spatial_cv_precip = sd(bio12, na.rm = T)/mean(bio12, na.rm = T),
            seasonal_var_temp = mean(bio4, na.rm = T),
            seasonal_var_precip = mean(bio15, na.rm = T),
            interannual_var_temp = mean(temp_interannual, na.rm = T),
            interannual_var_precip = mean(prec_interannual, na.rm = T),
            spatial_spread = mean(dist(cbind(latm, longm)), na.rm = T))

write.csv(birdrange_summstats, file.path(fp, 'bird_clim_summstats.csv'), row.names = FALSE)


# Attempt to get the area of a convex hull in km2

library(sp)

chull_area <- function(dat) {
  hpts <- chull(x = dat$longm, y = dat$latm)
  hpts <- c(hpts, hpts[1]) # This line may not be necessary.
  hcoords <- dat[hpts, c('longm','latm')]
  hpoly <- Polygon(hcoords, hole = F)
  return(area = as.numeric(hpoly@area)/1e6)
}

birdrange_hullareas <- vnsis %>%
  filter(!is.na(decimallatitude)) %>%
  cbind(data.frame(latm = bird_clim_sp@coords[,2], longm = bird_clim_sp@coords[,1])) %>%
  group_by(binomial) %>%
  do(data.frame(area = chull_area(.)))

write.csv(cbind(birdrange_summstats, area = birdrange_hullareas$area), file.path(fp, 'bird_clim_summstats.csv'), row.names = FALSE)
