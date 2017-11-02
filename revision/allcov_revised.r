# Script for revision: get covariates for the corrected coordinates.

#### COVARIATE SET 1

foraging <- read.delim(file.path(fp, 'BirdFuncDat.txt'), stringsAsFactors = FALSE)
foraging$binomial <- gsub(' ', '_', foraging$Scientific)
taxgrps <- with(foraging, data.frame(nontropical_binomial = binomial, order = IOCOrder, family = BLFamilyLatin, stringsAsFactors = FALSE))
sister_data <- left_join(sister_data, taxgrps)

# Is it passerine or not?
sister_data$songbird <- sister_data$order == 'Passeriformes'

#### COVARIATE SETS 2-6

bird_clim <- read.csv(file.path(fprev, 'vertnet_clim.csv'), stringsAsFactors = FALSE)

# Convert coordinates to a grid based CRS.
library(sp)
grid_crs <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
bird_clim_sp <- spTransform(SpatialPoints(coords=bird_clim[!is.na(bird_clim$decimallatitude),c('decimallongitude','decimallatitude')], proj4string = CRS('+proj=longlat +datum=WGS84')), CRSobj = grid_crs)

bird_clim$latm <- NA
bird_clim$longm <- NA
bird_clim$latm[!is.na(bird_clim$decimallatitude)] <- bird_clim_sp@coords[,2]
bird_clim$longm[!is.na(bird_clim$decimallongitude)] <- bird_clim_sp@coords[,1]

birdrange_summstats <- bird_clim %>%
  filter(!is.na(decimallatitude)) %>%
  group_by(binomial) %>%
  summarize(spatial_cv_temp = sd(bio1 + 273, na.rm = T)/mean(bio1 + 273, na.rm = T),
            spatial_cv_precip = sd(bio12, na.rm = T)/mean(bio12, na.rm = T),
            seasonal_var_temp = mean(bio4, na.rm = T),
            seasonal_var_precip = mean(bio15, na.rm = T),
            interannual_var_temp = mean(temp_interannual, na.rm = T),
            interannual_var_precip = mean(prec_interannual, na.rm = T),
            spatial_spread = mean(dist(cbind(latm, longm)), na.rm = T))

chull_area <- function(dat) {
  hpts <- chull(x = dat$longm, y = dat$latm)
  hpts <- c(hpts, hpts[1]) # This line may not be necessary.
  hcoords <- dat[hpts, c('longm','latm')]
  hpoly <- Polygon(hcoords, hole = F)
  return(area = as.numeric(hpoly@area)/1e6)
}

birdrange_hullareas <- bird_clim %>%
  filter(!is.na(decimallatitude)) %>%
  group_by(binomial) %>%
  do(data.frame(area = chull_area(.)))

#### COVARIATES 7-8

# Use BOTW polygons, matched with Vertnet names.

botw_cents <- mutate(botw_cents, binomial = gsub(pattern = '\\ ', '_', SCINAME))

for (i in 1:nrow(goodnames)) {
  idx <- which(botw_cents$binomial == goodnames$botw_name[i])
  botw_cents$binomial[idx] <- goodnames$vertnet_name[i] 
}


# Migrant status and range size.
rangedat <- function(x) {
  migrant_status <- 'none'
  if (any(x$SEASONAL %in% 2:3)) migrant_status <- 'partial'
  if (!any(x$SEASONAL == 1)) migrant_status <- 'obligate'
  data.frame(migrant_status = migrant_status, range_size = sum(x$Shape_Area[x$SEASONAL %in% 1:2]))
}

# Sum breeding and year-round ranges 
botw_ranges <- botw_cents %>%
  filter(PRESENCE %in% 1:3, ORIGIN %in% 1) %>% # Only where species is known present.
  group_by(binomial) %>%
  do(rangedat(.))

#### COVARIATES 9-10
foraging <- read.delim(file.path(fp, 'BirdFuncDat.txt'), stringsAsFactors = FALSE)
foraging$binomial <- gsub(' ', '_', foraging$Scientific)
lifehist <- read.csv(file.path(fp, 'Amniote_Database_Aug_2015.csv'), stringsAsFactors = FALSE)
lifehist$binomial <- paste(lifehist$genus, lifehist$species, sep = '_')

c(tropical_species, nontropical_species) %in% foraging$binomial
c(tropical_species, nontropical_species) %in% lifehist$binomial
dput(c(tropical_species, nontropical_species)[!c(tropical_species, nontropical_species) %in% lifehist$binomial])

wrong_names_lh <- c("Aimophila_sumichrasti", "Diuca_speculifera", "Euplectes_psammocromius", 
                    "Hylocharis_leucotis", "Parula_pitiayumi", "Pipilo_albicollis", 
                    "Serinus_sulphuratus", "Sterna_sandvicensis", "Aimophila_carpalis", 
                    "Hylocharis_xantusii", "Parula_americana", "Pipilo_fuscus", "Serinus_flaviventris", 
                    "Sterna_elegans")

corr_names_lh <- c('Peucaea_sumichrasti', 'Chionodacryon_speculiferum', 'Euplectes_psammacromius',
                   'Basilinna_leucotis', 'Setophaga_pitiayumi', 'Melozone_albicollis',
                   'Crithagra_sulphurata', 'Thalasseus_sandvicensis', 'Peucaea_carpalis',
                   'Basilinna_xantusii', 'Setophaga_americana', 'Melozone_fusca', 'Crithagra_flaviventris',
                   'Thalasseus_elegans')


sister_names_lh <- c(tropical_species, nontropical_species)
sister_names_lh[!c(tropical_species, nontropical_species) %in% lifehist$binomial] <- corr_names_lh

foraging <- foraging[,c("binomial", "Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", 
                        "Diet.Vunk", "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", 
                        "Diet.PlantO", "Diet.5Cat", "ForStrat.watbelowsurf", "ForStrat.wataroundsurf", 
                        "ForStrat.ground", "ForStrat.understory", "ForStrat.midhigh", 
                        "ForStrat.canopy", "ForStrat.aerial", "PelagicSpecialist", "Nocturnal", 
                        "BodyMass.Value")]
foraging_sisters <- left_join(data.frame(binomial =c(tropical_species,nontropical_species)), foraging)

lifehist <- lifehist[,c("binomial", "female_maturity_d", "litter_or_clutch_size_n", 
                        "litters_or_clutches_per_y", "adult_body_mass_g", "maximum_longevity_y", 
                        "gestation_d", "weaning_d", "birth_or_hatching_weight_g", "weaning_weight_g", 
                        "egg_mass_g", "incubation_d", "fledging_age_d", "longevity_y", 
                        "male_maturity_d", "inter_litter_or_interbirth_interval_y", "female_body_mass_g", 
                        "male_body_mass_g", "no_sex_body_mass_g", "egg_width_mm", "egg_length_mm", 
                        "fledging_mass_g", "adult_svl_cm", "male_svl_cm", "female_svl_cm", 
                        "birth_or_hatching_svl_cm", "female_svl_at_maturity_cm", "female_body_mass_at_maturity_g", 
                        "no_sex_svl_cm", "no_sex_maturity_d")]
lifehist_sisters <- left_join(data.frame(binomial = sister_names_lh), lifehist)

trait_sisters <- cbind(foraging_sisters, lifehist_sisters[,-1])
trait_sisters[trait_sisters==-999] <- NA

#### COVARIATE SET 11

# botw_grid was run on cluster.

fp <- '~/verts/mapoutput'
lapply(as.list(file.path(fp, dir(fp, pattern = 'Oct'))), load, .GlobalEnv)
overlap_all <- do.call('c', lapply(paste0('overlap', 1:250), get))

bird <- read.csv('~/verts/vertnet_coords.csv', stringsAsFactors = FALSE)
bird <- subset(bird, !is.na(decimallatitude) & !is.na(decimallongitude))

# Get genus for each of the data points
polygondat <- read.csv('~/verts/polygondat.csv', stringsAsFactors = FALSE)
bird$genus <- sapply(strsplit(bird$binomial, '_'), '[', 1)

# Get the richness and congeneric richness at each data point.
total_richness <- sapply(overlap_all, length)
congener_richness <- numeric(length(overlap_all))

bbin <- gsub('_', ' ', bird$binomial)

for (i in 1:length(overlap_all)) {
  spp_i <- polygondat[overlap_all[[i]], ]
  congener_richness[i] <- sum(bird$genus[i] == spp_i$genus & bbin[i] != spp_i$SCINAME)
}

richness_df <- data.frame(total_richness, congener_richness)
bird <- cbind(bird, richness_df)

write.csv(bird, '~/verts/bird_withrichness_Oct.csv', row.names = FALSE)

birdrichness <- read.csv(file.path(fprev, 'bird_withrichness_Oct.csv'), stringsAsFactors = FALSE)

birdrichness_mean <- birdrichness %>%
  group_by(binomial) %>%
  summarize(total_richness = mean(total_richness),
            congener_richness = mean(congener_richness))

#### COVARIATE SET 12

# Number of subpopulations from which each set of specimens was collected. Need to calculate.
# Round lat and long to nearest 0.5'. Anything which shares those coordinates is a subpopulation.
nsubpops <- vnsis_filter %>% 
  mutate(latround = plyr::round_any(decimallatitude, 1/120),
         lonround = plyr::round_any(decimallongitude, 1/120)) %>%
  group_by(binomial) %>%
  summarize(n = nrow(unique(cbind(latround, lonround))))

#### COVARIATE 13

bird_elev <- read.csv(file.path(fprev, 'bird_elev.csv'))
bird_clim <- left_join(bird_clim, bird_elev %>% rename(decimallatitude=lat, decimallongitude=lon))

bird_topovar <- bird_clim %>%
  mutate(elev0 = elev - min(elev, na.rm=T)) %>%
  group_by(binomial) %>%
  summarize(elevation_cv = sd(unique(elev0), na.rm = T)/mean(unique(elev0), na.rm = T))

#### COMPILE ALL COVARIATES

# also get distances in latitude, distance between the median points in great-circle distance, and phylogenetic distance.

# Difference in latitude of centroids, and distance in great circle.
deg2rad <- function(deg) return(deg*pi/180)

# Great circle distance: https://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.slc <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # Distance in km
}

# Phylogenetic distance
dist1 <- sisters$dist[match(tropical_species, sisters$sister1)]
dist2 <- sisters$dist[match(tropical_species, sisters$sister2)]
dist_df <- data.frame(tropical_binomial=tropical_species, nontropical_binomial=nontropical_species, dist_phy=pmin(dist1,dist2,na.rm=T))

dist_df <- left_join(dist_df, botw_filtered %>% 
                       rename(lat_trop = lat_centroid, lon_trop = lon_centroid, tropical_binomial = binomial) %>%
                       select(tropical_binomial, lat_trop, lon_trop))
dist_df <- left_join(dist_df, botw_filtered %>% 
                       rename(lat_nontrop = lat_centroid, lon_nontrop = lon_centroid, nontropical_binomial = binomial) %>%
                       select(nontropical_binomial, lat_nontrop, lon_nontrop))

dist_df <- mutate(dist_df, 
                  dlat = abs(lat_nontrop) - abs(lat_trop),
                  dist_slc = gcd.slc(lon_trop, lat_trop, lon_nontrop, lat_nontrop))

all_covariates <- trait_sisters %>%
  left_join(birdrange_summstats) %>%
  left_join(birdrange_hullareas) %>%
  left_join(botw_ranges) %>%
  left_join(birdrichness_mean) %>%
  left_join(nsubpops %>% rename(nsubpop = n)) %>%
  left_join(bird_topovar)

write.csv(all_covariates, file.path(fprev, 'all_covariates.csv'), row.names = FALSE)
write.csv(dist_df, file.path(fprev, 'all_distances.csv'), row.names = FALSE)

####################################
# addendum 1 Nov:
# Check correlations of some of the bio variables for a small side point in the MS.
with(bird_clim, cor(bio4, bio15, use = 'complete.obs'))
with(bird_clim, cor(bio12, bio15, use = 'complete.obs'))
