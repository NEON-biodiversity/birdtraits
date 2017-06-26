# Calculate all covariates for the bird sister pair analysis, and compile them into a data frame.

# Covariates needed:
# 1. order, family, whether or not it is passeriformes

  # 1. truesisters.r; line 318.

# 2. spatial CV of MAT and MAP in collection locations
# 3. seasonal CV of MAT and MAP in collection locations
# 4. interannual CV of MAT and MAP in collection locations
# 5. spatial spread of collection locations
# 6. area of convex hull around collection locations

  # 2-6: climate_summstats.r

# 7. whether or not it is a migrant
# 8. range size

  # 7-8. Run in getranges.r (produces botw_ranges.csv)

# 9. diet and foraging traits
# 10. life history traits

  # 9 & 10: gettraits.r

# 11. total and congeneric richness in collection locations

  # 11: must rerun botw_grid.r on cluster. Also additional code in truesisters.r

# 12. number of populations collected from

  # 12. truesisters.r; line 431

# 13. CV of elevations of collection locations

  # 13. extractsrtm.r


##################################################################################
# Load existing names.

fp <- 'C:/Users/Q/Dropbox/projects/verts'
sister_troptemp <- read.csv(file.path(fp, 'sister_troptemp_23Jun.csv'), stringsAsFactors = FALSE)

# Save coordinates to external file for extracting point level covariates.
birdcoords <- vnsis_flag2 %>% 
  filter(!outlierflag) %>%
  select(binomial, decimallatitude, decimallongitude)
write.csv(birdcoords, file = file.path(fp, 'birdcoords_23Jun.csv'), row.names = FALSE)

#### COVARIATE SET 1

foraging <- read.delim(file.path(fp, 'BirdFuncDat.txt'), stringsAsFactors = FALSE)
foraging$binomial <- gsub(' ', '_', foraging$Scientific)
taxgrps <- with(foraging, data.frame(sister1 = binomial, order = IOCOrder, family = BLFamilyLatin, stringsAsFactors = FALSE))
sister_troptemp <- left_join(sister_troptemp, taxgrps)

# Is it passerine or not?
sister_troptemp$songbird <- sister_troptemp$order == 'Passeriformes'

#### COVARIATE SETS 2-6

bird_clim <- read.csv(file.path(fp, 'bird_clim_23Jun.csv'))
birdcoords <- cbind(birdcoords, bird_clim)

library(sp)
grid_crs <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
bird_clim_sp <- spTransform(SpatialPoints(coords=birdcoords[!is.na(birdcoords$decimallatitude),c('decimallongitude','decimallatitude')], proj4string = CRS('+proj=longlat +datum=WGS84')), CRSobj = grid_crs)

birdrange_summstats <- birdcoords %>%
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

chull_area <- function(dat) {
  hpts <- chull(x = dat$longm, y = dat$latm)
  hpts <- c(hpts, hpts[1]) # This line may not be necessary.
  hcoords <- dat[hpts, c('longm','latm')]
  hpoly <- Polygon(hcoords, hole = F)
  return(area = as.numeric(hpoly@area)/1e6)
}

birdrange_hullareas <- birdcoords %>%
  filter(!is.na(decimallatitude)) %>%
  cbind(data.frame(latm = bird_clim_sp@coords[,2], longm = bird_clim_sp@coords[,1])) %>%
  group_by(binomial) %>%
  do(data.frame(area = chull_area(.)))


write.csv(cbind(birdrange_summstats, area = birdrange_hullareas$area), file.path(fp, 'bird_clim_summstats_23Jun.csv'), row.names = FALSE)

#### COVARIATE SETS 7-8

# Range sizes from the Birds of the World dataset. Calculated them in QGIS, now load into R
botw <- read.csv(file.path(fp,'botw_table.csv'), stringsAsFactors = FALSE)

# Should be able to get the migrant status out of this too.
botw <- mutate(botw, taxon = gsub(' ', '_', SCINAME))

taxa <- with(sister_troptemp, c(sister1,sister2))

botwmatch <- taxa %in% botw$taxon
#write.csv(sister_pair_long$taxon[!botwmatch], file = file.path(fp, 'botwnomatch.csv'), row.names = F)

# Reload the manually corrected name list
botwcorrect <- read.csv(file.path(fp, 'botwnomatch.csv'), stringsAsFactors = F)
for (i in 1:nrow(botwcorrect)) {
  idx <- botw$taxon == botwcorrect$correct[i]
  botw$taxon[idx] <- botwcorrect$x[i]
}

# Get the ranges.
botw_sisters <- subset(botw, taxon %in% taxa) 

rangedat <- function(x) {
  migrant_status <- 'none'
  if (any(x$SEASONAL %in% 2:3)) migrant_status <- 'partial'
  if (!any(x$SEASONAL == 1)) migrant_status <- 'obligate'
  data.frame(migrant_status = migrant_status, range_size = sum(x$Area_km2[x$SEASONAL %in% 1:2]))
}

botw_ranges <- botw_sisters %>%
  filter(PRESENCE == 1) %>% # Only where species is known present.
  group_by(taxon) %>%
  do(rangedat(.))

write.csv(botw_ranges, file = file.path(fp, 'botw_ranges_23Jun.csv'), row.names = FALSE)

#### COVARIATE SETS 9-10

sister_names2 <- with(sister_troptemp, c(sister1,sister2))

corr_names <- c('Ardenna_pacifica', 'Melozone_albicollis', 'Peucaea_sumichrasti', 'Setophaga_pitiayumi', 'Euplectes_psammacromius', 'Ardenna_bulleri', 'Melozone_fusca', 'Peucaea_carpalis', 'Setophaga_americana')

# Functional traits
foraging <- read.delim(file.path(fp, 'BirdFuncDat.txt'), stringsAsFactors = FALSE)
lifehist <- read.csv(file.path(fp, 'Amniote_Database_Aug_2015.csv'), stringsAsFactors = FALSE)

foraging$binomial <- gsub(' ', '_', foraging$Scientific)
sismatchfor <- sister_names2 %in% foraging$binomial
sister_names2[!sismatchfor] # The original (wrong) names match here.

foraging <- foraging[,c(41, 10:20, 24:31, 35,36)]
foraging_sisters <- left_join(data.frame(binomial = sister_names2), foraging)

lifehist$binomial <- paste(lifehist$genus, lifehist$species, sep = '_')
sismatchlh <- sister_names2 %in% lifehist$binomial
sister_names2[!sismatchlh]

corr_names_lh <- c('Peucaea_sumichrasti', 'Melaenornis_microrhynchus', 'Chionodacryon_speculiferum', 'Euplectes_psammacromius', 'Basilinna_leucotis', 'Setophaga_pitiayumi', 'Melozone_albicollis', 'Crithagra_sulphurata', 'Peucaea_carpalis', 'Melaenornis_mariquensis', 'Basilinna_xantusii', 'Niltava_davidi', 'Setophaga_americana', 'Melozone_fusca', 'Crithagra_flaviventris')

sister_names_lh <- sister_names2
sister_names_lh[!sismatchlh] <- corr_names_lh


lifehist_sisters <- left_join(data.frame(binomial = sister_names_lh), lifehist[,c(37, 8:36)])

trait_sisters <- cbind(foraging_sisters, lifehist_sisters[,-1])

write.csv(trait_sisters, file = file.path(fp,'trait_sisters_23Jun.csv'), row.names = FALSE)

#### COVARIATE SET 11

# Compile the outputs that were run in parallel on the cluster.
# For each data point, get the richness and the richness of confamilial and congeneric individuals.

fp <- '~/verts/mapoutput'
lapply(as.list(file.path(fp, dir(fp, pattern = '23Jun'))), load, .GlobalEnv)
overlap_all <- do.call('c', lapply(paste0('overlap', 1:50), get))

bird <- read.csv('~/verts/birdcoords_23Jun.csv', stringsAsFactors = FALSE)
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

write.csv(bird, '~/verts/bird_withrichness_23Jun.csv', row.names = FALSE)


#### COVARIATE SET 12

# Number of subpopulations from which each set of specimens was collected. Need to calculate.
# Round lat and long to nearest 0.5'. Anything which shares those coordinates is a subpopulation.
nsubpops <- vnbird_good %>% 
  mutate(latround = plyr::round_any(decimallatitude, 1/120),
         lonround = plyr::round_any(decimallongitude, 1/120)) %>%
  group_by(binomial) %>%
  summarize(n = nrow(unique(cbind(latround, lonround))))

write.csv(nsubpops, file = file.path(fp, 'nsubpops_23Jun.csv'), row.names = FALSE)

#### COVARIATE SET 13

# Elevations were extracted on cluster.

library(dplyr)
fp <- ('C:/Users/Q/Dropbox/projects/verts')
bird <- read.csv(file.path(fp, 'birdcoords_23Jun.csv'), stringsAsFactors = FALSE)
bird_elev <- read.csv(file.path(fp, 'bird_elev_23Jun.csv'))
bird <- left_join(bird, bird_elev %>% rename(decimallatitude=lat, decimallongitude=lon))

# Get CV of elevation. Ignore repeats. Add minimum elevation so that everything is greater or equal to zero
bird_topovar <- bird %>%
  mutate(elev = elev + 77) %>%
  group_by(binomial) %>%
  summarize(elevation_cv = sd(unique(elev), na.rm = T)/mean(unique(elev), na.rm = T))

write.csv(bird_topovar, file = file.path(fp, 'topo_23Jun.csv'), row.names = FALSE)


#################################################################################

# Join all covariates to the long and wide form dataframes.

bird_clim <- read.csv(file.path(fp, 'bird_clim_summstats_23Jun.csv'), stringsAsFactors = FALSE)
bird_range <- read.csv(file.path(fp, 'botw_ranges_23Jun.csv'), stringsAsFactors = FALSE)
bird_trait <- read.csv(file.path(fp, 'trait_sisters_23Jun.csv'), stringsAsFactors = FALSE)
bird_richness <- read.csv(file.path(fp, 'bird_withrichness_23Jun.csv'), stringsAsFactors = FALSE)
bird_nsubpop <- read.csv(file.path(fp, 'nsubpops_23Jun.csv'), stringsAsFactors = FALSE)
bird_topovar <- read.csv(file.path(fp, 'topo_23Jun.csv'), stringsAsFactors = FALSE)

bird_richness_summ <- bird_richness %>%
  group_by(binomial) %>%
  summarize(total_richness = median(total_richness, na.rm=T),
            congener_richness = median(congener_richness, na.rm=T))

bird_allcov <- bird_clim %>% 
  left_join(rename(bird_range, binomial=taxon)) %>%
  left_join(bird_trait) %>%
  left_join(rename(bird_nsubpop, nsubpop=n)) %>%
  left_join(bird_topovar) %>% 
  left_join(bird_richness_summ)

sister1cov <- sister2cov <- bird_allcov
names(sister1cov) <- paste0(names(sister1cov), '1')
names(sister2cov) <- paste0(names(sister2cov), '2')
names(sister1cov)[1] <- 'sister1'
names(sister2cov)[1] <- 'sister2'

sister_alldat <- sister_troptemp %>%
  left_join(sister1cov) %>%
  left_join(sister2cov)

sister_alldat[sister_alldat == -999] <- NA

write.csv(sister_alldat, file = file.path(fp, 'sister_alldat_25Jun.csv'), row.names = FALSE)
