# Entire code pipeline for bird manuscript
# Revision for Bio lett

# 1. Redo identification of sisters (using consensus tree rather than each individual tree) - data source: birdtree
# 2. Redo identification of temperate vs tropical species (using range polygons rather than collection locations) - data source: botw
# 3. Combine vertnet data, identified sisters, and botw classification of temperate vs tropical into one data frame. Note that species IDs may differ in each of the 3 places.
# 4. Do the same analysis as previous manuscript for the main question.
# 5. Apply phylogenetic correction to the paired t-test (make a tree with basal nodes of each sister pair)
# 6. Do a better model selection procedure than before: forward selection, use AICc, and check VIFs of final model.

fp <- 'C:/Users/Q/Dropbox/projects/verts'
fprev <- file.path(fp, 'manuscript/revision')

library(dplyr)

# 1. Identification of sisters with consensus tree

# (Done on cluster with 100 trees sampled randomly)

library(ape)
birdtree <- read.tree(file.path(fp, 'ericson_cons.tre'))
ericsondist <- cophenetic.phylo(birdtree)

identify_sisters <- function(mat) {
  
  sisters <- character(nrow(mat))
  # Modified 08 June: add the phylogenetic distance of sister pairs as a column for later analysis.
  dist_to_sister <- numeric(nrow(mat))
  
  for (i in 1:nrow(mat)) {
    x <- mat[i, -i]
    name_i <- dimnames(mat)[[1]][i]
    dist_to_sister[i] <- min(x)
    sisternames <- names(x)[x == dist_to_sister[i]] # Names of potential sisters.
    # Check which of these sisters has taxon i as a sister.
    for (j in 1:length(sisternames)) {
      j_index <- which(dimnames(mat)[[1]] == sisternames[j])
      y <- mat[j_index, -j_index]
      reciprocaldist_to_sister <- min(y)
      reciprocalsisternames <- names(y)[y == reciprocaldist_to_sister]
      if (name_i %in% reciprocalsisternames) sisters[[i]] <- sisternames[j]
    }
    
  }
  
  sisters_df <- data.frame(sister1 = dimnames(mat)[[1]], sister2 = sisters, dist = dist_to_sister)
  sisters_df <- subset(sisters_df, sister2 != '')
  sisters_df <- transform(sisters_df, sister1 = as.character(sister1), sister2 = as.character(sister2))
  
  # Get rid of duplicates.
  is_dup <- rep(FALSE, nrow(sisters_df))
  
  for (i in 1:nrow(sisters_df)) {
    sisterrows <- which(sisters_df$sister2 == sisters_df$sister1[i] & sisters_df$sister1 == sisters_df$sister2[i])
    if (length(sisterrows) > 0)
      if (sisterrows[1] < i)
        is_dup[i] <- TRUE
  }
  
  return(sisters_df[!is_dup,])
  
}

sisters <- identify_sisters(ericsondist)

# Remove duplicates
sisters$duplicate <- FALSE
for (i in 2:nrow(sisters)) {
  previous <- c(sisters$sister1[1:(i-1)], sisters$sister2[1:(i-1)])
  if (sisters$sister1[i] %in% previous | sisters$sister2[i] %in% previous) sisters$duplicate[i] <- TRUE
}

sisters <- subset(sisters, !duplicate)

# 2. Identification of nontropical and tropical species with range polygons

botw_cents <- read.csv(file.path(fprev, 'polygon_centroids.csv'), stringsAsFactors = FALSE)

# Metadata location http://datazone.birdlife.org/species/spcdistPOS
# Use only resident and breeding polygons, native, and presence status 1-3
botw_filtered <- botw_cents %>%
  filter(ORIGIN %in% 1, SEASONAL %in% 1:2, PRESENCE %in% 1:3) %>%
  mutate(realm = 'tropical') %>%
  group_by(SCINAME, realm) %>%
  summarize(lat_centroid = weighted.mean(x = X2, w = Shape_Area),
            lon_centroid = weighted.mean(x = X1, w = Shape_Area)) %>%
  as.data.frame

botw_filtered$realm[abs(botw_filtered$lat_centroid) > 23.5] <- 'nontropical'

# addendum 2 Nov:
# How many sister temperate-tropical pairs are there if we don't care if they have adequate numbers in the vertnet database?
botw_filtered <- botw_filtered %>% 
  mutate(binomial = gsub('\\ ', '_', SCINAME)) 
teh_pairs <- sisters %>%
  left_join(botw_filtered %>% rename(sister1=binomial, realm1=realm) %>% dplyr::select(sister1, realm1)) %>%
  left_join(botw_filtered %>% rename(sister2=binomial, realm2=realm) %>% dplyr::select(sister2, realm2)) %>%
  filter(realm1 != realm2)

# 3. Quality control of raw vertnet data, then join with sister status (#1) and tropical status (#2)

vnbird <- read.csv(file.path(fp, 'vertnet_birds_reduced.csv'), stringsAsFactors = FALSE)	

correctnames <- read.csv(file.path(fp, 'correctednames.csv'), stringsAsFactors = FALSE)
for (i in 1:nrow(correctnames)) {
  fix_rows <- vnbird$genus == correctnames$old_genus[i] & vnbird$specificepithet == correctnames$old_species[i]
  vnbird$genus[fix_rows] <- correctnames$new_genus[i]
  vnbird$specificepithet[fix_rows] <- correctnames$new_species[i]
}

# Get rid of captive and fossil
vnbird_qc1 <- filter(vnbird, !wascaptive, !isfossil) # Leave in invasive for now.

# Retain only records with mass
vnbird_qc2 <- filter(vnbird_qc1, !is.na(massing))

# Attempt to get rid of juvenile records
lifet <- table(vnbird_qc2$lifestage)

remove_names <- grep('juv|sub', names(lifet), ignore.case = TRUE, value = TRUE)

vnbird_qc3 <- filter(vnbird_qc2, !lifestage %in% remove_names)

lifet3 <- table(vnbird_qc3$lifestage)

remove_names3 <- grep('im|nat|larv|young|nestl|fledg|unoss|pullus|hatch|chick|down', names(lifet3), ignore.case = TRUE, value = TRUE)
vnbird_qc4 <- filter(vnbird_qc3, !lifestage %in% remove_names3)

vnbird <- vnbird_qc4
rm(vnbird_qc1, vnbird_qc2, vnbird_qc3, vnbird_qc4)

vnbird <- filter(vnbird, scientificname != '')

# Join it all together, then output so that outliers can be manually flagged.

sister_names <- with(sisters, c(sister1, sister2))
vnbird <- vnbird %>% mutate(binomial = paste(genus, sapply(strsplit(specificepithet, ' '), '[', 1), sep = '_'))

# Get rid of anything less than 10 mass measurements
nmass <- vnbird %>% group_by(binomial) %>% summarize(nmassmeas = sum(!is.na(massing)))
vnbird <- vnbird %>% left_join(nmass) %>% filter(nmassmeas >= 10)

vnbird$issister <- vnbird$binomial %in% sister_names
missingnames <- vnbird %>% filter(!issister) %>% group_by(binomial) %>% summarize(n = n()) %>% arrange(-n)

vnsis <- vnbird %>% filter(issister)

# Join identified sisters with the temperate/tropical distinction.
# Filter out only the pairs that are temperate+tropical
# Then manually flag outliers in just that small-ish dataset.

botw_filtered <- mutate(botw_filtered, binomial = gsub('\\ ', '_', SCINAME))

vnsis_spnames <- unique(vnsis$binomial)

badnames <- vnsis_spnames[!vnsis_spnames %in% botw_filtered$binomial]
# Must correct about 300 species names to make these match. Argh.
# write.csv(data.frame(vertnet_name = badnames), file=file.path(fprev, 'correctednames16oct.csv'), row.names = FALSE)

#z <- function(x) grep(x, botw_filtered$binomial, value = TRUE)

# Load manually corrected names.
# goodnames <- read.csv(file.path(fprev, 'correctednames_botw_vertnet.csv'), stringsAsFactors = FALSE)

# Correct so that most of the names are already fixed
# badnames2 <- full_join(goodnames, data.frame(vertnet_name=badnames))
# write.csv(badnames2, file=file.path(fprev, 'correctednames16oct.csv'), row.names = FALSE)

goodnames <- read.csv(file.path(fprev, 'correctednames16oct_botw_vertnet.csv'), stringsAsFactors = FALSE)


# Change the BOTW names to Vertnet names.
goodnames <- mutate(goodnames,
                    vertnet_species = sapply(strsplit(vertnet_name, '_'), '[', 2))
goodnames$vertnet_species[!goodnames$botw_species==''] <- goodnames$botw_species[!goodnames$botw_species=='']
goodnames <- mutate(goodnames,
                    botw_name = paste(botw_genus, vertnet_species, sep = '_'))

for (i in 1:nrow(goodnames)) {
  idx <- which(botw_filtered$binomial == goodnames$botw_name[i])
  botw_filtered$binomial[idx] <- goodnames$vertnet_name[i] 
}

vnsis_spnames[!vnsis_spnames %in% botw_filtered$binomial]


# Correct the "lumped" species in vertnet.
vnsis$binomial_corr <- vnsis$binomial
to_correct <- c("Corapipo_leucorrhoa", "Picoides_arcticus", "Carduelis_hornemanni", 
"Calandrella_cheleensis", "Serinus_whytii", "Butorides_striata", 
"Aphelocoma_insularis", "Himantopus_himantopus")
lumped_spp <- c('Corapipo_altera', 'Picoides_tridactylus', 'Carduelis_flammea', 'Calandrella_rufescens', 'Serinus_striolatus', 'Butorides_virescens', 'Aphelocoma_californica', 'Himantopus_mexicanus')

for (i in 1:length(to_correct)) {
	vnsis$binomial_corr[vnsis$binomial_corr == to_correct[i]] <- lumped_spp[i]
}

vnsis_spnames <- unique(vnsis$binomial_corr)

badnames <- vnsis_spnames[!vnsis_spnames %in% botw_filtered$binomial] # No bad names!
badnames

# Join with temperate/tropical distinction.
vnsis_filter <- left_join(vnsis,
                         botw_filtered %>% select(binomial, realm) %>% rename(binomial_corr = binomial))

# Get only temperate-tropical pairs from vnsis.

# Add sister's name as a column to the data frame.
vnsis_filter <- left_join(vnsis_filter, sisters %>% rename(binomial_corr=sister1, sister_binomial=sister2))
vnsis_filter <- left_join(vnsis_filter, sisters %>% rename(binomial_corr=sister2, sister_binomial2=sister1, dist2=dist))
vnsis_filter$sister_binomial[!is.na(vnsis_filter$sister_binomial2)]	<- vnsis_filter$sister_binomial2[!is.na(vnsis_filter$sister_binomial2)]	
vnsis_filter$dist[!is.na(vnsis_filter$dist2)] <- vnsis_filter$dist2[!is.na(vnsis_filter$dist2)]

vnsis_filter <- vnsis_filter %>% select(-sister_binomial2, -dist2)
vnsis_filter$sister_realm <- vnsis_filter$realm[match(vnsis_filter$sister_binomial, vnsis_filter$binomial_corr)]
vnsis_filter_troptemp <- filter(vnsis_filter, (realm=='tropical' & sister_realm=='nontropical') | (realm=='nontropical' & sister_realm=='tropical'))


bird_summ <- vnsis_filter_troptemp %>%
	group_by(binomial_corr, sister_binomial, realm, sister_realm) %>%
	summarize(n = n(),
			  lat = median(abs(decimallatitude), na.rm=T),
			  cv_logmass = sd(log10(massing))/mean(log10(massing)))
			  
nontropical_species <- with(bird_summ, c(binomial_corr[realm=='nontropical'], sister_binomial[sister_realm=='nontropical']))
tropical_species <- with(bird_summ, c(sister_binomial[sister_realm=='tropical'], binomial_corr[realm=='tropical']))

vnfactor <- vnsis_filter %>% group_by(binomial_corr) %>%
  summarize(median_mass = median(massing, na.rm = TRUE))
vnsis_filter <- vnsis_filter %>% left_join(vnfactor) %>% mutate(outlierflag = massing > median_mass * 10 | massing < median_mass / 10)
table(vnsis_filter$outlierflag)


#### Final QC of coordinates and mass.

# Export data and inspect the 100 species' records to see if any locations are wrong.
# Put together tropical and nontropical
#vnsis_export <- vnsis_filter %>% 
#  filter(!outlierflag) %>%
#  filter(binomial_corr %in% c(tropical_species, nontropical_species))

# write.csv(vnsis_export, file.path(fprev, 'vertnet_sister_data.csv'), row.names = FALSE)

# Load corrected coordinates, replace the corrected entries in the proper columns, then save this for covariate extraction and further analysis.
vnsis_filter <- read.csv(file.path(fprev, 'vertnet_sister_data_corrected.csv'), stringsAsFactors = FALSE)

vnsis_filter$decimallatitude[!is.na(vnsis_filter$latitude_corr)] <- vnsis_filter$latitude_corr[!is.na(vnsis_filter$latitude_corr)]
vnsis_filter$decimallongitude[!is.na(vnsis_filter$longitude_corr)] <- vnsis_filter$longitude_corr[!is.na(vnsis_filter$longitude_corr)]
vnsis_filter$massing[!is.na(vnsis_filter$mass_corr)] <- vnsis_filter$mass_corr[!is.na(vnsis_filter$mass_corr)]
vnsis_filter$decimallatitude[vnsis_filter$decimallatitude == -9999] <- NA
vnsis_filter$decimallongitude[vnsis_filter$decimallongitude == -9999] <- NA
vnsis_filter <- filter(vnsis_filter, invalid_flag != 'y')

# Reexport corrected.
vnsis_filter %>% 
  select(binomial, decimallatitude, decimallongitude) %>%
  arrange(binomial, decimallatitude, decimallongitude) %>%
  write.csv(file.path(fprev, 'vertnet_coords.csv'), row.names = FALSE)


nontropical_summ <- vnsis_filter %>% 
	filter(!outlierflag) %>%
	filter(binomial_corr %in% nontropical_species) %>% 
	group_by(binomial_corr) %>%
	summarize(nontropical_n = n(),
			  nontropical_lat = median(abs(decimallatitude), na.rm=T),
			  nontropical_cv_logmass = sd(log10(massing))/mean(log10(massing))) %>%
	rename(nontropical_binomial = binomial_corr)
	
tropical_summ <- vnsis_filter %>% 
	filter(!outlierflag) %>%
	filter(binomial_corr %in% tropical_species) %>% 
	group_by(binomial_corr) %>%
	summarize(tropical_n = n(),
			  tropical_lat = median(abs(decimallatitude), na.rm=T),
			  tropical_cv_logmass = sd(log10(massing))/mean(log10(massing))) %>%
	rename(tropical_binomial = binomial_corr)

sister_pairs <- data.frame(nontropical_binomial=nontropical_species, tropical_binomial=tropical_species)	

sister_data <- nontropical_summ %>% left_join(sister_pairs) %>% left_join(tropical_summ)

with(sister_data, t.test(nontropical_cv_logmass, tropical_cv_logmass, paired = TRUE, alternative = 'greater')) 




# Make a tree with the basal nodes of all the sisters to be used as a phylogenetic correction.

# Tree with only the 100 species in it:
sistertree <- drop.tip(birdtree, tip = birdtree$tip.label[!birdtree$tip.label %in% c(tropical_species, nontropical_species)])

# Find common ancestor nodes of each sister pair.
library(phytools)
sister_mrcas <- sapply(1:length(tropical_species), function(i) findMRCA(sistertree, tips = c(tropical_species[i], nontropical_species[i])))

# Collapse each sister node to a "star" phylogeny

sistertree_fixed <- sistertree
for (i in 1:length(tropical_species)) {
  idx <- which(sistertree_fixed$tip.label %in% c(tropical_species[i], nontropical_species[i]))
  sistertree_fixed$edge.length[which(sistertree_fixed$edge[,2] %in% idx)] <- 0
}
sistertree_fixed <- drop.tip(sistertree_fixed, tip = nontropical_species) # Better.

nontr_dat <- sister_data$nontropical_cv_logmass
names(nontr_dat) <- sister_data$tropical_binomial
tr_dat <- sister_data$tropical_cv_logmass
names(tr_dat) <- sister_data$tropical_binomial

#### plot showing the difference as a color.

ds <- nontr_dat - tr_dat
d_to_color <- function(ds) {
  d_color <- colorRampPalette(RColorBrewer::brewer.pal(9, 'RdBu'))(99)
  break_seq <- seq(-max(abs(ds)), max(abs(ds)), length=100)
  d_cut <- cut(ds, breaks=break_seq)
  d_color[match(d_cut, levels(d_cut))]
}

# Color scale
plot(sistertree, type = 'fan')
nodelabels(text = ' ', node = sister_mrcas, bg = d_to_color(ds), width=1, height=1)

# Binary
plot(sistertree, type = 'fan')
nodelabels(text = ' ', node = sister_mrcas, bg = c('red','blue')[(ds > 0) + 1], width=1, height=1)

nontr_dat_order <- nontr_dat[sistertree_fixed$tip.label]
tr_dat_order <- tr_dat[sistertree_fixed$tip.label]
phylosig(sistertree_fixed, nontr_dat_order-tr_dat_order, method='lambda', test=T) # lambda = 0.712

# Make sure the phylogenetically independent contrasts return non-null values.
sistertree_di <- multi2di(sistertree_fixed)
pic(nontr_dat_order - tr_dat_order, sistertree_di) 

phyl.pairedttest(tree = sistertree_fixed, x1 = 100*cbind(nontr_dat, tr_dat))


#### Here insert bootstrap subsample analysis.

nsim <- 999
ttests <- list()

vnbird_good <- vnsis_filter %>% 
	filter(!outlierflag, binomial_corr %in% c(tropical_species, nontropical_species)) %>%
	select(binomial_corr, massing) %>%
	rename(binomial=binomial_corr)

set.seed(28462)

pb <- txtProgressBar(0, nsim, style = 3)


for (sim in 1:nsim) {
  
  sister_subsample <- sister_data
  
  
  for (i in 1:nrow(sister_data)) {
    if (sister_data$nontropical_n[i] < sister_data$tropical_n[i]) {
      tropical_subsample <- log10(sample(vnbird_good$massing[vnbird_good$binomial == sister_data$tropical_binomial[i]], size = sister_data$nontropical_n[i], replace = FALSE))
      sister_subsample$tropical_cv_logmass[i] <- sd(tropical_subsample)/mean(tropical_subsample)
    }
    if (sister_data$tropical_n[i] < sister_data$nontropical_n[i]) {
      nontropical_subsample <- log10(sample(vnbird_good$massing[vnbird_good$binomial == sister_data$nontropical_binomial[i]], size = sister_data$tropical_n[i], replace = FALSE))
      sister_subsample$nontropical_cv_logmass[i] <- sd(nontropical_subsample)/mean(nontropical_subsample)
    }
  }
  
  ttests[[sim]] <- with(sister_subsample, t.test(nontropical_cv_logmass, tropical_cv_logmass, paired = TRUE, alternative = 'greater'))
  setTxtProgressBar(pb, sim)
  
}

close(pb)

ps <- c(.5,.025,.975)
quantile(sapply(ttests, function(x) as.numeric(x$estimate)), prob=ps) # 0.007819834 0.004301053 0.012284541
table(sapply(ttests, function(x) as.numeric(x$estimate) > 0)) # true in all cases.
quantile(sapply(ttests, function(x) as.numeric(x$stat)), prob=ps) # 1.923627 1.078500 2.811414 
quantile(sapply(ttests, function(x) as.numeric(x$p.val)), prob=ps) # 0.029325732 0.003231461 0.142339269 
sum( sapply(ttests, function(x) as.numeric(x$p.val)) < 0.05)/999 #  0.7167167


#### Here insert covariate analysis.

# Load covariates and calculate difference for the different pairs.
all_covariates <- read.csv(file.path(fprev, 'all_covariates.csv'), stringsAsFactors = FALSE)
all_distances <- read.csv(file.path(fprev, 'all_distances.csv'), stringsAsFactors = FALSE)

tropical_covariates <- all_covariates[match(tropical_species, all_covariates$binomial), ]
nontropical_covariates <- all_covariates[match(nontropical_species, all_covariates$binomial), ]
names(tropical_covariates) <- paste('tropical', names(tropical_covariates), sep = '_')
names(nontropical_covariates) <- paste('nontropical', names(nontropical_covariates), sep = '_')

all_sister_data <- left_join(sister_data, all_distances) %>%
  left_join(tropical_covariates) %>%
  left_join(nontropical_covariates)

# Export for supplement.
all_sister_data %>% select(-lat_trop, -lon_trop, -lat_nontrop, -lon_nontrop) %>%
  write.csv(file = file.path(fprev, 'archive/bird_sister_species_data.csv'), row.names = FALSE)

# Export raw body mass data too
write.csv(vnbird_good, file = file.path(fprev, 'archive/raw_body_mass_data.csv'), row.names = FALSE)

# Calculate differences for the different pairs.


predtypes <- c('Invertebrate','Omnivore','VertFishScav')

multregdat <- all_sister_data %>% mutate(spatial_temp = nontropical_spatial_cv_temp - tropical_spatial_cv_temp,
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
                                       predator = nontropical_Diet.5Cat %in% predtypes | tropical_Diet.5Cat %in% predtypes,
                                       migrant = nontropical_migrant_status %in% c('partial','obligate') | tropical_migrant_status %in% c('partial','obligate'),
                                       dist_phylo = dist_phy,
                                       dcv = nontropical_cv_logmass - tropical_cv_logmass) %>%
  dplyr::select(dcv, dlat, dist_phylo, dist_slc, spatial_temp, interann_temp, spatial_precip, interann_precip, rangesize, total_richness, congener_richness, n_pop, elev_cv, seasonal_temp, seasonal_precip, area, meanlogbodymass, predator, migrant) 

# Added 08 Dec. Look at how many have multiple predator and migrant types.
# Contingency table
with(all_sister_data, table(nontropical_Diet.5Cat %in% predtypes, tropical_Diet.5Cat %in% predtypes)) # most are both. 9 both not, 52 both yes, 7 opp.
with(all_sister_data, table(nontropical_migrant_status %in% c('partial','obligate'), tropical_migrant_status %in% c('partial','obligate'))) # Here we look at whether the NONtropical is a migrant, since no species pairs have the tropical being a migrant and the nontropical being one.

# Create new predator variable with 3 categories: both herbivore, both carn/omnivore, and mixed.
predator_mixed <- with(all_sister_data, xor(nontropical_Diet.5Cat %in% predtypes, tropical_Diet.5Cat %in% predtypes))
predator_both <- with(all_sister_data, nontropical_Diet.5Cat %in% predtypes & tropical_Diet.5Cat %in% predtypes)

multregdat <- mutate(multregdat, predator_mixed=predator_mixed, predator_both=predator_both) %>% dplyr::select(-predator) # Only run this line if you want to use the alternate (3-way) predator classification

multreg_complete <- filter(multregdat, complete.cases(multregdat)) %>% as.data.frame

multreg_pred <- multreg_complete[,-1]
cor(multreg_pred)
usdm::vif(multreg_pred) # Interannual temperature variation and seasonal temperature variation are closely related. However we're interested in which of these might be driving things so leave them both in for now.

lm_full <- lm(dcv ~ ., data = multreg_complete, na.action = 'na.pass')

# Old version of model selection. Used AICc all along.

library(MuMIn)
lm_dredge <- dredge(lm_full, m.lim = c(0,5))
lm_best <- subset(lm_dredge, delta < 2)

plot(lm_best)

lm_averaged <-model.avg(lm_best, revised.var = TRUE) 

lm_best1 <- lm(dcv ~ interann_temp + migrant + predator + seasonal_precip + spatial_temp, data = multreg_complete)


# Stepwise. Forward and backward. Does not use AICc.
# Comes up with the same result as the "naive" dredge approach.
min_model <- lm(dcv ~ 1, data = multreg_complete)
full_model <- lm(dcv ~ ., data = multreg_complete)
full_model_formula <- formula(full_model)
# These yield the same result:
step(min_model, direction = 'forward', scope = full_model_formula)
stepAIC(min_model, direction = 'forward', scope = full_model_formula)

# Try with AICc
source('~/GitHub/birdtraits/revision/stepaicc.r')
aicc_step_result <- stepAICc(min_model, direction = 'forward', scope = full_model_formula)

# Standardized coefficients
best_model_data <- multreg_complete %>%
  transmute(dcv = dcv,
            migrant = migrant,
            predator = predator,
            spatial_temp = spatial_temp/sd(spatial_temp),
            interann_temp = interann_temp/sd(interann_temp),
            seasonal_precip = seasonal_precip/sd(seasonal_precip))
best_model_std <- lm(dcv ~ ., data = best_model_data)

# Do this with two predator status variables
best_model_data <- multreg_complete %>%
  transmute(dcv = dcv,
            migrant = migrant,
            predator_both = predator_both,
            spatial_temp = spatial_temp/sd(spatial_temp),
            rangesize = rangesize/sd(rangesize))
best_model_std <- lm(dcv ~ ., data = best_model_data)


# Correlations
cor(multreg_complete[,c('seasonal_precip','interann_temp','spatial_temp')])

# Added 08 Dec. Get R^2 and VIF analysis for best model to add to MS.
summary(best_model_std)
usdm::vif(best_model_data[,-1])

#####
# Revised figures.
library(cowplot)

### MAP
# Recreate data.

mapdat_sister <- all_sister_data %>% mutate(d = tropical_cv_logmass - nontropical_cv_logmass)

heatramp <- colorRampPalette(RColorBrewer::brewer.pal(name='YlOrRd', n=9),bias=2,space="rgb")(50)
fillScale <- scale_fill_gradientn(colours = heatramp, name = 'Body mass variability', breaks=c(.05,.1,.15), labels=c('.05','.10','.15'))
linetypeScale <- scale_linetype_manual(name = '', values = c('dotted','solid'), labels = c('Tropical\nmore variable','Nontropical\nmore variable'))

worldMap <- borders('world', fill='gray75', color='black')
p_map <- ggplot() + worldMap + 
  geom_segment(aes(x=lon_trop, y=lat_trop, xend=lon_nontrop, yend=lat_nontrop, linetype = d>0), data=mapdat_sister) +
  geom_point(aes(x=lon_trop,y=lat_trop, fill=tropical_cv_logmass), data=mapdat_sister, pch=21) +
  geom_point(aes(x=lon_nontrop,y=lat_nontrop, fill=nontropical_cv_logmass), data=mapdat_sister, pch=21) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(breaks = c(-23.5, 23.5), labels = c('23.5° S', '23.5° N'), expand=c(0,0)) +
  fillScale + linetypeScale +
  theme(panel.background = element_rect(fill='skyblue'),
        panel.grid.major.y = element_line(color='black', size=1.1),
        panel.grid.major.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        legend.position=c(0,0), 
        legend.direction = 'horizontal',
        legend.justification = c(0,0),
        legend.text = element_text(size=10)) + 
  coord_map('albers',0,0, ylim = c(-60,60))


####################################################################################################
# Plot interaction plot.

ylabel <- expression(paste('CV of log'[10], ' body mass'))		

plotdat <- data.frame(realm = rep(c('tropical','nontropical'), each=nrow(sister_data)),
                      taxon = c(sister_data$tropical_binomial, sister_data$nontropical_binomial),
                      pairid = 1:nrow(sister_data),
                      cv_logmass = c(sister_data$tropical_cv_logmass, sister_data$nontropical_cv_logmass))

p_int <- ggplot(plotdat, aes(x = realm, y = cv_logmass)) +
  geom_line(aes(group = pairid), color = 'gray75') + 
  geom_point(aes(group = pairid)) +
  stat_summary(aes(group = 1), geom = 'line', color = 'red', fun.y = 'mean', size = 0.8) +
  stat_summary(geom = 'errorbar', color = 'red', fun.data = 'mean_se', width = 0.1) +
  scale_x_discrete(expand=c(0.1,0.1), name='Zone') +
  labs(y = ylabel) +
  panel_border(colour='black')

phist <- ggplot(sister_data, aes(x = nontropical_cv_logmass - tropical_cv_logmass)) +
  geom_histogram(bins = 15) +
  geom_vline(xintercept = 0.0092, color = 'red', lwd = 2) +
  geom_vline(xintercept = 0.0029, color = 'red', lwd = 0.5) +
  geom_vline(xintercept = 0, color = 'blue', lty = 3) +
  scale_y_continuous(limits = c(0,25), expand = c(0,0)) +
  labs(x = 'nontropical CV - tropical CV', y = 'count of species pairs') +
  panel_border(colour='black')

####################################################################################################
# Plot scatterplots of covariates that go into the main text.
# Must use the version that controls for other predictors.
# Added variable plots.

library(car)

avPlots(best_model_std)
avdat <- avPlots(best_model_std)
avdat <- lapply(avdat, function(x) {x <- as.data.frame(x); names(x) <- c('x','y'); x})

avreg1 <- lm(y ~ x, data = avdat$interann_temp)
avcoef1 <- avreg1$coefficients
avreg2 <- lm(y ~ x, data = avdat$spatial_temp)
avcoef2 <- avreg2$coefficients

hl <- geom_hline(yintercept = 0, color = 'gray80', lwd = 0.4)
vl <- geom_vline(xintercept = 0, color = 'gray80', lwd = 0.4)

y_scale <- scale_y_continuous(limits = c(-0.07, 0.11))

pcov1 <- ggplot(avdat$interann_temp, aes(x = x, y = y)) +
  hl + 
  geom_point() + stat_smooth(method='lm', se=F) +
  panel_border(colour='black') + 
  y_scale +
  labs(x = expression(paste('interannual ',Delta * CV[temperature])), y = expression(Delta * CV[bodymass])) +
  geom_text(data=data.frame(x=c(Inf, -Inf), y = c(Inf, -Inf), lab = c('Nontropical \nmore variable ', ' Tropical\n more variable')),
            aes(label=lab), hjust = c(1,0), vjust = c(1.1,-.5), size = 3)

pcov2 <- ggplot(avdat$spatial_temp, aes(x = x, y = y)) +
  hl + 
  geom_point() + stat_smooth(method='lm', se=F) +
  panel_border(colour='black') + 
  y_scale +
  labs(x = expression(paste('spatial ',Delta * CV[temperature])), y = expression(Delta * CV[bodymass])) +
  geom_text(data=data.frame(x=c(Inf, -Inf), y = c(Inf, -Inf), lab = c('Nontropical \nmore variable ', ' Tropical\n more variable')),
            aes(label=lab), hjust = c(1,0), vjust = c(1.1,-.5), size = 3)

##### Combine

smalllab <- theme(axis.title.x=element_text(size=10),
                  axis.title.y=element_text(size=10))

bottom_row <- plot_grid(p_int, pcov1, pcov2, labels = c('b', 'c', 'd'), align = 'h', rel_widths = c(1, 1.4, 1.4), ncol=3)
full_plot <- plot_grid(p_map + panel_border(colour='black') + theme(legend.position='bottom'), bottom_row, labels = c('a', ''), ncol = 1, rel_heights = c(1, 0.9))
ggsave(file.path(fprev, 'fig1_finaledit.png'), full_plot, height=6, width=8, dpi=400)


#### Supplemental figure.
####################################################################################################
# Plot scatterplots of covariates that go into the supplement.
# Edit 08 Dec. Generate added variable plots for these, as well as adding a line to the plots.

xaxisvars <- c('dlat', 'dist_phylo', 'dist_slc', 'spatial_temp','interann_temp','seasonal_temp','spatial_precip','interann_precip','seasonal_precip','rangesize','total_richness','congener_richness','n_pop','elev_cv','area', 'predator', 'migrant', 'meanlogbodymass')
xaxislabels <- c('Latitudinal distance', 'Phylogenetic distance', 'Geographic distance', 'Spatial CV of MAT', 'Interannual CV of MAT', 'Seasonal CV of MAT', 'Spatial CV of MAP', 'Interannual CV of MAP', 'Seasonal CV of MAP', 'Range size (log10 km2)', 'Total richness', 'Congener richness', 'Populations collected', 'CV of elevation', 'Collection area (log10 km2)', 'Predator status', 'Migrant status', 'Mean of log body mass')

# Create insane full model and do added variable regression to get data for plots
library(car)
lm_full <- lm(dcv ~ ., data = multreg_complete, na.action = 'na.pass')
avdatfull <- avPlots(lm_full)
avdatfull <- lapply(avdatfull, function(x) {x <- as.data.frame(x); names(x) <- c('x','y'); x})
xaxislabels <- c('Latitudinal distance', 'Phylogenetic distance', 'Geographic distance', 'Spatial CV of MAT', 'Interannual CV of MAT', 'Spatial CV of MAP', 'Interannual CV of MAP', 'Range size (log10 km2)', 'Total richness', 'Congener richness', 'Populations collected', 'CV of elevation', 'Seasonal CV of MAT', 'Seasonal CV of MAP', 'Collection area (log10 km2)', 'Mean of log body mass', 'Predator status', 'Migrant status')


scatterplots <- list()

for (i in 1:length(xaxislabels)) {
  p_i <- ggplot(avdatfull[[i]], aes(x, y)) +
    geom_point() +
    stat_smooth(method = 'lm', se = FALSE) +
    panel_border(colour = 'black') +
    labs(y = expression(Delta*CV[bodymass]), x = xaxislabels[i])
  scatterplots[[i]] <- p_i
}

p_all <- plot_grid(plotlist = scatterplots, ncol = 3, labels = letters[1:length(xaxislabels)])
ggsave(file.path(fprev, 'supplementalfig1.png'), p_all, height = 15*.85, width = 11*.9, dpi = 400)
       
#### Model selection plot

dimlabels <- c('(Intercept)', 'log10 collection area', 'Congener richness', 'Phylogenetic distance', 'Geographic distance', 'Latitudinal distance', 'CV of elevation', 'Interannual CV of MAP', 'Interannual CV of MAT', 'Mean of log body mass', 'Migrant status', 'Populations collected', 'Predator status', 'log10 range size', 'Seasonal CV of MAP', 'Seasonal CV of MAT', 'Spatial CV of MAP', 'Spatial CV of MAT', 'Total richness')
  
png(file.path(fprev, 'supplementalfig2.png'), height = 9, width = 8, res = 400, units = 'in')
 par(mar = c(1, 4.5, 12, 4))
 plot(lm_best, labels=dimlabels)
dev.off()


# Edit: add AICcs to supplemental plot. Don't plot the weird numbers.
lm_best <- subset(lm_dredge, delta < 2)
deltaaiccs <- lm_best$AICc - (min(lm_best$AICc))

lm_best_fixlabels <- lm_best
dimnames(lm_best_fixlabels)[[1]] <- round(deltaaiccs, 2)

plot(lm_best_fixlabels, labels=dimlabels)

png(file.path(fprev, 'supplementalfig2_withaiccs.png'), height = 9, width = 8, res = 400, units = 'in')
par(mar = c(1, 4.5, 12, 4.1))
plot(lm_best_fixlabels, labels=dimlabels)
dev.off()