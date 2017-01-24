# Find sister taxa in vertnet.

fp <- '/mnt/research/plz-lab/NEON/external_data/raw_external_data/vertnet'
vnmam <- read.csv(file.path(fp, 'vertnet_mammals_reduced.csv'), stringsAsFactors = FALSE)	

# Figure out how many sets of congeners there are.

library(dplyr)

congeners <- vnmam %>% 
	filter(genus != '', specificepithet != '') %>% 
	group_by(genus) %>% 
	summarize(n = length(unique(specificepithet)))

sister_taxa <- congeners$genus[congeners$n > 1]

# There are several hundred genera with multiple species.

# Next, find the median latitude of each species and the amount of itv in their body mass.

vert_itv <- vnmam %>% 
	filter(genus %in% sister_taxa, specificepithet != '') %>%
	mutate(logmass = log10(massing)) %>%
	group_by(genus, specificepithet) %>% 
	summarize(lat = median(decimallatitude, na.rm=T),
		      cv_length = sd(lengthinmm, na.rm=T)/mean(lengthinmm, na.rm=T),
			  cv_logmass = sd(logmass, na.rm=T)/mean(logmass, na.rm=T),
			  n_length = sum(!is.na(lengthinmm)),
			  n_logmass = sum(!is.na(logmass)))
			  
# Subset out the ones with adequate sample size (5 individuals)
vert_subset <- filter(vert_itv, n_logmass >= 5)

# Get rid of anything that has only one species in genus after that.
tgenus <- table(vert_subset$genus)
vert_subset <- filter(vert_subset, genus %in% names(tgenus)[tgenus > 1])
# Result: 226 genera.

# Plot.
ggplot(vert_subset, aes(x = abs(lat), y = cv_logmass, group = genus, color = genus)) +
	geom_point() + stat_smooth(method = lm, se = FALSE) +
	theme_minimal() + theme(legend.position = 'none')

##################################################
	
# Basic mixed model. Random intercept for each genus.
# Must use link function for cv_logmass.
# library(lme4)
# sister_lmm <- glmer(cv_logmass ~ abslat + (1|genus), data = mutate(vert_subset, abslat = abs(lat)), family = Gamma(link='log'))
# sister_lmm2 <- glmer(cv_logmass ~ abslat + (1|genus), data = mutate(vert_subset, abslat = abs(lat)), family = gaussian(link='log'))

# summary(sister_lmm)
# confint(sister_lmm, method='boot', nsim=999, parm='abslat')

##################################################

# Run as paired t-test.
# Select the most tropical and the most polar from each genus.

cv_pairs <- vert_subset %>%
	summarize(lat_tropical = min(abs(lat)),
			  lat_polar = max(abs(lat)),
			  cv_logmass_tropical = cv_logmass[which(lat == min(abs(lat)))[1]],
			  cv_logmass_polar = cv_logmass[which(lat == max(abs(lat)))[1]])

# Subset only those genera where the pair is separated by at least 10 degrees of latitude.
cv_pairs_10degrees <- filter(cv_pairs, lat_polar - lat_tropical > 10, !is.na(cv_logmass_tropical), !is.na(cv_logmass_polar))
# This retains 48 pairs of taxa.

# Run paired t-test. One-sided hypothesis because theory predicts that polar will have a greater CV than tropical.
# So we set alternative = 'less' meaning we expect x < y, or tropical < polar.
with(cv_pairs_10degrees, t.test(cv_logmass_tropical, cv_logmass_polar, paired = TRUE, alternative = 'less'))

# Try to subset the pairs that cross the tropical boundary
cv_pairs_tropic <- filter(cv_pairs, lat_polar > 23.5, lat_tropical < 23.5, !is.na(cv_logmass_tropical), !is.na(cv_logmass_polar))
# This retains 15 pairs of taxa.

with(cv_pairs_tropic, t.test(cv_logmass_tropical, cv_logmass_polar, paired = TRUE, alternative = 'less'))

write.csv(cv_pairs_tropic, file = '~/cv_pairs_mammals.csv', row.names = FALSE)


###############################################################
###############################################################
###############################################################

# Added 19 Jan: Test the hypothesis in all taxa.

fp <- '/mnt/research/plz-lab/NEON/external_data/raw_external_data/vertnet'
vnmam <- read.csv(file.path(fp, 'vertnet_mammals_reduced.csv'), stringsAsFactors = FALSE)	
vnbird <- read.csv(file.path(fp, 'vertnet_birds_reduced.csv'), stringsAsFactors = FALSE)	
vnamph <- read.csv(file.path(fp, 'vertnet_amphibians_reduced.csv'), stringsAsFactors = FALSE)	
vnrept <- read.csv(file.path(fp, 'vertnet_reptiles_reduced.csv'), stringsAsFactors = FALSE)	
vnfish <- read.csv(file.path(fp, 'vertnet_fishes_reduced.csv'), stringsAsFactors = FALSE)	

# Correct the mistake in the red-winged blackbird where some of the species names are spelled wrong.
# Also, Yuhina zantholeuca should not be in the genus Yuhina. It isn't even in the same family. Update that name.
# Also, Helmitheros vermivorum is a monotypic genus but some of the specimens are labeled vermivorus (bad Graeco-Latin grammar)
# vnbird$specificepithet[vnbird$genus == 'Agelaius' & grepl('phoe', vnbird$specificepithet)] <- 'phoeniceus'
# vnbird$genus[vnbird$genus == 'Yuhina' & vnbird$specificepithet == 'zantholeuca'] <- 'Erpornis'
# vnbird$specificepithet[vnbird$genus == 'Helmitheros' & vnbird$specificepithet == 'vermivorus'] <- 'vermivorum'

# 20 Jan 2017: updated name correction
correctnames <- read.csv('~/verts/correctednames.csv', stringsAsFactors = FALSE)
for (i in 1:nrow(correctnames)) {
	fix_rows <- vnbird$genus == correctnames$old_genus[i] & vnbird$specificepithet == correctnames$old_species[i]
	vnbird$genus[fix_rows] <- correctnames$new_genus[i]
	vnbird$specificepithet[fix_rows] <- correctnames$new_species[i]
}

library(dplyr)

get_sister <- function(vn) {

	congeners <- vn %>% 
		filter(genus != '', specificepithet != '') %>% 
		group_by(genus) %>% 
		summarize(n = length(unique(specificepithet)))

	sister_taxa <- congeners$genus[congeners$n > 1]

	# Next, find the median latitude of each species and the amount of itv in their body mass.

	vert_itv <- vn %>% 
		filter(genus %in% sister_taxa, specificepithet != '') %>%
		mutate(logmass = log10(massing)) %>%
		group_by(genus, specificepithet) %>% 
		summarize(lat = median(decimallatitude, na.rm=T),
				  cv_length = sd(lengthinmm, na.rm=T)/mean(lengthinmm, na.rm=T),
				  cv_logmass = sd(logmass, na.rm=T)/mean(logmass, na.rm=T),
				  n_length = sum(!is.na(lengthinmm)),
				  n_logmass = sum(!is.na(logmass)))
				  
	# Subset out the ones with adequate sample size (5 individuals)
	vert_subset <- filter(vert_itv, n_logmass >= 5)

	# Get rid of anything that has only one species in genus after that.
	tgenus <- table(vert_subset$genus)
	vert_subset <- filter(vert_subset, genus %in% names(tgenus)[tgenus > 1])
	
	return(vert_subset)

}

mam_subset <- get_sister(vnmam)
bird_subset <- get_sister(vnbird)
amph_subset <- get_sister(vnamph)
rept_subset <- get_sister(vnrept)
fish_subset <- get_sister(vnfish)

get_cv_pairs <- function(vns) {
	vns %>%
		summarize(species_tropical = specificepithet[which(lat == min(abs(lat)))[1]],
				  species_nontropical = specificepithet[which(lat == max(abs(lat)))[1]],
				  lat_tropical = min(abs(lat)),
				  lat_nontropical = max(abs(lat)),
				  cv_logmass_tropical = cv_logmass[which(lat == min(abs(lat)))[1]],
				  cv_logmass_nontropical = cv_logmass[which(lat == max(abs(lat)))[1]],
				  n_logmass_tropical = n_logmass[which(lat == min(abs(lat)))[1]],
				  n_logmass_nontropical = n_logmass[which(lat == max(abs(lat)))[1]]) %>%
		filter(lat_nontropical > 23.5, lat_tropical < 23.5, !is.na(cv_logmass_tropical), !is.na(cv_logmass_nontropical))
}

mam_pairs <- get_cv_pairs(mam_subset) # 15 valid pairs
bird_pairs <- get_cv_pairs(bird_subset) # 97 valid pairs (one was removed due to typo)
amph_pairs <- get_cv_pairs(amph_subset) # No valid pairs
rept_pairs <- get_cv_pairs(rept_subset) # Only 2 valid pairs
fish_pairs <- get_cv_pairs(fish_subset) # No valid pairs

with(mam_pairs, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less'))
with(bird_pairs, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less')) # Significant!!!!
# with(amph_pairs, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less'))
# with(rept_pairs, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less'))
# with(fish_pairs, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less'))

# Save all the pairs.
vert_pairs <- rbind(cbind(taxon = 'mammal', mam_pairs),
					cbind(taxon = 'bird', bird_pairs),
					cbind(taxon = 'reptile', rept_pairs))
					
with(vert_pairs, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less')) # This is also significant but driven by birds (98 of the 115 pairs)

write.csv(vert_pairs, file = '~/verts/cv_pairs_vertebrates.csv', row.names = FALSE)

# Make some plots to show the pairs.

library(ggplot2)

# Convert paired data frame to longer form data frame.
vert_pair_long <-
with(vert_pairs, data.frame(taxon = taxon, realm = c(rep('tropical', nrow(vert_pairs)), rep('nontropical', nrow(vert_pairs))), genus = genus,
							species = c(species_tropical, species_nontropical), lat = c(lat_tropical, lat_nontropical),
							cv_logmass = c(cv_logmass_tropical, cv_logmass_nontropical), n = c(n_logmass_tropical, n_logmass_nontropical)))

ylabel <- expression(paste('CV of log'[10], ' body mass'))							
							
pbird <- ggplot(filter(vert_pair_long, taxon == 'bird'), aes(x = realm, y = cv_logmass)) +
	geom_line(aes(group = genus), color = 'gray75') + geom_point(aes(group = genus)) +
	stat_summary(aes(group = 1), geom = 'line', color = 'red', fun.y = 'mean') +
	stat_summary(geom = 'pointrange', color = 'red', fun.data = 'mean_se') +
	theme_bw() + theme(panel.grid = element_blank()) +
	ggtitle('Body mass variation in bird sister taxa') + labs(y = ylabel)
	
pmam <- ggplot(filter(vert_pair_long, taxon == 'mammal'), aes(x = realm, y = cv_logmass)) +
	geom_line(aes(group = genus), color = 'gray75') + geom_point(aes(group = genus)) + 
	stat_summary(aes(group = 1), geom = 'line', color = 'red', fun.y = 'mean') +
	stat_summary(geom = 'pointrange', color = 'red', fun.data = 'mean_se') +
	theme_bw() + theme(panel.grid = element_blank()) +
	ggtitle('Body mass variation in mammal sister taxa') + labs(y = ylabel)
	
ggsave('~/verts/sistertaxa_birds.png', pbird, height=6, width=6, dpi=300)
ggsave('~/verts/sistertaxa_mammals.png', pmam, height=6, width=6, dpi=300)

# Histogram showing the differences are normally distributed

phist <- ggplot(filter(vert_pairs, taxon == 'bird'), aes(x = cv_logmass_tropical - cv_logmass_nontropical)) +
	geom_histogram(bins = 10) +
	geom_vline(xintercept = -0.049, color = 'red', lwd = 2) +
	geom_vline(xintercept = -0.031, color = 'red', lwd = 0.5) +
	geom_vline(xintercept = 0, color = 'blue', lty = 3) +
	theme_bw() + theme(panel.grid = element_blank()) +
	scale_y_continuous(limits = c(0,35), expand = c(0,0)) +
	labs(x = 'tropical CV - nontropical CV')
	
ggsave('~/verts/hist_meandiffs.png', phist, height=6, width=6, dpi=300)


############################################################

# Diagnostics: does sample size have any relation to the CV

ggplot(filter(vert_pair_long, taxon=='bird'), aes(x=n, y=cv_logmass, color=realm)) + geom_point() + theme_minimal() + scale_x_log10()

cvn_lm <- lm(cv_logmass ~ n, data = filter(vert_pair_long, taxon=='bird'))
cvlogn_lm <- lm(cv_logmass ~ I(log(n)), data = filter(vert_pair_long, taxon=='bird'))

# There is a small tendency for the CV to go up with log(n). R squared is 0.05. Is this a problem? Dunno.

# Correction for birds: subsample the more numerous member of each species pair until it equals the lesser one. Then CV is compared for that subsample.

get_sister_allrecords <- function(vn) {

	congeners <- vn %>% 
		filter(genus != '', specificepithet != '') %>% 
		group_by(genus) %>% 
		summarize(n = length(unique(specificepithet)))

	sister_taxa <- congeners$genus[congeners$n > 1]

	# Next, find the median latitude of each species and the amount of itv in their body mass.

	vert_itv <- vn %>% 
		filter(genus %in% sister_taxa, specificepithet != '') %>%
		mutate(logmass = log10(massing)) %>%
		group_by(genus, specificepithet) %>% 
		do(data.frame(lat = .$decimallatitude,
					  lon = .$decimallongitude,
					  length = .$lengthinmm,
					  logmass = .$logmass,
					  n_length = sum(!is.na(.$lengthinmm)),
					  n_logmass = sum(!is.na(.$logmass))))
	
	# Subset out the ones with adequate sample size (5 individuals)
	vert_subset <- filter(vert_itv, n_logmass >= 5)

	# Get rid of anything that has only one species in genus after that.
	tgenus <- vert_subset %>% ungroup %>% group_by(genus) %>% summarize(nspp = length(unique(specificepithet))) %>% filter(nspp > 1)
	vert_subset <- filter(vert_subset, genus %in% tgenus$genus)
	
	return(vert_subset)

}

get_cv_pairs_subsample <- function(dat) {
	lats <- dat %>% group_by(specificepithet) %>% summarize(lat = median(lat, na.rm = TRUE))
	lat_trop <- min(abs(lats$lat), na.rm = TRUE)
	lat_nontrop <- max(abs(lats$lat), na.rm = TRUE)
	if (lat_trop > 23.5 | lat_nontrop < 23.5) return(data.frame(species_tropical = NA, 
																species_nontropical = NA,
																lat_tropical = NA,
																lat_nontropical = NA,
																cv_logmass_tropical = NA,
																cv_logmass_nontropical = NA,
																n_logmass_tropical = NA,
																n_logmass_nontropical = NA,
																n_use = NA))
	sp_trop <- lats$specificepithet[which(lats$lat == min(abs(lats$lat)))[1]]
	sp_nontrop <- lats$specificepithet[which(lats$lat == max(abs(lats$lat)))[1]]
	masses_trop <- na.omit(dat$logmass[dat$specificepithet == sp_trop])
	masses_nontrop <- na.omit(dat$logmass[dat$specificepithet == sp_nontrop])
	n_use <- min(length(masses_trop), length(masses_nontrop))
	masses_trop_use <- sample(masses_trop, size = n_use, replace = FALSE)
	masses_nontrop_use <- sample(masses_nontrop, size = n_use, replace = FALSE)
	cv_trop <- sd(masses_trop_use)/mean(masses_trop_use)
	cv_nontrop <- sd(masses_nontrop_use)/mean(masses_nontrop_use)
	return(data.frame(species_tropical = sp_trop, 
					  species_nontropical = sp_nontrop,
					  lat_tropical = lat_trop,
					  lat_nontropical = lat_nontrop,
					  cv_logmass_tropical = cv_trop,
					  cv_logmass_nontropical = cv_nontrop,
					  n_logmass_tropical = length(masses_trop),
					  n_logmass_nontropical = length(masses_nontrop),
					  n_use = n_use))
}

bird_subset_all <- get_sister_allrecords(vnbird)

# Still significant for one possible subsample.
bird_pairs_subsample <- bird_subset_all %>% ungroup %>% group_by(genus) %>% do(get_cv_pairs_subsample(.)) %>% filter(!is.na(cv_logmass_tropical))
with(bird_pairs_subsample, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less')) 

# Save coordinates and measurements for the right genera
bird_allrecords <- bird_subset_all %>% filter(genus %in% bird_pairs_subsample$genus)
write.csv(bird_allrecords, '~/verts/bird_allrecords.csv', row.names = FALSE)

# Run many iterations to see what the distribution of mean differences is.
nsim <- 999
tstats <- numeric(nsim)
pvals <- numeric(nsim)
meandiffs <- numeric(nsim)

pb <- txtProgressBar(0, nsim, style = 3)
options(dplyr.show_progress = FALSE)

for (i in 1:nsim) {
	bird_pairs_subsample_i <- bird_subset_all %>% ungroup %>% group_by(genus) %>% do(get_cv_pairs_subsample(.)) %>% filter(!is.na(cv_logmass_tropical))
	ttest_i <- with(bird_pairs_subsample_i, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less')) 
	tstats[i] <- ttest_i$stat
	pvals[i] <- ttest_i$p.value
	meandiffs[i] <- ttest_i$estimate
	setTxtProgressBar(pb, i)
}

close(pb)

quantile(tstats, probs = c(0.5, 0.025, 0.975))
quantile(pvals, probs = c(0.5, 0.025, 0.975))
quantile(meandiffs, probs = c(0.5, 0.025, 0.975))

write.csv(data.frame(run=1:nsim, tstat=tstats, pval=pvals, meandiff=meandiffs), file = '~/verts/subsampletstats.csv', row.names=FALSE)


#####################################################################

# Extract climate layers from the measurement coordinates (there are >90000 so this may be a laborious undertaking)
bird_allrecords <- read.csv('~/verts/bird_allrecords.csv', stringsAsFactors = FALSE)

# About 20000 of the records do not have a coordinate.

# Load climate data.


#####################################################################

# Get phylogenetic information for all the birds in the target genera. This includes ones that didn't make the 97 list.

fp <- '/mnt/research/plz-lab/NEON/external_data/raw_external_data/vertnet'
vnbird <- read.csv(file.path(fp, 'vertnet_birds_reduced.csv'), stringsAsFactors = FALSE)	
#vnbird$specificepithet[vnbird$genus == 'Agelaius' & grepl('phoe', vnbird$specificepithet)] <- 'phoeniceus'
#vnbird$genus[vnbird$genus == 'Yuhina' & vnbird$specificepithet == 'zantholeuca'] <- 'Erpornis'
#vnbird$specificepithet[vnbird$genus == 'Helmitheros' & vnbird$specificepithet == 'vermivorus'] <- 'vermivorum'

correctnames <- read.csv('~/verts/correctednames.csv', stringsAsFactors = FALSE)
for (i in 1:nrow(correctnames)) {
	fix_rows <- vnbird$genus == correctnames$old_genus[i] & vnbird$specificepithet == correctnames$old_species[i]
	vnbird$genus[fix_rows] <- correctnames$new_genus[i]
	vnbird$specificepithet[fix_rows] <- correctnames$new_species[i]
}


library(dplyr)

bird_subset_all <- get_sister_allrecords(vnbird)
bird_pairs_subsample <- bird_subset_all %>% ungroup %>% group_by(genus) %>% do(get_cv_pairs_subsample(.)) %>% filter(!is.na(cv_logmass_tropical))
bird_allrecords <- bird_subset_all %>% filter(genus %in% bird_pairs_subsample$genus) # Down to 774 species.

# Retain only the genera that have at least one tropical and one nontropical member.

#generapairs <- bird_subset_all %>% summarize(lat = median(lat, na.rm=T)) %>% ungroup %>% group_by(genus) %>% do(data.frame(oneeach = any(abs(.$lat) < 23.5) & any(abs(.$lat) > 23.5)))
#pairedgenera <- filter(generapairs, oneeach)$genus
#bird_subset_all <- bird_subset_all %>% ungroup() %>% filter(genus %in% pairedgenera)

scinames <- with(bird_allrecords, unique(paste(genus, specificepithet)))
write.csv(data.frame(n=scinames), file='~/verts/birdscinames_correct.csv', row.names=FALSE)

#######################################################################

# After making manual corrections to names, get information for all the correctly classified individuals.
# Use this to get the latitude to see if we can find true sister taxa within genera that are temperate + tropical.

phylonames <- read.csv('~/verts/phylomatches_correct.csv', stringsAsFactors = FALSE)
bird_allrecords <- filter(bird_allrecords, specificepithet != 'maculatus,  ocai') %>% mutate(sciname = paste(genus, specificepithet))

name_to_use <- phylonames$n
name_to_use[!phylonames$match] <- phylonames$jetz_name[!phylonames$match]

bird_allrecords <- left_join(bird_allrecords, data.frame(sciname = phylonames$n, name_to_use = name_to_use))

bird_trop <- bird_allrecords %>% ungroup %>% group_by(name_to_use) %>% summarize(avglat = median(lat, na.rm=TRUE), tropical = abs(avglat) < 23.5)

library(ape)

# Read 100 trees containing the 774 species from birdtree.org (Ericson trees, randomly chosen subset)

erictree <- read.nexus('~/verts/ericson100.tre')
t1 <- erictree[[1]]
t1$tip.label <- gsub('_', ' ', t1$tip.label)

# Find sister species pairs that are temperate and tropical
# There are only 734 rows in distance matrix because some of the species got lumped.

# Use the tropical member as species 1 of the pair, and find any temperate sisters.

edist <- cophenetic(t1)
epairs <- list()
for (i in 1:nrow(edist)) {
  x <- edist[i, edist[i, ] > 0]
  sisters <- names(x)[which(x == min(x))]
  sp1_istrop <- bird_trop$tropical[bird_trop$name_to_use == dimnames(edist)[[1]][i]]
  sisters_istrop <- logical(0)
  for (j in 1:length(sisters)) {
	sisters_istrop[j] <- bird_trop$tropical[bird_trop$name_to_use == sisters[j]]
  }
  sisters <- sisters[!sisters_istrop]
  if (sp1_istrop & length(sisters) > 0) {
  epairs[[length(epairs) + 1]] <- sisters
  names(epairs)[[length(epairs)]] <- dimnames(edist)[[1]][i]
  }
}

# Write the pairs out in data frame

epairs <- lapply(epairs, function(x) x[1:6])
epairs_df <- do.call('rbind',epairs)
epairs_df <- as.data.frame(cbind(tropical_bird = dimnames(epairs_df)[[1]], epairs_df), stringsAsFactors=FALSE)
names(epairs_df)[2:7] <- paste('temperate_sister', 1:6, sep='_')

write.csv(epairs_df, file='~/verts/sisterpairs_by_phylogeny.csv', row.names = FALSE)

# Make a phylogeny with just these sisters
use_spp <- unique(na.omit(unlist(epairs_df)))
drop_spp <- t1$tip.label[!t1$tip.label %in% use_spp]
t1_drop <- drop.tip(t1, drop_spp)
plot(t1_drop)

# Find only single pairs.
epairs_single <- epairs_df[apply(!is.na(epairs_df),1,sum) == 2, 1:2]
temperatenames <- names(table(epairs_single[,2]))[table(epairs_single[,2]) == 1]
epairs_single <- epairs_single[epairs_single[,2] %in% temperatenames,] # There are 75 true pairs.
epairs_single <- data.frame(tropical = as.character(epairs_single[,1]), temperate=as.character(epairs_single[,2]), stringsAsFactors=FALSE)

# Test these things.
pairs_df <- list()
for (i in 1:nrow(epairs_single)) {
	sp_t <- epairs_single$tropical[i]
	sp_nont <- epairs_single$temperate[i]
	dat_t <- subset(bird_allrecords, name_to_use == sp_t)
	dat_nont <- subset(bird_allrecords, name_to_use == sp_nont)
	logmass_t <- na.omit(dat_t$logmass)
	logmass_nont <- na.omit(dat_nont$logmass)
	n_use <- min(length(logmass_t), length(logmass_nont))
	samp_t <- sample(logmass_t, size=n_use, replace=FALSE)
	samp_nont <- sample(logmass_nont, size=n_use, replace=FALSE)
	pairs_df[[i]] <- data.frame(sp_t = sp_t, sp_nont = sp_nont, mean_t = mean(samp_t), mean_nont = mean(samp_nont), n_t = length(logmass_t), n_nont = length(logmass_nont), n_use = n_use, cv_t = sd(samp_t)/mean(samp_t), cv_nont = sd(samp_nont)/mean(samp_nont))
	
	
}
pairs_df <- do.call('rbind', pairs_df)

# Run t-test
with(pairs_df, t.test(cv_t, cv_nont, paired = TRUE, alternative = 'less')) 


# Test with many iterations
nsim <- 999
ttests <- list()
pb <- txtProgressBar(0,nsim,style=3)

for (sim in 1:nsim) {
pairs_df <- list()
for (i in 1:nrow(epairs_single)) {
	sp_t <- epairs_single$tropical[i]
	sp_nont <- epairs_single$temperate[i]
	dat_t <- subset(bird_allrecords, name_to_use == sp_t)
	dat_nont <- subset(bird_allrecords, name_to_use == sp_nont)
	logmass_t <- na.omit(dat_t$logmass)
	logmass_nont <- na.omit(dat_nont$logmass)
	n_use <- min(length(logmass_t), length(logmass_nont))
	samp_t <- sample(logmass_t, size=n_use, replace=FALSE)
	samp_nont <- sample(logmass_nont, size=n_use, replace=FALSE)
	pairs_df[[i]] <- data.frame(sp_t = sp_t, sp_nont = sp_nont, mean_t = mean(samp_t), mean_nont = mean(samp_nont), n_t = length(logmass_t), n_nont = length(logmass_nont), n_use = n_use, cv_t = sd(samp_t)/mean(samp_t), cv_nont = sd(samp_nont)/mean(samp_nont))
	
	
}
pairs_df <- do.call('rbind', pairs_df)

ttests[[sim]] <- with(pairs_df, t.test(cv_t, cv_nont, paired = TRUE, alternative = 'less')) 

setTxtProgressBar(pb, sim)
}

close(pb)

ps <- c(.5,.025,.975)
quantile(sapply(ttests, function(x) as.numeric(x$estimate)), prob=ps)
quantile(sapply(ttests, function(x) as.numeric(x$stat)), prob=ps)
quantile(sapply(ttests, function(x) as.numeric(x$p.val)), prob=ps)

# Do not do any subsampling.
pairs_df <- list()
for (i in 1:nrow(epairs_single)) {
	sp_t <- epairs_single$tropical[i]
	sp_nont <- epairs_single$temperate[i]
	dat_t <- subset(bird_allrecords, name_to_use == sp_t)
	dat_nont <- subset(bird_allrecords, name_to_use == sp_nont)
	logmass_t <- na.omit(dat_t$logmass)
	logmass_nont <- na.omit(dat_nont$logmass)
	#n_use <- min(length(logmass_t), length(logmass_nont))
	#samp_t <- sample(logmass_t, size=n_use, replace=FALSE)
	#samp_nont <- sample(logmass_nont, size=n_use, replace=FALSE)
	pairs_df[[i]] <- data.frame(sp_t = sp_t, sp_nont = sp_nont, mean_t = mean(logmass_t), mean_nont = mean(logmass_nont), n_t = length(logmass_t), n_nont = length(logmass_nont), cv_t = sd(logmass_t)/mean(logmass_t), cv_nont = sd(logmass_nont)/mean(logmass_nont))
	
	
}
pairs_df <- do.call('rbind', pairs_df)

# Run t-test
with(pairs_df, t.test(cv_t, cv_nont, paired = TRUE, alternative = 'less')) 


# Use all the pairs, not just the single pairs.
pairs_df_all <- list()
for (i in 1:nrow(epairs_df)) {
	if (!is.na(epairs_df$temperate_sister_1[i])) {
		for (j in 1:sum(!is.na(epairs_df[i,2:7]))) {
			sp_t <- epairs_df$tropical_bird[i]
			sp_nont <- epairs_df[i, j+1]
			dat_t <- subset(bird_allrecords, name_to_use == sp_t)
			dat_nont <- subset(bird_allrecords, name_to_use == sp_nont)
			logmass_t <- na.omit(dat_t$logmass)
			logmass_nont <- na.omit(dat_nont$logmass)
			n_use <- min(length(logmass_t), length(logmass_nont))
			samp_t <- sample(logmass_t, size=n_use, replace=FALSE)
			samp_nont <- sample(logmass_nont, size=n_use, replace=FALSE)
			pairs_df_all[[length(pairs_df_all) + 1]] <- data.frame(sp_t = sp_t, sp_nont = sp_nont, mean_t = mean(samp_t), mean_nont = mean(samp_nont), n_t = length(logmass_t), n_nont = length(logmass_nont), n_use = n_use, cv_t = sd(samp_t)/mean(samp_t), cv_nont = sd(samp_nont)/mean(samp_nont), lat_t = median(dat_t$lat, na.rm = TRUE), lat_nont = median(dat_nont$lat, na.rm = TRUE), stringsAsFactors=FALSE)
		}
	}
}
pairs_df_all <- do.call('rbind', pairs_df_all)
with(pairs_df_all, t.test(cv_t, cv_nont, paired = TRUE, alternative = 'less')) 

write.csv(pairs_df_all, file = '~/verts/cv_pairs_truesisters.csv', row.names=FALSE)

strsplit(pairs_df_all$sp_t, ' ')
# Use only the most extreme difference pair in each genus.
pairs_df_all <- pairs_df_all %>% mutate(genus = sapply(strsplit(sp_t, ' '), '[', 1), latdiff = abs(lat_t - lat_nont)) 
genus_latdiffs <- pairs_df_all %>% group_by(genus) %>% summarize(maxlatdiff = max(latdiff))

pairs_df_extremes <- pairs_df_all %>% left_join(genus_latdiffs) %>% filter(latdiff == maxlatdiff)
with(pairs_df_extremes, t.test(cv_t, cv_nont, paired = TRUE, alternative = 'less')) 
