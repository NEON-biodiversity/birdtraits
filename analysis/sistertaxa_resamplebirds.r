library(dplyr)

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


fp <- '/mnt/research/plz-lab/NEON/external_data/raw_external_data/vertnet'
vnbird <- read.csv(file.path(fp, 'vertnet_birds_reduced.csv'), stringsAsFactors = FALSE)	
vnbird$specificepithet[vnbird$genus == 'Agelaius' & grepl('phoe', vnbird$specificepithet)] <- 'phoeniceus'
vnbird$genus[vnbird$genus == 'Yuhina' & vnbird$specificepithet == 'zantholeuca'] <- 'Erpornis'
vnbird$specificepithet[vnbird$genus == 'Helmitheros' & vnbird$specificepithet == 'vermivorus'] <- 'vermivorum'


bird_subset_all <- get_sister_allrecords(vnbird)

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