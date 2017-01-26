# Functions for sister taxa

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
