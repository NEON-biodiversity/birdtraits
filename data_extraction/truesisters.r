# New bird analysis: Start with phylogeny and get the sister pairs first. Then get the traits.

library(ape)
allbirds <- read.tree('~/verts/trees/AllBirdsEricson1.tre')

# Select 50 trees at random to use.
set.seed(6880535)
idx <- sample(1:1000, 50)
allbirds <- allbirds[idx]

# Calculate a cophenetic distance matrix.
ericsondist <- lapply(allbirds, cophenetic)



# Must sort the distance matrices so that they are all the same order.
roworder <- dimnames(ericsondist[[1]])[[1]]
ericsondist <- lapply(ericsondist, function(x) x[roworder, roworder])
# Mean branch lengths across the ten randomly sampled trees
ericsondistmean <- Reduce('+', ericsondist)/length(ericsondist)
#ericsondistmean <- apply(simplify2array(ericsondist), 1:2, mean) # This overflows! :-O

# A true sister pair is defined as taxa A and B such that the nearest neighbor of taxon A is taxon B, and the nearest neighbor of taxon B is taxon A and only taxon A.

sisters <- character(nrow(ericsondistmean))

for (i in 1:nrow(ericsondistmean)) {
	x <- ericsondistmean[i, -i]
	name_i <- dimnames(ericsondistmean)[[1]][i]
	dist_to_sister <- min(x)
	sisternames <- names(x)[x == dist_to_sister] # Names of potential sisters.
	# Check which of these sisters has taxon i as a sister.
	for (j in 1:length(sisternames)) {
		j_index <- which(dimnames(ericsondistmean)[[1]] == sisternames[j])
		y <- ericsondistmean[j_index, -j_index]
		reciprocaldist_to_sister <- min(y)
		reciprocalsisternames <- names(y)[y == reciprocaldist_to_sister]
		if (name_i %in% reciprocalsisternames) sisters[[i]] <- sisternames[j]
	}
	
}

sisters_df <- data.frame(sister1 = dimnames(ericsondistmean)[[1]], sister2 = sisters)
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

sisters_df <- sisters_df[!is_dup,]

write.csv(sisters_df, file = '~/verts/truesisters02feb.csv', row.names = FALSE)


####################################################


# Match this with the vertnet.
fp <- 'C:/Users/Q/Dropbox/projects/verts'
sisters <- read.csv(file.path(fp, 'truesisters02feb.csv'), stringsAsFactors = FALSE)

library(dplyr)

######################################################

# Load data and Quality control

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


#################################################################

sister_names <- with(sisters, c(sister1, sister2))
vnbird <- vnbird %>% mutate(binomial = paste(genus, sapply(strsplit(specificepithet, ' '), '[', 1), sep = '_'))

# Get rid of anything less than 5 mass measurements
nmass <- vnbird %>% group_by(binomial) %>% summarize(nmassmeas = sum(!is.na(massing)))
vnbird <- vnbird %>% left_join(nmass) %>% filter(nmassmeas >= 10)


sistersindata <- sister_names %in% vnbird$binomial
#sister_names[!sistersindata]

vnbird$issister <- vnbird$binomial %in% sister_names
missingnames <- vnbird %>% filter(!issister) %>% group_by(binomial) %>% summarize(n = n()) %>% arrange(-n)

# Get the information for each of the pairs.
bird_summ <- vnbird %>% filter(issister) %>% group_by(binomial) %>% summarize(n = n(), lat = median(decimallatitude, na.rm=T), cv_logmass = sd(log10(massing))/mean(log10(massing)))

# Save the summary information.
vnsis <- vnbird %>% filter(issister)
#write.csv(with(vnsis, data.frame(binomial, lat=decimallatitude, lon=decimallongitude, massing=massing)), file = file.path(fp, 'bird_coords_mass.csv'), row.names = FALSE)

sister_join <- left_join(sisters, with(bird_summ, data.frame(sister1=binomial, n1=n, lat1=lat, cv1=cv_logmass)))
sister_join <- left_join(sister_join, with(bird_summ, data.frame(sister2=binomial, n2=n, lat2=lat, cv2=cv_logmass)))
sister_join <- sister_join[complete.cases(sister_join), ]

sister_join <- mutate(sister_join, lat1 = abs(lat1), lat2 = abs(lat2))

# Sort them so that lat1 is always the lower latitude.
for (i in 1:nrow(sister_join)) {
  if (sister_join$lat1[i] > sister_join$lat2[i]) {
    tmp <- sister_join[i,]
    sister_join[i, ] <- tmp[c(2,1,6,7,8,3,4,5)]
  }
}

sister_join <- mutate(sister_join, dlat = lat2 - lat1)

##################################################################

# Analysis 1: Use the temperate and tropical pairs.
# NOTE: It's possible we may have to remove Butorides because the two species were lumped at one point.

sister_troptemp <- filter(sister_join, lat1 < 23.5 & lat2 > 23.5)

with(sister_troptemp, t.test(cv1, cv2, paired = TRUE, alternative = 'less')) # Significant!

# write.csv(sister_troptemp, file = file.path(fp, 'sister_troptemp.csv'), row.names = FALSE)

# Analysis 2: Use only ones with at least 10 degrees latitude separation.

#sister_10 <- filter(sister_join, dlat > 10)
#with(sister_10, t.test(cv1, cv2, paired = TRUE, alternative = 'less')) # Significant!


#################################################################

# Analysis 1b: Use the temperate and tropical pairs, do a subsample, and do many iterations.

nsim <- 999
ttests <- list()
pb <- txtProgressBar(0, nsim, style = 3)

for (sim in 1:nsim) {
  
  sister_subsample <- sister_troptemp
  
  
  for (i in 1:nrow(sister_troptemp)) {
    if (sister_troptemp$n1[i] < sister_troptemp$n2[i]) {
      logmass2_subsample <- log10(sample(vnbird$massing[vnbird$binomial == sister_troptemp$sister2[i]], size = sister_troptemp$n1[i], replace = FALSE))
      sister_subsample$cv2[i] <- sd(logmass2_subsample)/mean(logmass2_subsample)
    }
    if (sister_troptemp$n2[i] < sister_troptemp$n1[i]) {
      logmass1_subsample <- log10(sample(vnbird$massing[vnbird$binomial == sister_troptemp$sister1[i]], size = sister_troptemp$n2[i], replace = FALSE))
      sister_subsample$cv1[i] <- sd(logmass1_subsample)/mean(logmass1_subsample)
    }
  }
  
  ttests[[sim]] <- with(sister_subsample, t.test(cv1, cv2, paired = TRUE, alternative = 'less'))
  setTxtProgressBar(pb, sim)
  
}

close(pb)

ps <- c(.5,.025,.975)
quantile(sapply(ttests, function(x) as.numeric(x$estimate)), prob=ps)
quantile(sapply(ttests, function(x) as.numeric(x$stat)), prob=ps)
quantile(sapply(ttests, function(x) as.numeric(x$p.val)), prob=ps)

# Figures of true sisters.


library(ggplot2)

# Convert paired data frame to longer form data frame.
sister_pair_long <-
  with(sister_troptemp, data.frame(taxon = c(sister1, sister2), 
                                   genus = sapply(strsplit(sister1, '_'), '[', 1), 
                                   realm = c(rep('tropical', nrow(sister_troptemp)), rep('nontropical', nrow(sister_troptemp))),
                                   lat = c(lat1, lat2),
                                   cv_logmass = c(cv1, cv2), 
                                   n = c(n1, n2)))

ylabel <- expression(paste('CV of log'[10], ' body mass'))							

pbird <- ggplot(sister_pair_long, aes(x = realm, y = cv_logmass)) +
  geom_line(aes(group = genus), color = 'gray75') + geom_point(aes(group = genus)) +
  stat_summary(aes(group = 1), geom = 'line', color = 'red', fun.y = 'mean', size = 1.5) +
  stat_summary(geom = 'pointrange', color = 'red', fun.data = 'mean_se') +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(y = ylabel)

phist <- ggplot(sister_troptemp, aes(x = cv1 - cv2)) +
  geom_histogram(bins = 15) +
  geom_vline(xintercept = -0.020, color = 'red', lwd = 2) +
  geom_vline(xintercept = -0.0092, color = 'red', lwd = 0.5) +
  geom_vline(xintercept = 0, color = 'blue', lty = 3) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,53), expand = c(0,0)) +
  labs(x = 'tropical CV - nontropical CV', y = 'count of species pairs')

ggsave(file.path(fp, 'vertnet_results/birdinteractionplot03feb.png'), pbird, height=4.5, width=4.5, dpi=300)
ggsave(file.path(fp, 'vertnet_results/birdhistogram03feb.png'), phist, height=4.5, width=4.5, dpi=300)


######################################################################

# Compare the CVs with different variables
trait_sisters <- read.csv(file.path(fp, 'trait_sisters.csv'), stringsAsFactors = FALSE)
sister_pair_long <- left_join(sister_pair_long, trait_sisters %>% rename(taxon = binomial))

# These don't seem to have any relationship.
plot(cv_logmass ~ Diet.Inv, data = sister_pair_long)
boxplot(cv_logmass ~ Diet.5Cat, data = sister_pair_long)
plot(cv_logmass ~ log10(BodyMass.Value), data = sister_pair_long)

# Look at whether the *difference* is due to some characteristic of the pair.

names(trait_sisters) <- paste(names(trait_sisters), 'trop', sep = '_')

sister_tt_traits <- sister_troptemp %>% left_join(trait_sisters %>% rename(sister1 = binomial_trop))

names(trait_sisters) <- gsub('trop', 'temp', names(trait_sisters))

sister_tt_traits <- sister_tt_traits %>% left_join(trait_sisters %>% rename(sister2 = binomial_temp))

library(ggplot2)

sister_tt_traits[sister_tt_traits == -999] <- NA

sister_tt_traits <- sister_tt_traits %>%
  mutate(meanlogbodymass = (log10(BodyMass.Value_temp) + log10(BodyMass.Value_temp))/2,
         meanlongevity = (maximum_longevity_y_trop + maximum_longevity_y_temp)/2
  )

ggplot(sister_tt_traits, aes(x = meanlogbodymass, y = cv1 - cv2)) +
  geom_point() + theme_minimal()
# Doesn't appear that the difference between temperate and tropical variability changes with body size.

ggplot(sister_tt_traits, aes(x = meanlongevity, y = cv1 - cv2)) +
  geom_point() + theme_minimal()
# Not enough of the species have a longevity for both members of the pair to say anything about that.

# Is the difference related to the functional group of the tropical species?
ggplot(sister_tt_traits, aes(x = Diet.5Cat_trop, y = cv1 - cv2)) +
  geom_boxplot() + theme_minimal()

ggplot(sister_tt_traits, aes(x = Diet.5Cat_temp, y = cv1 - cv2)) +
  geom_boxplot() + theme_minimal()

# Possibly the larger species that are carnivores tend to have a bigger discrepancy between the two.

###############################################################

# Extract order and family.'

foraging <- read.delim(file.path(fp, 'BirdFuncDat.txt'), stringsAsFactors = FALSE)
foraging$binomial <- gsub(' ', '_', foraging$Scientific)
taxgrps <- with(foraging, data.frame(sister1 = binomial, order = IOCOrder, family = BLFamilyLatin))
sister_troptemp <- left_join(sister_troptemp, taxgrps)

# Is it passerine or not?
sister_troptemp$songbird <- sister_troptemp$order == 'Passeriformes'
ggplot(sister_troptemp, aes(x = songbird, y = cv1 - cv2)) + geom_boxplot() + theme_minimal()
# It seems like the difference between tropical and temperate is bigger in non-songbirds.

###############################################################

# See whether the climate cv has anything to do with the actual cv.

birdrange_summstats <- read.csv(file.path(fp, 'bird_clim_summstats.csv'), stringsAsFactors = FALSE)
sister_pair_long <- left_join(sister_pair_long, birdrange_summstats %>% rename(taxon = binomial))

ggplot(sister_pair_long, aes(x = spatial_cv_temp, y = cv_logmass, color = realm)) +
  geom_point() + theme_minimal()
ggplot(sister_pair_long, aes(x = spatial_cv_precip, y = cv_logmass, color = realm)) +
  geom_point() + theme_minimal()
ggplot(sister_pair_long, aes(x = seasonal_var_temp, y = cv_logmass, color = realm)) +
  geom_point() + theme_minimal()
ggplot(sister_pair_long, aes(x = seasonal_var_precip, y = cv_logmass, color = realm)) +
  geom_point() + theme_minimal()
ggplot(sister_pair_long, aes(x = interannual_var_temp, y = cv_logmass, color = realm)) +
  geom_point() + theme_minimal()
ggplot(sister_pair_long, aes(x = interannual_var_precip, y = cv_logmass, color = realm)) +
  geom_point() + theme_minimal()
ggplot(sister_pair_long, aes(x = spatial_spread, y = cv_logmass, color = realm)) +
  geom_point() + theme_minimal()

# Is the difference in CVs in the species pair related to the difference in the different climate variables?
names(birdrange_summstats) <- paste(names(birdrange_summstats), '1', sep = '')
sister_tt_clim <- sister_troptemp %>% left_join(birdrange_summstats %>% rename(sister1 = binomial1))
names(birdrange_summstats) <- gsub('1', '2', names(birdrange_summstats))
sister_tt_clim <- sister_tt_clim %>% left_join(birdrange_summstats %>% rename(sister2 = binomial2))


ggplot(sister_tt_clim, aes(x = seasonal_var_temp1 - seasonal_var_temp2, y = cv1 - cv2)) +
  geom_point() + theme_minimal()
# The difference between the variabilities might not be due to the difference in seasonal variation

ggplot(sister_tt_clim, aes(x = interannual_var_temp1 - interannual_var_temp2, y = cv1 - cv2)) +
  geom_point() + theme_minimal()
ggplot(sister_tt_clim, aes(x = spatial_spread1 - spatial_spread2, y = cv1 - cv2)) +
  geom_point() + theme_minimal()


####################################################

# Load the range size and migrant status data and compare it.

botw_ranges <- read.csv(file.path(fp, 'botw_ranges.csv'), stringsAsFactors = F)

sister_pair_long <- left_join(sister_pair_long, botw_ranges)

names(botw_ranges) <- paste(names(botw_ranges), '1', sep = '')
sister_tt_range <- sister_troptemp %>% left_join(botw_ranges %>% rename(sister1 = taxon1))
names(botw_ranges) <- gsub('1', '2', names(botw_ranges))
sister_tt_range <- sister_tt_range %>% left_join(botw_ranges %>% rename(sister2 = taxon2))

# Replace one missing migrant status
sister_tt_range$migrant_status2[sister_tt_range$sister2 == 'Butorides_virescens'] <- 'partial'
sister_pair_long$migrant_status[sister_pair_long$taxon == 'Butorides_virescens'] <- 'partial'


# Look at effects of migrant status on CV.
with(sister_pair_long, table(realm, migrant_status))

ggplot(sister_pair_long, aes(x = migrant_status, y = cv_logmass, color = realm)) + 
  geom_jitter() + theme_minimal()
ggplot(sister_pair_long, aes(x = migrant_status, y = cv_logmass)) + 
  geom_boxplot() + theme_minimal()

# Look at effects of range size on CV
ggplot(sister_pair_long, aes(x = log10(range_size), y = cv_logmass, color = realm)) +
  geom_point() + theme_minimal()
ggplot(sister_tt_range, aes(x = log10(range_size1) - log10(range_size2), y = cv1 - cv2)) +
  geom_point() + theme_minimal()
ggplot(sister_pair_long, aes(x = realm, y = log10(range_size))) + geom_boxplot() + theme_minimal() # Might be different, but not really.

##############################################

# Save all the covariates together so that we can easily do analysis.

sister_pair_long <- left_join(sister_pair_long, taxgrps %>% rename(taxon = sister1))
sister_pair_long$songbird <- sister_pair_long$order == 'Passeriformes'

sister_alldat <- cbind(sister_tt_clim, sister_tt_range[,13:16], sister_tt_traits[,10:111])
save(sister_alldat, sister_pair_long, file = file.path(fp, 'sistercovariates.r'))
save(sister_alldat, sister_pair_long, vnsis, file = 'bird_sistertaxa_raw.RData')
