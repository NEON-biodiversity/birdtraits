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
  mutate(realm = 'tropical')

botw_filtered$realm[abs(botw_filtered$X2) > 23.5] <- 'nontropical'

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
vnsis <- left_join(vnsis,
                         botw_filtered %>% select(binomial, realm) %>% rename(binomial_corr = binomial))

# Get only temperate-tropical pairs from vnsis.
bird_realms <- vnsis %>%
	group_by(binomial_corr) %>%
	summarize(nrealm = length(unique(realm)))
	
# Get rid of species that have breeding and year-round range in multiple realms
tworealmspp <- bird_realms$binomial_corr[bird_realms$nrealm==2]

vnsis_filter <- vnsis %>%
	filter(!binomial_corr %in% tworealmspp)
	
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
sistertree_coll <- sistertree
for (i in 1:50) {
	mrca_node <- findMRCA(sistertree_coll, tips = c(tropical_species[i], nontropical_species[i]))
	sistertree_coll$edge.length[which(sistertree_coll$edge[,1] == mrca_node)] <- 0
}
sistertree_coll <- drop.tip(sistertree_coll, tip = nontropical_species)
sistertree_coll_reroot <- reroot(sistertree_coll, node.number = 1, position = 10)

nontr_dat <- sister_data$nontropical_cv_logmass
names(nontr_dat) <- sister_data$tropical_binomial
tr_dat <- sister_data$tropical_cv_logmass
names(tr_dat) <- sister_data$tropical_binomial

phyl.pairedttest(tree = sistertree_coll_reroot, x1 = nontr_dat, x2 = tr_dat)
phyl.pairedttest(tree = sistertree_coll, x1 = nontr_dat, x2 = tr_dat, lambda = 1, fixed = TRUE)
phyl.pairedttest(tree = drop.tip(sistertree, tip = nontropical_species), x1 = nontr_dat, x2 = tr_dat)
phyl.pairedttest(tree = drop.tip(sistertree, tip = nontropical_species), x1 = cbind(nontr_dat,tr_dat), lambda = 0, fixed = TRUE)

sistertree_justtrop <- drop.tip(sistertree, tip = nontropical_species)
nontr_dat <- nontr_dat[sistertree_justtrop$tip.label]
tr_dat <- tr_dat[sistertree_justtrop$tip.label]
phyl.pairedttest(tree = multi2di(sistertree_justtrop), x1=tr_dat, x2=nontr_dat)

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
quantile(sapply(ttests, function(x) as.numeric(x$estimate)), prob=ps) # 0.008022447 0.003522459 0.013900431
table(sapply(ttests, function(x) as.numeric(x$estimate) > 0)) # true in all cases.
quantile(sapply(ttests, function(x) as.numeric(x$stat)), prob=ps) # 1.7348267 0.8049204 2.7027429 
quantile(sapply(ttests, function(x) as.numeric(x$p.val)), prob=ps) # 0.044529173 0.004711851 0.212380702 
sum( sapply(ttests, function(x) as.numeric(x$p.val)) < 0.05)/999 #  0.5505506


#### Here insert covariate analysis.