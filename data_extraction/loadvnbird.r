# Load vn bird locally.
library(dplyr)

fp <- 'C:/Users/Q/Dropbox/projects/verts'
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

bird_subset_all <- get_sister_allrecords(vnbird)
bird_pairs_subsample <- bird_subset_all %>% ungroup %>% group_by(genus) %>% do(get_cv_pairs_subsample(.)) %>% filter(!is.na(cv_logmass_tropical))
bird_allrecords <- bird_subset_all %>% filter(genus %in% bird_pairs_subsample$genus) # Down to 774 species.

phylonames <- read.csv(file.path(fp, '/phylomatches_correct.csv'), stringsAsFactors = FALSE)
# Remove the one hybrid towhee species from the list.
bird_allrecords <- filter(bird_allrecords, specificepithet != 'maculatus,  ocai') %>% mutate(sciname = paste(genus, specificepithet))

name_to_use <- phylonames$n
name_to_use[!phylonames$match] <- phylonames$jetz_name[!phylonames$match]

bird_allrecords <- left_join(bird_allrecords, data.frame(sciname = phylonames$n, name_to_use = name_to_use))

bird_trop <- bird_allrecords %>% ungroup %>% group_by(name_to_use) %>% summarize(avglat = median(lat, na.rm=TRUE), tropical = abs(avglat) < 23.5)

library(ape)

# Read 100 trees containing the 774 species from birdtree.org (Ericson trees, randomly chosen subset)

erictree <- read.nexus(file.path(fp, 'ericson100.tre'))
t1 <- erictree[[1]]
t1$tip.label <- gsub('_', ' ', t1$tip.label)
