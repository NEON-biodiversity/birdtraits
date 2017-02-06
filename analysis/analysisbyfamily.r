# Complete bird pair analysis. From scratch. The old code was buggy. This should be fully reproducible.

fp <- 'C:/Users/Q/Dropbox/projects/verts'

# Source functions.
source('analysis/sisterfns.r')

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

#################################################

# Remove all species with less than 5 mass measurements.

# Create genus + species scientific name to get rid of subspecies classification
vnbird <- vnbird %>% mutate(binomial = paste(genus, specificepithet))

nmass <- vnbird %>% group_by(binomial) %>% summarize(nmassmeas = sum(!is.na(massing)))
vnbird <- vnbird %>% left_join(nmass) %>% filter(nmassmeas >= 5)

# Get rid of species with only one in a family.
vnbird$familyclean <- sapply(strsplit(vnbird$family, ' '), '[', 1)

n_infamily <- vnbird %>% group_by(familyclean) %>% summarize(nspp = length(unique(binomial)))
vnbird <- vnbird %>% left_join(n_infamily) %>% filter(nspp > 1)

# Remove non georeferenced entries
vnbird <- filter(vnbird, !is.na(decimallatitude))

# Still need to get rid of a lot of things. Try to get rid of families that do not have at least one temperate and one tropical member.

latbysp <- vnbird %>% group_by(familyclean, binomial) %>% summarize(mlat = median(decimallatitude, na.rm=T))
latboth <- latbysp %>% summarize(hasboth = any(abs(mlat) > 23.5) & any(abs(mlat) < 23.5))

vnbird <- vnbird %>% left_join(latboth) %>% filter(hasboth)
# This is now down to 217k records with 4569 species.

##########################################################


# Do analysis on this.

# Get rid of bird species that are listed as being in more than one family.
nfampersp <- vnbird %>% group_by(binomial) %>% summarize(nfam = length(unique(familyclean)))
vnbirdclean <- vnbird %>% left_join(nfampersp) %>% filter(nfam == 1)

cv_all <- vnbirdclean %>% group_by(familyclean, binomial) %>% 
  summarize(latitude = abs(median(decimallatitude)),
            cv_logmass = sd(log10(massing))/mean(log10(massing)),
            n_logmass = length(massing))

library(lme4)

allfamlmer <- lmer(cv_logmass ~ latitude + (1|familyclean), data = cv_all %>% filter(n_logmass >= 10))
summary(allfamlmer)
confint(allfamlmer, method = 'boot', nsim = 999, parm = 'latitude')


library(MuMIn)
r.squaredGLMM(allfamlmer)


# Export everything in the cleaned up by-family thing so we can get a big tree with the ~3000 species.
# There are a few too many unfortunately.
write.csv(cv_all$binomial, file = file.path(fp, 'namesforphylobyfam.csv'), row.names = FALSE)


# Use the slimmed down version
# bird_subset_all <- get_sister_allrecords(vnbird)
# bird_pairs_subsample <- bird_subset_all %>% ungroup %>% group_by(genus) %>% do(get_cv_pairs_subsample(.)) %>% filter(!is.na(cv_logmass_tropical))
# bird_allrecords <- bird_subset_all %>% filter(genus %in% bird_pairs_subsample$genus) # Down to 415 species.
# 
# write.table(unique(bird_allrecords$name_to_use), file = file.path(fp, 'phylo_namelist_02feb.csv'), row.names=F)
