# Comparison of sister taxa methods.

source('analysis/sisterfns.r')
source('data_extraction/loadvnbird.r')

# Single subsample; most extreme pair in each genus by latitude
with(bird_pairs_subsample, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less')) 

# No subsample
bird_subset <- get_sister(vnbird)
bird_pairs <- get_cv_pairs(bird_subset) # 93 valid pairs 
with(bird_pairs, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less')) # Significant!!!!

# Test the hypothesis that latitude is related to itv.
bird_cvs <- bird_allrecords %>% summarize(lat = median(lat, na.rm=T),
                                          cv_logmass = sd(logmass, na.rm=T)/mean(logmass, na.rm=T))

library(lme4)
bird_lmer <- lmer(cv_logmass ~ lat + (1|genus), data = bird_cvs)
summary(bird_lmer)
confint(bird_lmer, method = 'boot', parm = 'lat', nsim = 999)
