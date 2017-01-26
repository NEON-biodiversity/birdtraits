# Comparison of sister taxa methods.

source('analysis/sisterfns.r')
source('data_extraction/loadvnbird.r')

# Single subsample; most extreme pair in each genus by latitude
with(bird_pairs_subsample, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less')) 

# No subsample
bird_subset <- get_sister(vnbird)
bird_pairs <- get_cv_pairs(bird_subset) # 93 valid pairs 
with(bird_pairs, t.test(cv_logmass_tropical, cv_logmass_nontropical, paired = TRUE, alternative = 'less')) # Significant!!!!
