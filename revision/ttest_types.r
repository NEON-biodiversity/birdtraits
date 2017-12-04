# Unpaired t test with only sister species

with(sister_data, t.test(nontropical_cv_logmass, tropical_cv_logmass, paired = FALSE, alternative = 'greater')) # t131.69 = 1.8071, p = 0.0365

# Unpaired t test with ALL species

allspp <- unique(vnbird$binomial)
table(allspp %in% birdtree$tip.label)

allcvs <- left_join(vnbird, botw_filtered %>% select(binomial, realm)) %>%
  filter(binomial %in% birdtree$tip.label) %>%
  group_by(binomial, realm) %>%
  summarize(cv_logmass = sd(log10(massing), na.rm = TRUE)/mean(log10(massing), na.rm = TRUE)) %>%
  filter(!is.na(realm))

tropical_cvs <- with(allcvs, setNames(cv_logmass[which(realm=='tropical')], binomial[which(realm=='tropical')]))
nontropical_cvs <- with(allcvs, setNames(cv_logmass[which(realm=='nontropical')], binomial[which(realm=='nontropical')]))

length(tropical_cvs); length(nontropical_cvs)

t.test(nontropical_cvs, tropical_cvs, paired = FALSE, alternative = 'greater') # t2186.7 = 8.3315, p < 2.2e-16

# Unpaired t test with only sister species, phylogenetically corrected

sistertree
unpaired_data <- with(sister_data, 
                      data.frame(species = c(nontropical_binomial, tropical_binomial),
                                 zone = rep(c('nontropical', 'tropical'), each = length(nontropical_binomial)),
                                 cvlogmass = c(nontropical_cv_logmass, tropical_cv_logmass)))

library(phylolm)
dimnames(unpaired_data)[[1]] <- unpaired_data$species
plm <- phylolm(cvlogmass ~ zone, data = unpaired_data, phy = sistertree, model = 'BM')
plm <- phylolm(cvlogmass ~ zone, data = unpaired_data, phy = sistertree, model = 'OUrandomRoot')

cor.test(picx, picy)

# Unpaired t test with ALL species, phylogenetically corrected
allcvs_df <- as.data.frame(allcvs)
dimnames(allcvs_df)[[1]] <- allcvs_df$binomial
birdtree_sub <- drop.tip(birdtree, tip = birdtree$tip.label[!birdtree$tip.label %in% allcvs$binomial])
plm_all <- phylolm(cv_logmass ~ realm, data = allcvs_df, phy = birdtree_sub, model = 'BM')
plm_all <- phylolm(cv_logmass ~ realm, data = allcvs_df, phy = birdtree_sub, model = 'OUrandomRoot')


# Paired t test with only sister species ***

# Paired t test with only sister species, phylogenetically corrected