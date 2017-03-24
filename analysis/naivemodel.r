# Statistical tests for covariates (multiple regression model)

fp <- 'C:/Users/Q/Dropbox/projects/verts'
load(file.path(fp, 'sistercovariates.r'))

library(dplyr)

# Multiple regression model on difference with various predictors.

sister_alldat$d <- sister_alldat$cv1 - sister_alldat$cv2
lmdat <- sister_alldat %>% 
  mutate(logarea1 = log10(area1), logarea2 = log10(area2)) %>%
  select(sister1, sister2, d, songbird, spatial_spread1, spatial_spread2, migrant_status1, migrant_status2, range_size1, range_size2, logarea1, logarea2, meanlogbodymass, Diet.5Cat_trop, Diet.5Cat_temp, spatial_cv_temp1, spatial_cv_temp2, seasonal_var_temp1, seasonal_var_temp2, interannual_var_temp1, interannual_var_temp2) %>%
  filter(complete.cases(.))

biglm <- lm(d ~ songbird + spatial_spread1 + spatial_spread2 + migrant_status1 + migrant_status2 + range_size1 + range_size2 + logarea1 + logarea2 + meanlogbodymass + Diet.5Cat_trop + Diet.5Cat_temp + spatial_cv_temp1 + spatial_cv_temp2 + seasonal_var_temp1 + seasonal_var_temp2 + interannual_var_temp1 + interannual_var_temp2, data = lmdat, na.action = 'na.pass')

library(MuMIn)
bigdredge <- dredge(biglm, m.lim = c(0,5))
subset(bigdredge, delta < 5)

# Explanatory variables: area from which specimens were taken, interannual and seasonal variability in temperature, and the average body mass.

lmdat2 <- sister_alldat %>% 
  mutate(logarea1 = log10(area1+1), logarea2 = log10(area2+1), lograngesize1 = log10(range_size1), lograngesize2 = log10(range_size2)) %>%
  dplyr::select(sister1, sister2, d, logarea1, logarea2, meanlogbodymass, spatial_cv_temp1, spatial_cv_temp2, interannual_var_temp1, interannual_var_temp2, spatial_cv_precip1, spatial_cv_precip2, interannual_var_precip1, interannual_var_precip2, lograngesize1, lograngesize2, total_richness1, total_richness2, congener_richness1, congener_richness2, n_subpops1, n_subpops2, elevation_cv1, elevation_cv2, seasonal_var_temp1, seasonal_var_temp2, seasonal_var_precip1, seasonal_var_precip2) %>%
  filter(complete.cases(.))

betterlm <- lm(d ~ ., data = lmdat2[,-(1:2)], na.action = 'na.pass')
dredge(betterlm)
betterlm2 <- lm(d ~ I(interannual_var_temp1 - interannual_var_temp2) + I(interannual_var_precip1 - interannual_var_precip2) + meanlogbodymass + I(spatial_cv_temp1 - spatial_cv_temp2) + I(spatial_cv_precip1 - spatial_cv_precip2) + I(logarea1 - logarea2) + I(lograngesize1 - lograngesize2) + I(total_richness1 - total_richness2) + I(congener_richness1 - congener_richness2) + I(n_subpops1 - n_subpops2) + I(elevation_cv1 - elevation_cv2), data = lmdat2)
summary(betterlm2)
confint(betterlm2)


#################################

# Model for individual species

lmdatlong <- sister_pair_long %>%
  mutate(logmass = log10(BodyMass.Value), logarea = log10(area + 1), lograngesize = log10(range_size)) %>%
  select(taxon, realm, songbird, lat, cv_logmass, Diet.5Cat, logmass, spatial_cv_temp, seasonal_var_temp, interannual_var_temp, spatial_spread, logarea, migrant_status, lograngesize) %>%
  filter(complete.cases(.))

biglonglm <- lm(cv_logmass ~ songbird + spatial_spread + migrant_status + lograngesize + logarea + logmass + Diet.5Cat + spatial_cv_temp + seasonal_var_temp + interannual_var_temp, data = lmdatlong, na.action = 'na.pass')

longdredge <- dredge(biglonglm)
subset(longdredge, delta < 2)

longlm2 <- lm(cv_logmass ~ logarea + logmass + interannual_var_temp + songbird + Diet.5Cat, data = lmdatlong)
summary(longlm2)
# CV goes up when area of collection goes up, goes down with increasing body size, goes up with increasing interannual temperature variability, goes up with increasing range size, and is lower in songbirds.


################################

# Plot.

library(ggplot2)

# Annotation for plots showing where tropical is more variable.
hl <- geom_hline(yintercept = 0, color = 'indianred3')
vl <- geom_vline(xintercept = 0, color = 'indianred3')
#hlab <- annotate('text', x = c(-Inf, -Inf), y = c(Inf, -Inf), label = c('tropical more variable','temperate more variable'), alpha = 0.6)

pdf(file.path(fp, 'vertnet_results/covariatescatterplots.pdf'), height = 6, width = 6)

ggplot(lmdat2, aes(x = interannual_var_temp1 - interannual_var_temp2, y = d)) +
  hl + vl +
  geom_point() + stat_smooth(method='lm') +
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is correlated\nwith reduced interannual temperature variability') +
  labs(x = expression(atop(paste('interannual ',Delta * CV[temperature]), 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))

ggplot(lmdat2, aes(x = interannual_var_precip1 - interannual_var_precip2, y = d)) +
  hl + vl +
  geom_point() + 
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is NOT correlated\nwith reduced interannual precipitation variability') +
  labs(x = expression(atop(paste('interannual ',Delta * CV[rainfall]), 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))

ggplot(lmdat2, aes(x = seasonal_var_temp1 - seasonal_var_temp2, y = d)) +
  hl + vl +
  geom_point() + stat_smooth(method='lm') +
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is correlated\nwith reduced SEASONAL temperature variability') +
  labs(x = expression(atop(paste('seasonal ',Delta * CV[temperature]), 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))

ggplot(lmdat2, aes(x = seasonal_var_precip1 - seasonal_var_precip2, y = d)) +
  hl + vl +
  geom_point() + geom_point() + stat_smooth(method='lm') +
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is correlated\nwith INCREASED SEASONAL precipitation variability') +
  labs(x = expression(atop(paste('seasonal ',Delta * CV[rainfall]), 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))

ggplot(lmdat2, aes(x = meanlogbodymass, y = d)) +
  hl +
  geom_point() + stat_smooth(method='lm') +
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is correlated\nwith larger body sizes') +
  labs(x = expression(paste(log[10],' body mass')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))

p1 <-
ggplot(lmdat2, aes(x = spatial_cv_temp1 - spatial_cv_temp2, y = d)) +
  hl + vl +
  geom_point() + 
  theme_minimal() + 
  #ggtitle('Reduced variability in tropical species is NOT correlated\nwith reduced spatial temperature variability') +
  labs(x = expression(atop(paste('spatial ',Delta * CV[temperature]), 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))
#ggsave(file.path(fp, 'vertnet_results/temp_spatial_scatter.png'), height=6, width=6, dpi=400)

ggplot(lmdat2, aes(x = spatial_cv_precip1 - spatial_cv_precip2, y = d)) +
  hl + vl +
  geom_point() + 
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is NOT correlated\nwith reduced spatial precipitation variability') +
  labs(x = expression(atop(paste('spatial ',Delta * CV[rainfall]), 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))

p2 <-
ggplot(lmdat2, aes(x = logarea1 - logarea2, y = d)) +
  hl + vl +
  geom_point() + 
  theme_minimal() + 
  #ggtitle('Reduced variability in tropical species is NOT correlated\nwith reduced area over which specimens were collected') +
  labs(x = expression(atop(Delta * log[10]*area, 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))
#ggsave(file.path(fp, 'vertnet_results/area_collected_scatter.png'), height=6, width=6, dpi=400)

p3 <-
ggplot(lmdat2, aes(x = lograngesize1 - lograngesize2, y = d)) +
  hl + vl +
  geom_point() + 
  theme_minimal() + 
  #ggtitle('Reduced variability in tropical species is NOT correlated\nwith reduced tropical range size') +
  labs(x = expression(atop(Delta * log[10]*rangesize, 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))
#ggsave(file.path(fp, 'vertnet_results/range_size_scatter.png'), height=6, width=6, dpi=400)

p4 <-
ggplot(lmdat2, aes(x = total_richness1 - total_richness2, y = d)) +
  hl + vl +
  geom_point() + 
  theme_minimal() + 
  #ggtitle('Reduced variability in tropical species is NOT correlated\nwith increased tropical richness') +
  labs(x = expression(atop(Delta * richness, 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))
#ggsave(file.path(fp, 'vertnet_results/total_richness_scatter.png'), height=6, width=6, dpi=400)

p5 <-
ggplot(lmdat2, aes(x = congener_richness1 - congener_richness2, y = d)) +
  hl + vl +
  geom_point() + 
  theme_minimal() + 
  #ggtitle('Reduced variability in tropical species is NOT correlated\nwith increased tropical congener richness') +
  labs(x = expression(atop(Delta * congeners, 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))
#ggsave(file.path(fp, 'vertnet_results/congener_richness_scatter.png'), height=6, width=6, dpi=400)

p6 <-
ggplot(lmdat2, aes(x = n_subpops1 - n_subpops2, y = d)) +
  hl + vl +
  geom_point() + 
  theme_minimal() + 
  #ggtitle('Reduced variability in tropical species is NOT correlated\nwith number of subpopulations collected from') +
  labs(x = expression(atop(Delta * subpopulations, 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))
#ggsave(file.path(fp, 'vertnet_results/number_subpops_scatter.png'), height=6, width=6, dpi=400)

p7 <-
ggplot(lmdat2, aes(x = elevation_cv1 - elevation_cv2, y = d)) +
  hl + vl +
  geom_point() + 
  theme_minimal() + 
  #ggtitle('Reduced variability in tropical species is NOT correlated\nwith reduced variability in elevation') +
  labs(x = expression(atop(Delta * CV[elevation], 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass])))
#ggsave(file.path(fp, 'vertnet_results/elevation_range_scatter.png'), height=6, width=6, dpi=400)



dev.off()

library(gridExtra)
png(file.path(fp, 'vertnet_results/other_covariates_scatter.png'), height=15, width=8.5, res=400, units='in')
grid.arrange(grobs = list(p1,p2,p3,p4,p5,p6,p7), nrow=4)
dev.off()

###########################

# Does migrant status modify any of these relationships
sister_pair_long$migrantbinary <- sister_pair_long$migrant_status %in% c('obligate','partial')
ggplot(sister_pair_long, aes(x = interaction(realm, migrantbinary), y = cv_logmass)) + 
  geom_boxplot() + theme_minimal()
realmbymigrant <- lm(cv_logmass ~ realm * migrantbinary, data = sister_pair_long)
summary(realmbymigrant)


lmpairwise <- lm(d ~ migrant_status1 + migrant_status2, data = sister_alldat)
summary(lmpairwise)
