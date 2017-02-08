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
  select(sister1, sister2, d, logarea1, logarea2, meanlogbodymass, spatial_cv_temp1, spatial_cv_temp2, interannual_var_temp1, interannual_var_temp2, lograngesize1, lograngesize2) %>%
  filter(complete.cases(.))

betterlm <- lm(d ~ ., data = lmdat2[,-(1:2)], na.action = 'na.pass')
dredge(betterlm)
betterlm2 <- lm(d ~ I(interannual_var_temp1 - interannual_var_temp2) + meanlogbodymass + I(spatial_cv_temp1 - spatial_cv_temp2) + I(logarea1 - logarea2) + I(lograngesize1 - lograngesize2), data = lmdat2)
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

pdf(file.path(fp, 'vertnet_results/covariatescatterplots.pdf'), height = 6, width = 6)

ggplot(lmdat2, aes(x = interannual_var_temp1 - interannual_var_temp2, y = d)) +
  geom_point() + stat_smooth(method='lm') +
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is correlated\nwith reduced interannual temperature variability') +
  labs(x = expression(paste('interannual ',Delta * CV[temperature])), y = expression(Delta * CV[bodymass]))

ggplot(lmdat2, aes(x = meanlogbodymass, y = d)) +
  geom_point() + stat_smooth(method='lm') +
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is correlated\nwith larger body sizes') +
  labs(x = expression(paste(log[10],' body mass')), y = expression(Delta * CV[bodymass]))

ggplot(lmdat2, aes(x = spatial_cv_temp1 - spatial_cv_temp2, y = d)) +
  geom_point() + 
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is NOT correlated\nwith reduced spatial temperature variability') +
  labs(x = expression(paste('spatial ',Delta * CV[temperature])), y = expression(Delta * CV[bodymass]))

ggplot(lmdat2, aes(x = logarea1 - logarea2, y = d)) +
  geom_point() + 
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is NOT correlated\nwith reduced area over which specimens were collected') +
  labs(x = expression(Delta * log[10]*area), y = expression(Delta * CV[bodymass]))

ggplot(lmdat2, aes(x = lograngesize1 - lograngesize2, y = d)) +
  geom_point() + 
  theme_minimal() + 
  ggtitle('Reduced variability in tropical species is NOT correlated\nwith reduced tropical range size') +
  labs(x = expression(Delta * log[10]*rangesize), y = expression(Delta * CV[bodymass]))

dev.off()