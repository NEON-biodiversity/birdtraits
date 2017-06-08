# All loading, processing, and plotting of data for bird trait main figure (and results, not including bootstraps)
# 27 April 2017
# 1 May 2017: multiple regression with VIF variable selection added.

library(dplyr)
library(cowplot)

fp <- 'C:/Users/Q/Dropbox/projects/verts'
sisters <- read.csv(file.path(fp, 'truesisters08jun.csv'), stringsAsFactors = FALSE)

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
indices <- pmin(match(vnbird$binomial, sisters$sister1), match(vnbird$binomial, sisters$sister2), na.rm=T)
vnbird$sisterdist <- sisters$dist[indices]
missingnames <- vnbird %>% filter(!issister) %>% group_by(binomial) %>% summarize(n = n()) %>% arrange(-n)

# Get the information for each of the pairs.
bird_summ <- vnbird %>% 
  filter(issister) %>% 
  group_by(binomial) %>% 
  summarize(n = n(), 
            lat = median(decimallatitude, na.rm=T),
            lon = median(decimallongitude, na.rm=T),
            cv_logmass = sd(log10(massing))/mean(log10(massing)),
            sisterdist = sisterdist[1])

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

# Get rid of orioles. They are an outlier when it comes to climate.
sister_join <- filter(sister_join, !grepl('Oriolus',sister1))

##################################################################

# Analysis 1: Use the temperate and tropical pairs.
# NOTE: It's possible we may have to remove Butorides because the two species were lumped at one point.

sister_troptemp <- filter(sister_join, lat1 < 23.5 & lat2 > 23.5)
with(sister_troptemp, t.test(cv2, cv1, paired = TRUE, alternative = 'greater')) # Significant!

# Output: t_99 = -3.115, p = 0.0012, cimax = -0.0096, d = -0.021

# Analysis 1b: Use the temperate and tropical pairs, do a subsample, and do many iterations.

nsim <- 999
ttests <- list()

set.seed(28462)

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
  
  ttests[[sim]] <- with(sister_subsample, t.test(cv2, cv1, paired = TRUE, alternative = 'greater'))
  setTxtProgressBar(pb, sim)
  
}

close(pb)

ps <- c(.5,.025,.975)
quantile(sapply(ttests, function(x) as.numeric(x$estimate)), prob=ps) #-0.01736805 -0.02457452 -0.01105836 
quantile(sapply(ttests, function(x) as.numeric(x$stat)), prob=ps) #-2.474535 -3.238879 -1.662695 
quantile(sapply(ttests, function(x) as.numeric(x$p.val)), prob=ps) #0.0075206051 0.0008169745 0.0497684158
sum( sapply(ttests, function(x) as.numeric(x$p.val)) < 0.05)/999 # 97.5%

# Convert paired data frame to longer form data frame.
sister_pair_long <-
  with(sister_troptemp, data.frame(taxon = c(sister1, sister2), 
                                   genus = sapply(strsplit(sister1, '_'), '[', 1), 
                                   realm = c(rep('tropical', nrow(sister_troptemp)), rep('nontropical', nrow(sister_troptemp))),
                                   lat = c(lat1, lat2),
                                   cv_logmass = c(cv1, cv2), 
                                   n = c(n1, n2)))


# FIGURES -----------------------------------------------------------------

load(file.path(fp, 'sistercovariates.r'))
sister_alldat$d <- sister_alldat$cv2 - sister_alldat$cv1
sister_alldat <- filter(sister_alldat, !grepl('Oriolus', sister1))
sister_pair_long <- filter(sister_pair_long, !grepl('Oriolus', genus))

lmdat <- sister_alldat %>% 
  mutate(logarea1 = log10(area1), logarea2 = log10(area2)) %>%
  dplyr::select(sister1, sister2, d, songbird, spatial_spread1, spatial_spread2, migrant_status1, migrant_status2, range_size1, range_size2, logarea1, logarea2, meanlogbodymass, Diet.5Cat_trop, Diet.5Cat_temp, spatial_cv_temp1, spatial_cv_temp2, seasonal_var_temp1, seasonal_var_temp2, interannual_var_temp1, interannual_var_temp2) %>%
  filter(complete.cases(.))
lmdat2 <- sister_alldat %>% 
  mutate(logarea1 = log10(area1+1), logarea2 = log10(area2+1), lograngesize1 = log10(range_size1), lograngesize2 = log10(range_size2)) %>%
  dplyr::select(sister1, sister2, d, logarea1, logarea2, meanlogbodymass, spatial_cv_temp1, spatial_cv_temp2, interannual_var_temp1, interannual_var_temp2, spatial_cv_precip1, spatial_cv_precip2, interannual_var_precip1, interannual_var_precip2, lograngesize1, lograngesize2, total_richness1, total_richness2, congener_richness1, congener_richness2, n_subpops1, n_subpops2, elevation_cv1, elevation_cv2, seasonal_var_temp1, seasonal_var_temp2, seasonal_var_precip1, seasonal_var_precip2) %>%
  filter(complete.cases(.))

mapdat <- read.csv(file.path(fp, 'mapdata.csv'))
sister_pair_long <- left_join(sister_pair_long, mapdat %>% rename(taxon=binomial) %>% select(taxon,lon))

mapdat$lat_temp <- mapdat$lat_trop <- mapdat$lat
mapdat$lon_temp <- mapdat$lon_trop <- mapdat$lon
mapdat$lat_temp[abs(mapdat$lat_temp) < 23.5] <- NA
mapdat$lon_temp[abs(mapdat$lat_temp) < 23.5] <- NA
mapdat$lat_trop[abs(mapdat$lat_trop) > 23.5] <- NA
mapdat$lon_trop[abs(mapdat$lat_trop) > 23.5] <- NA

sister_alldat <- left_join(sister_alldat, 
                           mapdat %>% rename(sister1=binomial) %>% dplyr::select(sister1, lat_trop,lon_trop))
sister_alldat <- left_join(sister_alldat, 
                           mapdat %>% rename(sister2=binomial) %>% dplyr::select(sister2, lat_temp,lon_temp))

############################################
# World Map

heatramp <- colorRampPalette(RColorBrewer::brewer.pal(name='YlOrRd', n=9),bias=2,space="rgb")(50)
fillScale <- scale_fill_gradientn(colours = heatramp, name = 'Body mass variability')
colorScale <- scale_color_manual(values = c('indianred','forestgreen')) # Line colors
linetypeScale <- scale_linetype_manual(name = '', values = c('dotted','solid'), labels = c('Tropical\nmore variable','Nontropical\nmore variable'))

worldMap <- borders('world', fill='gray75', color='black')
p_map <- ggplot() + worldMap + 
  geom_segment(aes(x=lon_temp, y=lat_temp, xend=lon_trop, yend=lat_trop, linetype = d>0), data=sister_alldat) +
  geom_point(aes(x=lon_trop,y=lat_trop, fill=cv1), data=sister_alldat, pch=21) +
  geom_point(aes(x=lon_temp,y=lat_temp, fill=cv2), data=sister_alldat, pch=21) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(breaks = c(-23.5, 23.5), labels = c('23.5° S', '23.5° N'), expand=c(0,0)) +
  fillScale + linetypeScale +
  theme(panel.background = element_rect(fill='skyblue'),
        panel.grid.major.y = element_line(color='black', size=1.1),
        panel.grid.major.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        legend.position=c(0,0), 
        legend.direction = 'horizontal',
        legend.justification = c(0,0)) + 
  coord_map('albers',0,0, ylim = c(-60,60))


#############################################
# Interaction and Histogram

ylabel <- expression(paste('CV of log'[10], ' body mass'))		

plotdat <- data.frame(realm = rep(c('tropical','nontropical'), each=nrow(sister_troptemp)),
                      taxon = c(sister_troptemp$sister1, sister_troptemp$sister2),
                      pairid = 1:nrow(sister_troptemp),
                      cv_logmass = c(sister_troptemp$cv1, sister_troptemp$cv2))

p_int <- ggplot(plotdat, aes(x = realm, y = cv_logmass)) +
  geom_line(aes(group = pairid), color = 'gray75') + 
  geom_point(aes(group = pairid)) +
  stat_summary(aes(group = 1), geom = 'line', color = 'red', fun.y = 'mean', size = 1.5) +
  stat_summary(geom = 'pointrange', color = 'red', fun.data = 'mean_se') +
  scale_x_discrete(expand=c(0.1,0.1), name='Realm') +
  labs(y = ylabel) +
  panel_border(colour='black')

phist <- ggplot(sister_troptemp, aes(x = cv1 - cv2)) +
  geom_histogram(bins = 15) +
  geom_vline(xintercept = -0.020, color = 'red', lwd = 2) +
  geom_vline(xintercept = -0.0092, color = 'red', lwd = 0.5) +
  geom_vline(xintercept = 0, color = 'blue', lty = 3) +
  scale_y_continuous(limits = c(0,53), expand = c(0,0)) +
  labs(x = 'tropical CV - nontropical CV', y = 'count of species pairs') +
  panel_border(colour='black')

##########################################
# Covariates

hl <- geom_hline(yintercept = 0, color = 'indianred3')
vl <- geom_vline(xintercept = 0, color = 'indianred3')

lmtemp <- lm(d ~ I(seasonal_var_temp2 - seasonal_var_temp1), data=lmdat2)
lmsize <- lm(d ~ meanlogbodymass, data=lmdat2)

pcov1 <- ggplot(lmdat2 %>% mutate(ivt = interannual_var_temp2 - interannual_var_temp1), aes(x = ivt, y = d)) +
  hl + 
  geom_point() + stat_smooth(method='lm', se=F) +
  panel_border(colour='black') + 
  labs(x = expression(paste('yearly ',Delta * CV[temperature])), y = expression(Delta * CV[bodymass])) +
  geom_text(data=data.frame(ivt=c(Inf, -2.5), d = c(Inf, -Inf), lab = c('Nontropical\nmore variable', 'Tropical\nmore variable')),
            aes(label=lab), hjust = c(1,0), vjust = c(1,-.5))
  #ggtitle('Reduced variability in tropical species is correlated\nwith reduced SEASONAL temperature variability') +
 # labs(x = expression(atop(paste('seasonal ',Delta * CV[temperature]), 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass]))) 

pcov2 <- ggplot(lmdat2 %>% mutate(bmass = 10^meanlogbodymass), aes(x = bmass, y = d)) +
  hl + scale_x_log10() +
  geom_point() + stat_smooth(method='lm', se=F) +
  panel_border(colour='black') + 
  labs(x = 'Body mass (g)', y = expression(Delta * CV[bodymass])) +
  geom_text(data=data.frame(bmass=c(Inf, 3), d = c(Inf, -Inf), lab = c('Nontropical\nmore variable', 'Tropical\nmore variable')),
            aes(label=lab), hjust = c(1,0), vjust = c(1,-.5))
  #ggtitle('Reduced variability in tropical species is correlated\nwith larger body sizes') +
  #labs(x = 'Body mass (g)', y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass]))) 

#########################################
# Assemble Plot

smalllab <- theme(axis.title.x=element_text(size=10),
                  axis.title.y=element_text(size=10))

bottom_row <- plot_grid(p_int, pcov1, pcov2, labels = c('b', 'c', 'd'), align = 'h', rel_widths = c(1, 1.4, 1.4), ncol=3)
full_plot <- plot_grid(p_map + panel_border(colour='black') + theme(legend.position='bottom'), bottom_row, labels = c('a', ''), ncol = 1, rel_heights = c(1, 0.9))
ggsave(file.path(fp, 'vertnet_results/multipanelplot1.png'), full_plot, height=6, width=8, dpi=400)

##########################################
# 1 May: multiple regression

multregdat <- lmdat2 %>% mutate(spatial_temp = spatial_cv_temp2-spatial_cv_temp1,
                                interann_temp = interannual_var_temp2-interannual_var_temp1,
                                spatial_precip = spatial_cv_precip2-spatial_cv_precip1,
                                interann_precip = interannual_var_precip2-interannual_var_precip1,
                                rangesize = lograngesize2-lograngesize1,
                                total_richness = total_richness2 - total_richness1,
                                congener_richness = congener_richness2 - congener_richness1,
                                n_pop = n_subpops2 - n_subpops1,
                                elev_cv = elevation_cv2 - elevation_cv1,
                                seasonal_temp = seasonal_var_temp2 - seasonal_var_temp1,
                                seasonal_precip = seasonal_var_precip2 - seasonal_var_precip1,
                                area = logarea2-logarea1) %>%
  dplyr::select(d, spatial_temp, interann_temp, spatial_precip, interann_precip, rangesize, total_richness, congener_richness, n_pop, elev_cv, seasonal_temp, seasonal_precip, area, meanlogbodymass)

library(usdm)
vif(multregdat[,-1])
cor(multregdat[,-1])

multregdat2 <- multregdat %>% dplyr::select(-interann_temp, -interann_precip)
multreg <- lm(d ~ ., data = multregdat2)

multregall <- lm(d ~ ., data = multregdat, na.action = 'na.pass')
library(MuMIn)
multregdredge <- dredge(multregall)
subset(multregdredge, delta < 2)

multregbest <- lm(d ~ interann_temp + meanlogbodymass + seasonal_precip + spatial_temp, data = multregdat)
summary(multregbest)

# Check the spatial temp one
summary(lm(d ~ spatial_temp, data = multregdat))
###########################################
# 1 May: supplemental figures
# Edited 6 June: clean up figure.

library(cowplot)

xaxisvars <- c('spatial_temp','interann_temp','seasonal_temp','spatial_precip','interann_precip','seasonal_precip','rangesize','total_richness','congener_richness','n_pop','elev_cv','area')
xaxislabels <- c('Spatial CV of MAT', 'Interannual CV of MAT', 'Seasonal CV of MAT', 'Spatial CV of MAP', 'Interannual CV of MAP', 'Seasonal CV of MAP', 'Range size (log10 km2)', 'Total richness', 'Congener richness', 'Populations collected', 'CV of elevation', 'Collection area (log10 km2)')

scatterplots <- list()

for (i in 1:length(xaxislabels)) {
  p_i <- ggplot(multregdat, aes_string(x = xaxisvars[i], y = 'd')) +
    geom_point() +
    panel_border(colour='black') +
    labs(y = expression(Delta*CV[bodymass]), x = xaxislabels[i])
  scatterplots[[i]] <- p_i
}

p_all <- plot_grid(plotlist = scatterplots, ncol = 3, labels = letters[1:length(xaxislabels)])
ggsave('C:/Users/Q/Dropbox/projects/verts/vertnet_results/supplementalscatterplots.png', p_all, height = 15*.85, width = 11*.85, dpi = 400)
