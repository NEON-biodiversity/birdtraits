# All code for bird analysis (once sisters are identified)
# Cleaned up 22 June 2017
# QDR

fp <- 'C:/Users/Q/Dropbox/projects/verts'
library(dplyr)

####################################################################################################
# Load data that has been flagged for bad values and outliers that are due to typos, and clean it up.

vnsis_flag <- read.csv('C:/Users/Q/Dropbox/projects/verts/bird_records_flagged.csv', stringsAsFactors = FALSE)
table(vnsis_flag$flag)

# Remove invalid, not_native, and not_wild. If typo or absent, use the corrected numbers I put in. If listed as opposite sign, change the sign and replace.

vnsis_flag <- filter(vnsis_flag, !flag %in% c('invalid','not_native','not_wild'))
latflag1 <- vnsis_flag$flag %in% c('lat_typo','latlong_typo','latlong_absent')
latflag2 <- vnsis_flag$flag %in% c('lat_opp','latlong_opp')
vnsis_flag$decimallatitude[latflag1] <- as.numeric(vnsis_flag$correct_lat[latflag1])
vnsis_flag$decimallatitude[latflag2] <- -vnsis_flag$decimallatitude[latflag2]

lonflag1 <- vnsis_flag$flag %in% c('long_typo','latlong_typo','latlong_absent')
lonflag2 <- vnsis_flag$flag %in% c('long_opp','latlong_opp')
vnsis_flag$decimallongitude[lonflag1] <- as.numeric(vnsis_flag$correct_lon[lonflag1])
vnsis_flag$decimallongitude[lonflag2] <- -vnsis_flag$decimallongitude[lonflag2]

# 22 June: flag those records that are greater than or less than 10x the median value for the species.
# 10 times is a lot for things to be varying and probably indicates a mistake.

vnfactor <- vnsis_flag %>% group_by(binomial) %>%
  summarize(median_mass = median(massing, na.rm = TRUE))
vnsis_flag2 <- vnsis_flag %>% left_join(vnfactor) %>% mutate(outlierflag = massing > median_mass * 10 | massing < median_mass / 10)
table(vnsis_flag2$outlierflag) # 1833 outliers flagged if 5x, 1178 flagged if 10x


####################################################################################################
# Join the data with the sister list and reorganize by pairs.

sister_summary <- read.csv(file.path(fp, 'sister_summary.csv'), stringsAsFactors = FALSE)

# Get only those that are supported in a good percentage of the trees.
sister_goodpairs <- subset(sister_summary, proportion >= 0.1)
length(unique(sister_goodpairs$sister2))

# Remove duplicates
sister_goodpairs$duplicate <- FALSE
for (i in 2:nrow(sister_goodpairs)) {
  previous <- c(sister_goodpairs$sister1[1:(i-1)], sister_goodpairs$sister2[1:(i-1)])
  if (sister_goodpairs$sister1[i] %in% previous | sister_goodpairs$sister2[i] %in% previous) sister_goodpairs$duplicate[i] <- TRUE
}

sister_goodpairs <- subset(sister_goodpairs, !duplicate)

bird_summ <- vnsis_flag2 %>% 
  filter(!outlierflag) %>% 
  group_by(binomial) %>% 
  summarize(n = n(), 
            lat = median(abs(decimallatitude), na.rm=T),
            lon = median(decimallongitude, na.rm=T),
            cv_logmass = sd(log10(massing))/mean(log10(massing)))

sister_join <- left_join(sister_goodpairs, with(bird_summ, data.frame(sister1=binomial, n1=n, lat1=lat, lon1=lon, cv1=cv_logmass)))
sister_join <- left_join(sister_join, with(bird_summ, data.frame(sister2=binomial, n2=n, lat2=lat, lon2=lon, cv2=cv_logmass)))
sister_join <- sister_join[complete.cases(sister_join), ]

# Sort them so that lat1 is always the lower latitude.
for (i in 1:nrow(sister_join)) {
  if (sister_join$lat1[i] > sister_join$lat2[i]) {
    tmp <- sister_join[i,]
    sister_join[i, c('sister1','n1','lat1','lon1','cv1')] <- tmp[c('sister2','n2','lat2','lon2','cv2')]
    sister_join[i, c('sister2','n2','lat2','lon2','cv2')] <- tmp[c('sister1','n1','lat1','lon1','cv1')]
  }
}

sister_join <- mutate(sister_join, dlat = lat2 - lat1, dcv = cv2 - cv1)

####################################################################################################
# Use only tropical-nontropical pairs and run t-test.

# Get rid of orioles. They are an outlier when it comes to climate.
sister_join <- filter(sister_join, !grepl('Oriolus',sister1))

# With all quality controls in place, we now have 95 sister pairs.
sister_troptemp <- filter(sister_join, lat1 < 23.5 & lat2 > 23.5)

with(sister_troptemp, t.test(cv2, cv1, paired = TRUE, alternative = 'greater')) 


####################################################################################################
# Run bootstrapped t-tests.

nsim <- 999
ttests <- list()

vnbird_good <- vnsis_flag2 %>% filter(!outlierflag)

set.seed(28462)

pb <- txtProgressBar(0, nsim, style = 3)

for (sim in 1:nsim) {
  
  sister_subsample <- sister_troptemp
  
  
  for (i in 1:nrow(sister_troptemp)) {
    if (sister_troptemp$n1[i] < sister_troptemp$n2[i]) {
      logmass2_subsample <- log10(sample(vnbird_good$massing[vnbird_good$binomial == sister_troptemp$sister2[i]], size = sister_troptemp$n1[i], replace = FALSE))
      sister_subsample$cv2[i] <- sd(logmass2_subsample)/mean(logmass2_subsample)
    }
    if (sister_troptemp$n2[i] < sister_troptemp$n1[i]) {
      logmass1_subsample <- log10(sample(vnbird_good$massing[vnbird_good$binomial == sister_troptemp$sister1[i]], size = sister_troptemp$n2[i], replace = FALSE))
      sister_subsample$cv1[i] <- sd(logmass1_subsample)/mean(logmass1_subsample)
    }
  }
  
  ttests[[sim]] <- with(sister_subsample, t.test(cv2, cv1, paired = TRUE, alternative = 'greater'))
  setTxtProgressBar(pb, sim)
  
}

close(pb)

ps <- c(.5,.025,.975)
quantile(sapply(ttests, function(x) as.numeric(x$estimate)), prob=ps) # 0.006610572 0.003098649 0.010572700 
table(sapply(ttests, function(x) as.numeric(x$estimate) > 0)) # true in all cases.
quantile(sapply(ttests, function(x) as.numeric(x$stat)), prob=ps) # 1.8915261 0.8951932 2.7166665 
quantile(sapply(ttests, function(x) as.numeric(x$p.val)), prob=ps) # 0.031157276 0.004070005 0.186736880
sum( sapply(ttests, function(x) as.numeric(x$p.val)) < 0.05)/999 #  0.6896897

####################################################################################################
# Gather all covariates and run multiple regression with them. (include the new covariates)

# Calculation of geographical distance within sister pairs.
deg2rad <- function(deg) return(deg*pi/180)

# Great circle distance: https://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.slc <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # Distance in km
}

sister_troptemp <- sister_troptemp %>%
  mutate(dist_slc = gcd.slc(lon1, lat1, lon2, lat2))

# Convert paired data frame to longer form data frame.
sister_pair_long <-
  with(sister_troptemp, data.frame(taxon = c(sister1, sister2), 
                                   genus = sapply(strsplit(sister1, '_'), '[', 1), 
                                   realm = c(rep('tropical', nrow(sister_troptemp)), rep('nontropical', nrow(sister_troptemp))),
                                   lat = c(lat1, lat2),
                                   cv_logmass = c(cv1, cv2), 
                                   n = c(n1, n2),
                                   dist_phy = meandist,
                                   dist_geo = dist_slc))

write.csv(sister_troptemp, file = file.path(fp,'sister_troptemp_23Jun.csv'), row.names = FALSE)

#########################################
load(file.path(fp, 'sistercovariates.r'))

# Recreate covariate data frame with the new reduced list of sister pairs.
# On second thought it is probably better to go back and recalculate them all since some are dependent on coordinates that I corrected.
sister_alldat2 <- subset(sister_alldat, sister1 %in% sister_troptemp$sister1)
sister_troptemp$sister2 %in% sister_alldat2$sister2
sister_troptemp <- sister_troptemp[order(sister_troptemp$sister1),]
sister_alldat2 <- sister_alldat2[order(sister_alldat2$sister1),]

all(sister_troptemp$sister1==sister_alldat2$sister1)
all(sister_troptemp$sister2==sister_alldat2$sister2)

sister_alldat_final <- cbind(sister_troptemp, sister_alldat2)

####################################################################################################
# Plot map.

####################################################################################################
# Plot interaction plot.

####################################################################################################
# Plot scatterplots of covariates that go into the main text.

####################################################################################################
# Plot scatterplots of covariates that go into the supplement.