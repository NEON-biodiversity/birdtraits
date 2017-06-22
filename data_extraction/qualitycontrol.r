# Additional quality control of vertnet.

# Check years of specimen collection.
yrs <- as.numeric(vnsis$year)
goodyrs <- yrs[yrs > 1800 & yrs <= 2016]
hist(goodyrs)
table(goodyrs)

# Correspondence of records with climate data that we have.
y1 <- table(goodyrs >= 1979 & goodyrs <= 2013)
y2 <- table(goodyrs >= 1950 & goodyrs <= 2016)
y1/length(goodyrs)
y2/length(goodyrs)

# Check coordinates of vertnet records to see if they are OK. Do this first by plotting them and then see if there are any ways to speed up the process.
# Need to check if there are any valid rows, also need to check if the points are in the middle of the ocean or just add a bigger buffer.

library(cowplot)
library(dplyr)
spnames <- unique(vnsis$binomial)

write.csv(spnames, file.path(fp, 'allmaps_pagenumber.csv'))

pb <- txtProgressBar(0, length(spnames), style = 3)
count <- 0
pdf('~/figs/allmaps_vertnetcoords.pdf')
  for (i in spnames) {
    count <- count + 1
    setTxtProgressBar(pb, count)
    data_i <- vnsis %>% filter(binomial == i) %>% select(decimallongitude, decimallatitude) %>% filter(complete.cases(.))
    if (nrow(data_i) > 0) {
      xlim_i <- range(data_i$decimallongitude) + c(-5,5)
      ylim_i <- range(data_i$decimallatitude) + c(-5,5)
      p <- ggplot(subset(vnsis, binomial == i), aes(x = decimallongitude, y = decimallatitude)) +
        borders('world', xlim=xlim_i, ylim=ylim_i) +
        coord_map() +
        geom_point(color = 'red') +
        ggtitle(i)
    }
    else {
      p <- ggplot(data.frame(x=1,y=1,label='no data')) + geom_text(aes(x,y,label=label)) + ggtitle(i)
    }
    print(p)
  }
dev.off()
close(pb)

####
# 21 June: after flagging records, find any more bad outliers.
vnsis_flag <- read.csv('C:/Users/Q/Dropbox/projects/verts/bird_records_flagged.csv', stringsAsFactors = FALSE)

library(dplyr)
vnsd <- vnsis_flag %>% group_by(binomial) %>% 
  summarize(sd_mass = sd(massing), 
            mass_hi = mean(massing) + 3*sd_mass, 
            mass_lo = mean(massing) - 3*sd_mass)
vnsis_flag <- vnsis_flag %>% left_join(vnsd) %>% mutate(outlierflag = massing > mass_hi | massing < mass_lo)
table(vnsis_flag$outlierflag) # approx 1780 "outliers"

write.csv(vnsis_flag, file = 'C:/Users/Q/Dropbox/projects/verts/bird_records_flagged.csv', row.names = FALSE)

####
# 22 June: correct the flagged lat and long entries.

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

# 22 June: flag those records that are greater than or less than 5x the median value for the species.
# 5 times is a lot for things to be varying and probably indicates a mistake.

vnfactor <- vnsis_flag %>% group_by(binomial) %>%
  summarize(median_mass = median(massing, na.rm = TRUE))
vnsis_flag2 <- vnsis_flag %>% left_join(vnfactor) %>% mutate(outlierflag = massing > median_mass * 10 | massing < median_mass / 10)
table(vnsis_flag2$outlierflag) # 1833 outliers flagged if 5x, 1178 flagged if 10x

# Run analysis without the outliers.
bird_summ <- vnsis_flag2 %>% 
  filter(!outlierflag) %>% 
  group_by(binomial) %>% 
  summarize(n = n(), 
            lat = median(abs(decimallatitude), na.rm=T),
            lon = median(decimallongitude, na.rm=T),
            cv_logmass = sd(log10(massing))/mean(log10(massing)),
            sisterdist = sisterdist[1])

fp <- 'C:/Users/Q/Dropbox/projects/verts'
sisters <- read.csv(file.path(fp, 'truesisters08jun.csv'), stringsAsFactors = FALSE)

sister_join <- left_join(sisters, with(bird_summ, data.frame(sister1=binomial, n1=n, lat1=lat, lon1=lon, cv1=cv_logmass)))
sister_join <- left_join(sister_join, with(bird_summ, data.frame(sister2=binomial, n2=n, lat2=lat, lon2=lon, cv2=cv_logmass)))
sister_join <- sister_join[complete.cases(sister_join), ]

#sister_join <- mutate(sister_join, lat1abs = abs(lat1), lat2abs = abs(lat2))

# Sort them so that lat1 is always the lower latitude.
for (i in 1:nrow(sister_join)) {
  if (sister_join$lat1[i] > sister_join$lat2[i]) {
    tmp <- sister_join[i,]
    sister_join[i, ] <- tmp[c(2,1,3,8,9,10,11,4,5,6,7)]
  }
}

sister_join <- mutate(sister_join, dlat = lat2 - lat1)

# Get rid of orioles. They are an outlier when it comes to climate.
sister_join <- filter(sister_join, !grepl('Oriolus',sister1))

# With all quality controls in place, we now have 95 sister pairs.
sister_troptemp <- filter(sister_join, lat1 < 23.5 & lat2 > 23.5) %>%
  mutate(dcv = cv2 - cv1)

with(sister_troptemp, t.test(cv2, cv1, paired = TRUE, alternative = 'greater')) 

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
quantile(sapply(ttests, function(x) as.numeric(x$estimate)), prob=ps) 
quantile(sapply(ttests, function(x) as.numeric(x$stat)), prob=ps) 
quantile(sapply(ttests, function(x) as.numeric(x$p.val)), prob=ps) 
sum( sapply(ttests, function(x) as.numeric(x$p.val)) < 0.05)/999 