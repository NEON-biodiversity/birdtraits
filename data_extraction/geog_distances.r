# Latitude as continuous predictor (rather than t-test)
# 8 June 2017

vnsis <- vnbird %>% filter(issister)
#write.csv(with(vnsis, data.frame(binomial, lat=decimallatitude, lon=decimallongitude, massing=massing)), file = file.path(fp, 'bird_coords_mass.csv'), row.names = FALSE)

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


##################################################################

# Analysis 1: Use the temperate and tropical pairs.
# NOTE: It's possible we may have to remove Butorides because the two species were lumped at one point.

sister_troptemp <- filter(sister_join, lat1 < 23.5 & lat2 > 23.5) %>%
  mutate(dcv = cv2 - cv1)

ggplot(sister_troptemp, aes(x=dlat, y=dcv)) + geom_point()
lmlat <- lm(dcv ~ dlat, data = sister_troptemp)
summary(lmlat) # Not significant.

sister_join <- mutate(sister_join, dcv = cv2 - cv1)
ggplot(sister_join, aes(x=dlat, y=dcv)) + geom_point()
lmlat <- lm(dcv ~ dlat, data = sister_join)
summary(lmlat) # Not significant.

# Use "great circle" distance between the points, in kilometers, just to see.

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

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(pmin(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

sister_troptemp <- sister_troptemp %>%
  mutate(dist_slc = gcd.slc(lon1, lat1, lon2, lat2),
         dist_hf = gcd.hf(lon1, lat1, lon2, lat2)) # The two distances are practically identical.
