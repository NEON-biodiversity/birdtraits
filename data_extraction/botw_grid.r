# Load the Birds of the World shape files that I converted with QGIS
# If necessary, combine them into a single file.
# Go across all grid cells in the world and calculate richness.
# Make sure the seasonality is correct.
# Also calculate richness by genus if I can.

n_tasks <- 50

fp <- '/mnt/research/aquaxterra/DATA/raw_data/bird_traits/BOTWshape'

library(maptools)
library(rgdal)

# b1 <- readOGR(dsn = fp, layer = 'BOTW00001to03000')
# b2 <- readOGR(dsn = fp, layer = 'BOTW03001to06000')
# b3 <- readOGR(dsn = fp, layer = 'BOTW06001to09000')
# b4 <- readOGR(dsn = fp, layer = 'BOTW09001to14000')
# b5 <- readOGR(dsn = fp, layer = 'BOTW14001to18111')

# Define projection and read the BOTW polygons in
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
bpoly1 <- readShapePoly(file.path(fp, 'BOTW00001to03000.shp'), proj4string = crswgs84)
bpoly2 <- readShapePoly(file.path(fp, 'BOTW03001to06000.shp'), proj4string = crswgs84)
bpoly3 <- readShapePoly(file.path(fp, 'BOTW06001to09000.shp'), proj4string = crswgs84)
bpoly4 <- readShapePoly(file.path(fp, 'BOTW09001to14000.shp'), proj4string = crswgs84)
bpoly5 <- readShapePoly(file.path(fp, 'BOTW14001to18111.shp'), proj4string = crswgs84)

bpoly <- rbind(bpoly1, bpoly2, bpoly3, bpoly4, bpoly5, makeUniqueIDs = TRUE)

polygondat <- bpoly@data[,c('SCINAME','PRESENCE','ORIGIN','SEASONAL')]
polygondat$genus <- sapply(strsplit(as.character(polygondat$SCINAME), ' '), '[', 1)
write.csv(polygondat, file = '~/verts/polygondat.csv', row.names = FALSE)

x <- as(bpoly, 'SpatialPolygons')

# Load points
bird <- read.csv('~/verts/bird_coords_mass.csv', stringsAsFactors = FALSE)
bird <- subset(bird, !is.na(lat) & !is.na(lon))

# Get slice of bird so that the command can be parallelized.
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
rowidx <- round(seq(0,nrow(bird),length.out=n_tasks + 1))
rowidxmin <- rowidx[task]+1
rowidxmax <- rowidx[task+1]

bird_task <- bird[rowidxmin:rowidxmax,]
bird_spatial <- SpatialPoints(cbind(x = bird_task$lon, y = bird_task$lat), proj4string = crswgs84)
objname <- paste0('overlap', task)

assign(objname, over(bird_spatial, x, returnList = TRUE))
save(list = objname, file = paste0('~/verts/mapoutput/', objname, '.R'))


#################################################################

# Compile the outputs that were run in parallel on the cluster.
# For each data point, get the richness and the richness of confamilial and congeneric individuals.

fp <- '~/verts/mapoutput'
lapply(as.list(file.path(fp, dir(fp, pattern = '23Jun'))), load, .GlobalEnv)
overlap_all <- do.call('c', lapply(paste0('overlap', 1:50), get))

bird <- read.csv('~/verts/birdcoords_23Jun.csv', stringsAsFactors = FALSE)
bird <- subset(bird, !is.na(decimallatitude) & !is.na(decimallongitude))

# Get genus for each of the data points
polygondat <- read.csv('~/verts/polygondat.csv', stringsAsFactors = FALSE)
bird$genus <- sapply(strsplit(bird$binomial, '_'), '[', 1)

# Get the richness and congeneric richness at each data point.
total_richness <- sapply(overlap_all, length)
congener_richness <- numeric(length(overlap_all))

bbin <- gsub('_', ' ', bird$binomial)

for (i in 1:length(overlap_all)) {
	spp_i <- polygondat[overlap_all[[i]], ]
	congener_richness[i] <- sum(bird$genus[i] == spp_i$genus & bbin[i] != spp_i$SCINAME)
}

richness_df <- data.frame(total_richness, congener_richness)
bird <- cbind(bird, richness_df)

write.csv(bird, '~/verts/bird_withrichness_23Jun.csv', row.names = FALSE)
