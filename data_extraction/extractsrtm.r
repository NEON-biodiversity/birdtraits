# Get elevations from SRTM for all the bird points (lots of SRTM tiles will be downloaded!!!)
# Now that the data are downloaded we can actually get the elevations out.
# This might take a really long time.

setwd('/mnt/research/plz-lab/SRTM/')

# Load points
bird <- read.csv('~/verts/birdcoords_23Jun.csv', stringsAsFactors = FALSE)
bird <- subset(bird, !is.na(decimallatitude) & !is.na(decimallongitude))

# Get coordinates for all points that need elevations
#bird <- transform(bird, latr = round(lat), lonr = round(lon))
birdcoords <- unique(cbind(bird$decimallatitude, bird$decimallongitude))

library(raster)

# Another way of doing it is loading each tile as a raster, then getting the elevation of each coordinate in there.

rasternames <- dir(pattern = '*.tif')
elevs <- numeric(nrow(birdcoords))

for (i in rasternames) {
	ri <- raster(i)
	exti <- extent(ri)
	isinbox <- birdcoords[,2] >= exti@xmin & birdcoords[,2] <= exti@xmax & birdcoords[,1] >= exti@ymin & birdcoords[,1] <= exti@ymax
	if (sum(isinbox) > 0) {
		elevsi <- extract(x = ri, y = cbind(x = birdcoords[isinbox,2], y = birdcoords[isinbox,1]))
		elevsi[is.na(elevsi)] <- 0 # Distinguish true missing values from sea level values.
		elevs[isinbox] <- elevsi
	}
	print(i)
}

elevs <- data.frame(lat = birdcoords[,1], lon = birdcoords[,2], elev = elevs)
write.csv(elevs, file = '~/verts/bird_elev.csv', row.names = FALSE)

