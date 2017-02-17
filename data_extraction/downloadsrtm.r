# Get elevations from SRTM for all the bird points (lots of SRTM tiles will be downloaded!!!)

setwd('/mnt/research/plz-lab/SRTM/')

# Load points
bird <- read.csv('~/verts/bird_coords_mass.csv', stringsAsFactors = FALSE)
bird <- subset(bird, !is.na(lat) & !is.na(lon))

# Find tiles that need to be downloaded
bird <- transform(bird, latr = round(lat), lonr = round(lon))
tilecoords <- unique(cbind(bird$latr, bird$lonr))

library(raster)

for (i in 1:nrow(tilecoords)) {
	try(getData('SRTM', lon = tilecoords[i,2], lat = tilecoords[i,1]), TRUE)
	print(i)
}
