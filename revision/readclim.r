# Extract all the climate for the coords.
# Run this remotely.

fprefix <- '/mnt/research/plz-lab/NEON/external_data/raw_external_data/bioclim/CHELSA_bio'
fsuffix <- '_1979-2013_V1_1.tif'

# Load coords
vncoords <- read.csv('~/verts/vertnet_coords.csv', stringsAsFactors = FALSE)
vncoords_good <- vncoords[!is.na(vncoords$decimallatitude), c('decimallongitude','decimallatitude')]

library(raster)

vnbioclim <- list()

for (i in 1:19) {
  x <- raster(paste0(fprefix, i, fsuffix))
  vnbioclim[[i]] <- extract(x, vncoords_good)
}

x <- raster('/mnt/research/plz-lab/NEON/external_data/raw_external_data/bioclim/CHELSA_prec_interannual_1979-2013_V1_1.tif')
prec_inter <- extract(x, vncoords_good)
x <- raster('/mnt/research/plz-lab/NEON/external_data/raw_external_data/bioclim/CHELSA_temp_interannual_1979-2013_V1_1.tif')
temp_inter <- extract(x, vncoords_good)


vnbioclim <- do.call('cbind', vnbioclim)
vnbioclim <- cbind(vnbioclim, prec_inter, temp_inter)
vnbioclim_mat <- matrix(NA, ncol=21, nrow=nrow(vncoords))
vnbioclim_mat[!is.na(vncoords$decimallatitude), ] <- vnbioclim

vnbioclim_mat <- as.data.frame(vnbioclim_mat)
names(vnbioclim_mat) <- c(paste0('bio',1:19), 'prec_interannual', 'temp_interannual')

write.csv(cbind(vncoords, vnbioclim_mat), file = '~/verts/vertnet_clim.csv', row.names = FALSE)
