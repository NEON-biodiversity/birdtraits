# Extract climate for sister taxa records.

library(raster)
bird <- read.csv('~/verts/birdcoords_23Jun.csv', stringsAsFactors = FALSE)
fp <- '/mnt/research/plz-lab/NEON/external_data/raw_external_data/bioclim'
bird_coord <- data.frame(x = bird$decimallongitude, y = bird$decimallatitude)

bird_clim <- list()

for (i in dir(fp)) {
	print(i)
	ri <- raster(file.path(fp, i))
	bird_clim[[length(bird_clim) + 1]] <- extract(ri, bird_coord, method = 'simple')
}

# Format into data frame.

col_names <- sapply(strsplit(dir(fp), '_'), '[', 2)
bird_clim <- do.call('cbind', bird_clim)
bird_clim <- as.data.frame(bird_clim)
names(bird_clim) <- col_names
bird_clim <- bird_clim[, c(paste0('bio',1:19), 'prec','temp')]
names(bird_clim)[20:21] <- c('prec_interannual', 'temp_interannual')

write.csv(bird_clim, file = '~/verts/bird_clim_23Jun.csv', row.names = FALSE)
