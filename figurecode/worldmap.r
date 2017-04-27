# Map

fp <- 'C:/Users/Q/Dropbox/projects/verts'
load(file.path(fp, 'sistercovariates.r'))
mapdat <- read.csv(file.path(fp, 'mapdata.csv'))

library(dplyr)

sister_pair_long <- left_join(sister_pair_long, mapdat %>% rename(taxon=binomial) %>% select(taxon,lon))

mapdat$lat_temp <- mapdat$lat_trop <- mapdat$lat
mapdat$lon_temp <- mapdat$lon_trop <- mapdat$lon
mapdat$lat_temp[abs(mapdat$lat_temp) < 23.5] <- NA
mapdat$lon_temp[abs(mapdat$lat_temp) < 23.5] <- NA
mapdat$lat_trop[abs(mapdat$lat_trop) > 23.5] <- NA
mapdat$lon_trop[abs(mapdat$lat_trop) > 23.5] <- NA

sister_alldat <- left_join(sister_alldat, 
                           mapdat %>% rename(sister1=binomial) %>% select(sister1, lat_trop,lon_trop))
sister_alldat <- left_join(sister_alldat, 
                           mapdat %>% rename(sister2=binomial) %>% select(sister2, lat_temp,lon_temp))

# Create map with colored points for the different taxa. Color should correspond to the variability. Lines connecting the pairs will be red if tropical has higher var, blue if temperate

library(ggplot2)

heatramp <- colorRampPalette(RColorBrewer::brewer.pal(name='YlOrRd', n=9),bias=2,space="rgb")(50)
fillScale <- scale_fill_gradientn(colours = heatramp, name = 'Body mass variability')
colorScale <- scale_color_manual(values = c('indianred','forestgreen')) # Line colors

worldMap <- borders('world', fill='gray75', color='black')
p_map <- ggplot() + worldMap + 
  geom_segment(aes(x=lon_temp, y=lat_temp, xend=lon_trop, yend=lat_trop, color = d<0), data=sister_alldat) +
  geom_point(aes(x=lon_trop,y=lat_trop, fill=cv1), data=sister_alldat, pch=21) +
  geom_point(aes(x=lon_temp,y=lat_temp, fill=cv2), data=sister_alldat, pch=21) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(breaks = c(-23.5, 23.5), labels = c('23.5° S', '23.5° N'), expand=c(0,0)) +
  fillScale + colorScale +
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

ggsave(file.path(fp,'vertnet_results/worldmap.png'), height=4, width=8, dpi=400)
