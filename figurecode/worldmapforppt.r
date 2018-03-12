worldMap <- borders('world', fill='beige', color='black')
p_map <- ggplot() + worldMap + 
  geom_segment(aes(x=lon_trop, y=lat_trop, xend=lon_nontrop, yend=lat_nontrop, linetype = d>0), data=mapdat_sister, color = 'gray40') +
  geom_point(aes(x=lon_trop,y=lat_trop, fill=tropical_cv_logmass), data=mapdat_sister, pch=21) +
  geom_point(aes(x=lon_nontrop,y=lat_nontrop, fill=nontropical_cv_logmass), data=mapdat_sister, pch=21) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(breaks = c(-23.5, 23.5), labels = c('23.5° S', '23.5° N'), expand=c(0,0)) +
  fillScale + linetypeScale +
  theme(panel.background = element_rect(fill='skyblue'),
        panel.grid.major.y = element_line(color='black', size=1.1),
        panel.grid.major.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color='white'),
        legend.direction = 'horizontal',
        legend.justification = c(0,0),
        legend.text = element_text(size=10, color = 'white'),
        legend.title = element_text(color = 'white'),
        legend.position = 'bottom',
        plot.background = element_rect(fill='black', color = 'black')) + 
  coord_map('albers',0,0, ylim = c(-60,60)) +
  panel_border(colour='black') +
  guides(linetype = guide_legend(override.aes=list(color = 'white')))
ggsave('C:/Users/Q/Dropbox/presentations/sesync2018/birdmap.png', p_map, height = 5, width = 8, dpi = 400)
