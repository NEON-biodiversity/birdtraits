# Figures for multipanel plot (main figure of current biology manuscript)

# a. map
# b. interaction plot
# c, d. covariate plots

library(cowplot)

### Interaction plot OR histogram
ylabel <- expression(paste('CV of log'[10], ' body mass'))		

plotdat <- data.frame(realm = rep(c('tropical','nontropical'), each=nrow(sister_troptemp)),
                      taxon = c(sister_troptemp$sister1, sister_troptemp$sister2),
                      pairid = 1:nrow(sister_troptemp),
                      cv_logmass = c(sister_troptemp$cv1, sister_troptemp$cv2))

p_int <- ggplot(plotdat, aes(x = realm, y = cv_logmass)) +
  geom_line(aes(group = pairid), color = 'gray75') + 
  geom_point(aes(group = pairid)) +
  stat_summary(aes(group = 1), geom = 'line', color = 'red', fun.y = 'mean', size = 1.5) +
  stat_summary(geom = 'pointrange', color = 'red', fun.data = 'mean_se') +
  scale_x_discrete(expand=c(0.1,0.1), name='Realm') +
  labs(y = ylabel) +
  panel_border(colour='black')

phist <- ggplot(sister_troptemp, aes(x = cv1 - cv2)) +
  geom_histogram(bins = 15) +
  geom_vline(xintercept = -0.020, color = 'red', lwd = 2) +
  geom_vline(xintercept = -0.0092, color = 'red', lwd = 0.5) +
  geom_vline(xintercept = 0, color = 'blue', lty = 3) +
  scale_y_continuous(limits = c(0,53), expand = c(0,0)) +
  labs(x = 'tropical CV - nontropical CV', y = 'count of species pairs') +
  panel_border(colour='black')

### Covariate plots
hl <- geom_hline(yintercept = 0, color = 'indianred3')
vl <- geom_vline(xintercept = 0, color = 'indianred3')


pcov1 <- ggplot(lmdat2, aes(x = seasonal_var_temp1 - seasonal_var_temp2, y = d)) +
  hl + vl +
  geom_point() + stat_smooth(method='lm', se=F) +
  panel_border(colour='black') + 
  #ggtitle('Reduced variability in tropical species is correlated\nwith reduced SEASONAL temperature variability') +
  labs(x = expression(atop(paste('seasonal ',Delta * CV[temperature]), 'temperate more variable <--> tropical more variable')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass]))) 

pcov2 <- ggplot(lmdat2, aes(x = meanlogbodymass, y = d)) +
  hl +
  geom_point() + stat_smooth(method='lm', se=F) +
  panel_border(colour='black') + 
  #ggtitle('Reduced variability in tropical species is correlated\nwith larger body sizes') +
  labs(x = expression(paste(log[10],' body mass')), y = expression(atop('temperate more variable <--> tropical more variable', Delta * CV[bodymass]))) 


###
# Arrange plots with cowplot

smalllab <- theme(axis.title.x=element_text(size=10),
                  axis.title.y=element_text(size=10))

bottom_row <- plot_grid(p_int, pcov1 + smalllab, pcov2 + smalllab, labels = c('b', 'c', 'd'), align = 'h', rel_widths = c(1, 1.4, 1.4), ncol=3)
full_plot <- plot_grid(p_map + panel_border(colour='black'), bottom_row, labels = c('a', ''), ncol = 1, rel_heights = c(1.2, 1))
ggsave(file.path(fp, 'vertnet_results/multipanelplot1.png'), full_plot, height=8, width=8, dpi=400)
