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
table(vnsis_flag$outlierflag)

write.csv(vnsis_flag, file = 'C:/Users/Q/Dropbox/projects/verts/bird_records_flagged.csv', row.names = FALSE)
