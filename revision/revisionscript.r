# Code to do suggested methods revisions for Biology Letters
# 6 Oct 2017

# To do:

# 1. Check whether any of the species are cosmopolitan (neither just temp or just trop) and if so, get rid.
# 2. Try out the phylogenetic correction using the basal node of each sister pair.
# 3. Model selection: check VIF, do forward not backward, and do AICc not AIC.

fprev <- 'C:/Users/Q/Dropbox/projects/verts/manuscript/revision'

# 1. Check cosmopolitan species.

# Modify the below code to only plot species from the 78 pairs

# Load vnsis and load the species names.

vn_toplot <- vnsis_flag2 %>% 
  filter(!outlierflag, binomial %in% with(sister_troptemp, c(sister1, sister2)))

library(cowplot)
spnames <- unique(vn_toplot$binomial)

pb <- txtProgressBar(0, length(spnames), style = 3)
count <- 0
pdf(file.path(fprev,'allmaps_vertnetcoords.pdf'))
for (i in spnames) {
  count <- count + 1
  setTxtProgressBar(pb, count)
  data_i <- vn_toplot %>% filter(binomial == i) %>% select(decimallongitude, decimallatitude) %>% filter(complete.cases(.))
  if (nrow(data_i) > 0) {
    xlim_i <- range(data_i$decimallongitude) + c(-5,5)
    ylim_i <- range(data_i$decimallatitude) + c(-5,5)
    p <- ggplot(subset(vn_toplot, binomial == i), aes(x = decimallongitude, y = decimallatitude)) +
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

