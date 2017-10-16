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

####
# 10 Oct.
# Use centroids calculated on hpcc to check whether classification based on centroids of year-round or breeding polygons
# matches classification based on median latitude of specimens.

library(dplyr)

fprev <- 'C:/Users/Q/Dropbox/projects/verts/manuscript/revision'
botw_cents <- read.csv(file.path(fprev, 'polygon_centroids.csv'), stringsAsFactors = FALSE)

# Metadata location http://datazone.birdlife.org/species/spcdistPOS
# Use only resident and breeding polygons, native, and presence status 1-3
botw_filtered <- botw_cents %>%
  filter(ORIGIN %in% 1, SEASONAL %in% 1:2, PRESENCE %in% 1:3) %>%
  mutate(realm = 'tropical')

botw_filtered$realm[abs(botw_filtered$X2) > 23.5] <- 'nontropical'


#####################

####################################################################################################
# Load data that has been flagged for bad values and outliers that are due to typos, and clean it up.

vnsis_flag <- read.csv('C:/Users/Q/Dropbox/projects/verts/bird_records_flagged.csv', stringsAsFactors = FALSE)
table(vnsis_flag$flag)

# Remove invalid, not_native, and not_wild. If typo or absent, use the corrected numbers I put in. If listed as opposite sign, change the sign and replace.

vnsis_flag <- filter(vnsis_flag, !flag %in% c('invalid','not_native','not_wild'))
latflag1 <- vnsis_flag$flag %in% c('lat_typo','latlong_typo','latlong_absent')
latflag2 <- vnsis_flag$flag %in% c('lat_opp','latlong_opp')
vnsis_flag$decimallatitude[latflag1] <- as.numeric(vnsis_flag$correct_lat[latflag1])
vnsis_flag$decimallatitude[latflag2] <- -vnsis_flag$decimallatitude[latflag2]

lonflag1 <- vnsis_flag$flag %in% c('long_typo','latlong_typo','latlong_absent')
lonflag2 <- vnsis_flag$flag %in% c('long_opp','latlong_opp')
vnsis_flag$decimallongitude[lonflag1] <- as.numeric(vnsis_flag$correct_lon[lonflag1])
vnsis_flag$decimallongitude[lonflag2] <- -vnsis_flag$decimallongitude[lonflag2]

# 22 June: flag those records that are greater than or less than 10x the median value for the species.
# 10 times is a lot for things to be varying and probably indicates a mistake.

vnfactor <- vnsis_flag %>% group_by(binomial) %>%
  summarize(median_mass = median(massing, na.rm = TRUE))
vnsis_flag2 <- vnsis_flag %>% left_join(vnfactor) %>% mutate(outlierflag = massing > median_mass * 10 | massing < median_mass / 10)
table(vnsis_flag2$outlierflag) # 1833 outliers flagged if 5x, 1178 flagged if 10x

# Try to match botw with vertnet.
botw_filtered <- mutate(botw_filtered, binomial = gsub('\\ ', '_', SCINAME))

vnsis_spnames <- unique(vnsis_flag2$binomial)

badnames <- vnsis_spnames[!vnsis_spnames %in% botw_filtered$binomial]
# Must correct about 300 species names to make these match. Argh.
#write.csv(data.frame(vertnet_name = badnames), file=file.path(fprev, 'correctednames.csv'), row.names = FALSE)

#z <- function(x) grep(x, botw_filtered$binomial, value = TRUE)

# Load manually corrected names.
goodnames <- read.csv(file.path(fprev, 'correctednames_botw_vertnet.csv'), stringsAsFactors = FALSE)

# Change the BOTW names to Vertnet names.
goodnames <- mutate(goodnames,
                    vertnet_species = sapply(strsplit(vertnet_name, '_'), '[', 2))
goodnames$vertnet_species[!goodnames$botw_species==''] <- goodnames$botw_species[!goodnames$botw_species=='']
goodnames <- mutate(goodnames,
                    botw_name = paste(botw_genus, vertnet_species, sep = '_'))

for (i in 1:nrow(goodnames)) {
  idx <- which(botw_filtered$binomial == goodnames$botw_name[i])
  botw_filtered$binomial[idx] <- goodnames$vertnet_name[i] 
}

vnsis_spnames[!vnsis_spnames %in% botw_filtered$binomial]

# Correct the names in vertnet if needed ***TO*** the BOTW name.
# This only refers to lumped species
vnsis_flag2$binomial[vnsis_flag2$binomial == 'Aphelocoma_insularis'] <- 'Aphelocoma_californica'
vnsis_flag2$binomial[vnsis_flag2$binomial == 'Carduelis_hornemanni'] <- 'Carduelis_flammea'
vnsis_flag2$binomial[vnsis_flag2$binomial == 'Calandrella_cheleensis'] <- 'Calandrella_rufescens'
vnsis_flag2$binomial[vnsis_flag2$binomial == 'Serinus_whytii'] <- 'Serinus_striolatus'
vnsis_flag2$binomial[vnsis_flag2$binomial == 'Butorides_striata'] <- 'Butorides_virescens'

vnsis_spnames <- unique(vnsis_flag2$binomial)
vnsis_spnames[!vnsis_spnames %in% botw_filtered$binomial]

# they all match!!!

# Join.
vnsis_flag2 <- left_join(vnsis_flag2,
                         botw_filtered %>% select(binomial, realm))

# Find the sister pairs w/at least 10 specimens apiece.

# First redo the sisters with a CONSENSUS tree!!!
# (Done on cluster with 100 trees sampled randomly)

library(ape)
birdtree <- read.tree(file.path(fp, 'ericson_cons.tre'))
ericsondist <- cophenetic.phylo(birdtree)

# Function to identify sisters from a distance matrix, then run on the mean, upperci, and lowerci.

identify_sisters <- function(mat) {
  
  sisters <- character(nrow(mat))
  # Modified 08 June: add the phylogenetic distance of sister pairs as a column for later analysis.
  dist_to_sister <- numeric(nrow(mat))
  
  for (i in 1:nrow(mat)) {
    x <- mat[i, -i]
    name_i <- dimnames(mat)[[1]][i]
    dist_to_sister[i] <- min(x)
    sisternames <- names(x)[x == dist_to_sister[i]] # Names of potential sisters.
    # Check which of these sisters has taxon i as a sister.
    for (j in 1:length(sisternames)) {
      j_index <- which(dimnames(mat)[[1]] == sisternames[j])
      y <- mat[j_index, -j_index]
      reciprocaldist_to_sister <- min(y)
      reciprocalsisternames <- names(y)[y == reciprocaldist_to_sister]
      if (name_i %in% reciprocalsisternames) sisters[[i]] <- sisternames[j]
    }
    
  }
  
  sisters_df <- data.frame(sister1 = dimnames(mat)[[1]], sister2 = sisters, dist = dist_to_sister)
  sisters_df <- subset(sisters_df, sister2 != '')
  sisters_df <- transform(sisters_df, sister1 = as.character(sister1), sister2 = as.character(sister2))
  
  # Get rid of duplicates.
  is_dup <- rep(FALSE, nrow(sisters_df))
  
  for (i in 1:nrow(sisters_df)) {
    sisterrows <- which(sisters_df$sister2 == sisters_df$sister1[i] & sisters_df$sister1 == sisters_df$sister2[i])
    if (length(sisterrows) > 0)
      if (sisterrows[1] < i)
        is_dup[i] <- TRUE
  }
  
  return(sisters_df[!is_dup,])
  
}

sisters_df <- identify_sisters(ericsondist)


#####################################################################

# In the interim, do this with the old sisters.
fp <- 'C:/Users/Q/Dropbox/projects/verts'

sister_summary <- read.csv(file.path(fp, 'sister_summary.csv'), stringsAsFactors = FALSE)

# Get only those that are supported in a good percentage of the trees.
sister_goodpairs <- subset(sister_summary, proportion >= 0.1)
length(unique(sister_goodpairs$sister2))

# Remove duplicates
sister_goodpairs$duplicate <- FALSE
for (i in 2:nrow(sister_goodpairs)) {
  previous <- c(sister_goodpairs$sister1[1:(i-1)], sister_goodpairs$sister2[1:(i-1)])
  if (sister_goodpairs$sister1[i] %in% previous | sister_goodpairs$sister2[i] %in% previous) sister_goodpairs$duplicate[i] <- TRUE
}

sister_goodpairs <- subset(sister_goodpairs, !duplicate)

bird_summ <- vnsis_flag2 %>% 
  filter(!outlierflag) %>% 
  group_by(binomial, realm) %>% 
  summarize(n = n(), 
            lat = median(abs(decimallatitude), na.rm=T),
            lon = median(decimallongitude, na.rm=T),
            cv_logmass = sd(log10(massing))/mean(log10(massing)))

sister_join <- left_join(sister_goodpairs, with(bird_summ, data.frame(sister1=binomial, n1=n, lat1=lat, lon1=lon, cv1=cv_logmass, realm1=realm)))
sister_join <- left_join(sister_join, with(bird_summ, data.frame(sister2=binomial, n2=n, lat2=lat, lon2=lon, cv2=cv_logmass, realm2=realm)))
sister_join <- sister_join[complete.cases(sister_join), ]

# Sort them so that lat1 is always the lower latitude.
for (i in 1:nrow(sister_join)) {
  if (sister_join$lat1[i] > sister_join$lat2[i]) {
    tmp <- sister_join[i,]
    sister_join[i, c('sister1','n1','lat1','lon1','cv1')] <- tmp[c('sister2','n2','lat2','lon2','cv2')]
    sister_join[i, c('sister2','n2','lat2','lon2','cv2')] <- tmp[c('sister1','n1','lat1','lon1','cv1')]
  }
}

sister_join <- mutate(sister_join, dlat = lat2 - lat1, dcv = cv2 - cv1)

####################################################################################################
# Use only tropical-nontropical pairs and run t-test.

# Get rid of orioles. They are an outlier when it comes to climate.
sister_join <- filter(sister_join, !grepl('Oriolus',sister1))

# With all quality controls in place, we now have 95 sister pairs.
sister_troptemp <- filter(sister_join, realm1 == 'tropical' & realm2 == 'nontropical')

with(sister_troptemp, t.test(cv2, cv1, paired = TRUE, alternative = 'greater')) 
