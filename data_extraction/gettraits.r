# Look at covariates of result from paired t-test

sister_names2 <- with(sister_troptemp, c(sister1,sister2))

corr_names <- c('Ardenna_pacifica', 'Melozone_albicollis', 'Peucaea_sumichrasti', 'Setophaga_pitiayumi', 'Euplectes_psammacromius', 'Ardenna_bulleri', 'Melozone_fusca', 'Peucaea_carpalis', 'Setophaga_americana')

# Functional traits
foraging <- read.delim(file.path(fp, 'BirdFuncDat.txt'), stringsAsFactors = FALSE)
lifehist <- read.csv(file.path(fp, 'Amniote_Database_Aug_2015.csv'), stringsAsFactors = FALSE)

foraging$binomial <- gsub(' ', '_', foraging$Scientific)
sismatchfor <- sister_names2 %in% foraging$binomial
sister_names2[!sismatchfor] # The original (wrong) names match here.

foraging <- foraging[,c(41, 10:20, 24:31, 35,36)]
foraging_sisters <- left_join(data.frame(binomial = sister_names2), foraging)

lifehist$binomial <- paste(lifehist$genus, lifehist$species, sep = '_')
sismatchlh <- sister_names2 %in% lifehist$binomial
sister_names2[!sismatchlh]

corr_names_lh <- c('Melaenornis_microrhynchus', 'Gymnoris_dentata', 'Crithagra_sulphurata', 'Chionodacryon_speculiferum', 'Melozone_albicollis', 'Peucaea_sumichrasti', 'Setophaga_pitiayumi', 'Euplectes_psammacromius', 'Basilinna_leucotis', 'Niltava_davidi', 'Melaenornis_mariquensis', 'Gymnoris_xanthocollis', 'Crithagra_flaviventris', 'Melozone_fusca', 'Peucaea_carpalis', 'Setophaga_americana', 'Basilinna_xantusii')

sister_names_lh <- sister_names2
sister_names_lh[!sismatchlh] <- corr_names_lh

lifehist_sisters <- left_join(data.frame(binomial = sister_names_lh), lifehist[,c(37, 8:36)])

trait_sisters <- cbind(foraging_sisters, lifehist_sisters[,-1])

write.csv(trait_sisters, file = file.path(fp,'trait_sisters.csv'), row.names = FALSE)
