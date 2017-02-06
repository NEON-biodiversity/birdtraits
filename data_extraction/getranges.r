# Get geographic ranges for birds from the Clements checklist.

checklist <- read.csv(file.path(fp, 'eBird-Clements-Checklist-v2016-10-August-2016.csv'), stringsAsFactors = FALSE)

sister_names2 <- with(sister_troptemp, c(sister1, sister2))
checklist$Scientific.name <- gsub(' ', '_', checklist$Scientific.name)

sismatch <- sister_names2 %in% checklist$Scientific.name
wrong_names <- sister_names2[!sismatch]

corr_names <- c('Ardenna_pacifica', 'Melozone_albicollis', 'Peucaea_sumichrasti', 'Setophaga_pitiayumi', 'Euplectes_psammacromius', 'Ardenna_bulleri', 'Melozone_fusca', 'Peucaea_carpalis', 'Setophaga_americana')

for (i in 1:length(wrong_names)) sister_troptemp[sister_troptemp == wrong_names[i]] <- corr_names[i]

sister_names2 <- with(sister_troptemp, c(sister1, sister2))
sismatch <- sister_names2 %in% checklist$Scientific.name

left_join(data.frame(Scientific.name = sister_names2), checklist %>% select(Scientific.name, English.name, Range))

##################################################################

# Range sizes from the Birds of the World dataset. Calculated them in QGIS, now load into R
botw <- read.csv(file.path(fp,'botw_table.csv'), stringsAsFactors = FALSE)

# Should be able to get the migrant status out of this too.
botw <- mutate(botw, taxon = gsub(' ', '_', SCINAME))

botwmatch <- sister_pair_long$taxon %in% botw$taxon
#write.csv(sister_pair_long$taxon[!botwmatch], file = file.path(fp, 'botwnomatch.csv'), row.names = F)

# Reload the manually corrected name list
botwcorrect <- read.csv(file.path(fp, 'botwnomatch.csv'), stringsAsFactors = F)
for (i in 1:nrow(botwcorrect)) {
  idx <- botw$taxon == botwcorrect$correct[i]
  botw$taxon[idx] <- botwcorrect$x[i]
}

# Get the ranges.
botw_sisters <- subset(botw, taxon %in% sister_pair_long$taxon) # Ranges for 201 species are here.

# Get migrant status and total range size of resident+breeding range
# Codes for presence are 1 extant, 2 prob extant, 3 poss extant, 4 poss extinct, 5 extinct, 6 uncertain
# Codes for season are 1 resident, 2 breeding, 3 nonbreeding or winter, 4 passage, 5 uncertain

rangedat <- function(x) {
  migrant_status <- 'none'
  if (any(x$SEASONAL %in% 2:3)) migrant_status <- 'partial'
  if (!any(x$SEASONAL == 1)) migrant_status <- 'obligate'
  data.frame(migrant_status = migrant_status, range_size = sum(x$Area_km2[x$SEASONAL %in% 1:2]))
}

botw_ranges <- botw_sisters %>%
  filter(PRESENCE == 1) %>% # Only where species is known present.
  group_by(taxon) %>%
  do(rangedat(.))

write.csv(botw_ranges, file = file.path(fp, 'botw_ranges.csv'), row.names = FALSE)
