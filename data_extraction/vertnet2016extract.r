# Extract trait and georeference data from Vertnet
# Author: QDR
# Project: NEON ITV
# Created: 09 Jan 2017

# Modified 18 Jan: Extract other taxa.

fp <- '/mnt/research/plz-lab/NEON/external_data/raw_external_data/vertnet'

#vnet <- read.csv(file.path(fp, 'vertnet_20151028_for_traits_extraction.csv'), stringsAsFactors = FALSE) # 11gb so takes time to load
vnmam <- read.csv(file.path(fp, 'vertnet_latest_mammals_09Jan2017.csv'), stringsAsFactors = FALSE) #8gb

library(dplyr)

# Get rid of the rows that have neither length nor mass
# Get rid of some of the columns so that I can reload this dataframe later without consuming so much ram.
vnmam <- vnmam %>% 
	filter(haslength | hasmass) %>%
	select(scientificname, vernacularname, lengthinmm, massing, lifestage, sex, reproductivecondition, class, order, family, genus, subgenus, specificepithet, infraspecificepithet, organismid, organismname, wascaptive, wasinvasive, isfossil, day, month, year, eventdate, eventtime, eventid, country, decimallatitude, decimallongitude, coordinateprecision, coordinateuncertaintyinmeters, geodeticdatum, minimumelevationinmeters, maximumelevationinmeters, locality, locationid, occurrenceid, occurrencestatus, individualcount, samplingeffort, typestatus, datasetid)

write.csv(vnmam, file = file.path(fp, 'vertnet_mammals_reduced.csv'), row.names = FALSE)	
	
# In future, probably get rid of captives, invasives, and fossils.

# Extract the other taxa: birds, amphibians, fishes, reptiles.

########## birds

vnbird <- read.csv(file.path(fp, 'vertnet_latest_birds.csv'), stringsAsFactors = FALSE) #11gb

# Get rid of the rows that have neither length nor mass
# Get rid of some of the columns so that I can reload this dataframe later without consuming so much ram.
vnbird <- vnbird %>% 
	filter(haslength | hasmass) %>%
	select(scientificname, vernacularname, lengthinmm, massing, lifestage, sex, reproductivecondition, class, order, family, genus, subgenus, specificepithet, infraspecificepithet, organismid, organismname, wascaptive, wasinvasive, isfossil, day, month, year, eventdate, eventtime, eventid, country, decimallatitude, decimallongitude, coordinateprecision, coordinateuncertaintyinmeters, geodeticdatum, minimumelevationinmeters, maximumelevationinmeters, locality, locationid, occurrenceid, occurrencestatus, individualcount, samplingeffort, typestatus, datasetid)

write.csv(vnbird, file = file.path(fp, 'vertnet_birds_reduced.csv'), row.names = FALSE)	

########### amphibians

vnamph <- read.csv(file.path(fp, 'vertnet_latest_amphibians.csv'), stringsAsFactors = FALSE)

# Get rid of the rows that have neither length nor mass
# Get rid of some of the columns so that I can reload this dataframe later without consuming so much ram.
vnamph <- vnamph %>% 
	filter(haslength | hasmass) %>%
	select(scientificname, vernacularname, lengthinmm, massing, lifestage, sex, reproductivecondition, class, order, family, genus, subgenus, specificepithet, infraspecificepithet, organismid, organismname, wascaptive, wasinvasive, isfossil, day, month, year, eventdate, eventtime, eventid, country, decimallatitude, decimallongitude, coordinateprecision, coordinateuncertaintyinmeters, geodeticdatum, minimumelevationinmeters, maximumelevationinmeters, locality, locationid, occurrenceid, occurrencestatus, individualcount, samplingeffort, typestatus, datasetid)

write.csv(vnamph, file = file.path(fp, 'vertnet_amphibians_reduced.csv'), row.names = FALSE)	

########### reptiles

vnrept <- read.csv(file.path(fp, 'vertnet_latest_reptiles.csv'), stringsAsFactors = FALSE)

# Get rid of the rows that have neither length nor mass
# Get rid of some of the columns so that I can reload this dataframe later without consuming so much ram.
vnrept <- vnrept %>% 
	filter(haslength | hasmass) %>%
	select(scientificname, vernacularname, lengthinmm, massing, lifestage, sex, reproductivecondition, class, order, family, genus, subgenus, specificepithet, infraspecificepithet, organismid, organismname, wascaptive, wasinvasive, isfossil, day, month, year, eventdate, eventtime, eventid, country, decimallatitude, decimallongitude, coordinateprecision, coordinateuncertaintyinmeters, geodeticdatum, minimumelevationinmeters, maximumelevationinmeters, locality, locationid, occurrenceid, occurrencestatus, individualcount, samplingeffort, typestatus, datasetid)

write.csv(vnrept, file = file.path(fp, 'vertnet_reptiles_reduced.csv'), row.names = FALSE)	


########### fishes

vnfish <- read.csv(file.path(fp, 'vertnet_latest_fishes.csv'), stringsAsFactors = FALSE)

# Get rid of the rows that have neither length nor mass
# Get rid of some of the columns so that I can reload this dataframe later without consuming so much ram.
vnfish <- vnfish %>% 
	filter(haslength | hasmass) %>%
	select(scientificname, vernacularname, lengthinmm, massing, lifestage, sex, reproductivecondition, class, order, family, genus, subgenus, specificepithet, infraspecificepithet, organismid, organismname, wascaptive, wasinvasive, isfossil, day, month, year, eventdate, eventtime, eventid, country, decimallatitude, decimallongitude, coordinateprecision, coordinateuncertaintyinmeters, geodeticdatum, minimumelevationinmeters, maximumelevationinmeters, locality, locationid, occurrenceid, occurrencestatus, individualcount, samplingeffort, typestatus, datasetid, waterbody, minimumdepthinmeters, maximumdepthinmeters)

write.csv(vnfish, file = file.path(fp, 'vertnet_fishes_reduced.csv'), row.names = FALSE)
