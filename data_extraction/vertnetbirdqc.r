# QC on vnbird

library(dplyr)

# Get rid of captive and fossil
vnbird_qc1 <- filter(vnbird, !wascaptive, !isfossil) # Leave in invasive for now.

# Retain only records with mass
vnbird_qc2 <- filter(vnbird_qc1, !is.na(massing))

# Attempt to get rid of juvenile records
lifet <- table(vnbird_qc2$lifestage)

#grep('ad', names(lifet), ignore.case = TRUE, value = TRUE)
#grep('sub', names(lifet), ignore.case = TRUE, value = TRUE)
remove_names <- grep('juv|sub', names(lifet), ignore.case = TRUE, value = TRUE)

vnbird_qc3 <- filter(vnbird_qc2, !lifestage %in% remove_names)

lifet3 <- table(vnbird_qc3$lifestage)

remove_names3 <- grep('im|nat|larv|young|nestl|fledg|unoss|pullus|hatch|chick|down', names(lifet3), ignore.case = TRUE, value = TRUE)
vnbird_qc4 <- filter(vnbird_qc3, !lifestage %in% remove_names3)
