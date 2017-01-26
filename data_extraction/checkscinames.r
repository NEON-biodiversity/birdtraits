# Vertnet birds taxize

library(taxize)
scinames <- read.csv('C:/Users/Q/Dropbox/projects/verts/birdscinames.csv', stringsAsFactors = FALSE)

checknames <- gnr_resolve(names = scinames$n) # This takes a fair amount of time since there are ~5000 to resolve.
write.csv(checknames, file = 'C:/Users/Q/Dropbox/projects/verts/checkednames.csv', row.names = FALSE)

# Isolate problematic names
table(checknames$score) # Most check out.
subset(checknames, score < 0.988)[,1:3]

head(subset(checknames, score < 0.988)[,1:3],150)

###############################################################
# The bird names have been manually corrected elsewhere.
scinames <- read.csv('C:/Users/Q/Dropbox/projects/verts/birdscinames_correct.csv', stringsAsFactors = FALSE)

# Try to get all synonyms
allsyns <- synonyms(scinames$n, db = 'nbn')

# Compare the list of scientific names from vertnet to the master taxonomy used by Jetz et al to make the phylogeny.

jetz_tax <- read.csv('C:/Users/Q/Dropbox/projects/verts/Jetz2012_SI/2012-03-04206D-master_taxonomy.csv', stringsAsFactors = FALSE)

# Which ones match
match1 <- scinames$n %in% jetz_tax$Scientific # There are 128 in our list that aren't in Jetz's list :-(

scinames$n[!match1]

match2 <- jetz_tax$Scientific %in% scinames$n

# Write the names to a file so I can edit them.
write.csv(data.frame(n = scinames$n, match = match1), file = 'C:/Users/Q/Dropbox/projects/verts/phylomatches.csv', row.names=FALSE)

# Also try to pull the synonyms to see if any of them are in the phylogeny.
syntable <- read.csv('C:/Users/Q/Dropbox/projects/aquaxterra/Supplemental_Table_6_Aves_Synonyms.csv', stringsAsFactors = FALSE)
syns <- with(syntable, data.frame(name1 = paste(genus, species), name2 = paste(synonymous_genus, synonymous_species)))
syns <- unique(syns)


library(dplyr)
scinames <- scinames %>% rename(name1=n) %>% left_join(syns)
#scinames$name2[scinames$name1==scinames$name2] <- NA

match_either <- scinames$name1 %in% jetz_tax$Scientific | scinames$name2 %in% jetz_tax$Scientific
scinames$name1[!match_either]

# flatten synonym table
synflat <- scinames %>% group_by(name1) %>% do(name2 = unique(.$name2))
synall <- sapply(synflat$name2, function(x) c(as.character(x), rep(NA,12))[1:13])
synflatall <- cbind(as.character(synflat$name1), t(synall))

match_any <- apply(synflatall, 1, function(x) any(x %in% jetz_tax$Scientific)) # There are still 529 that don't match.

match_df <- data.frame(match=match_any, synflatall[,1:5])

#########################################

# 23 Jan: all mismatches should now be fixed. Put them all together.
scinames_matched <- read.csv('C:/Users/Q/Dropbox/projects/verts/phylomatches_correct.csv', stringsAsFactors = FALSE)
name_to_use <- scinames_matched$n
name_to_use[!scinames_matched$match] <- scinames_matched$jetz_name[!scinames_matched$match]

# Export list of names to use
write.csv(data.frame(n = name_to_use), file = 'C:/Users/Q/Dropbox/projects/verts/phylo_namelist.csv', row.names=FALSE)

#########################################

# Now, load the bird phylogeny and try to fit all the species into it. There will probably be some that do not match.

library(ape)

# Read 100 trees containing the 774 species from birdtree.org (Ericson trees, randomly chosen subset)

erictree <- read.nexus('C:/Users/Q/Dropbox/projects/verts/ericson100.tre')

# Find sister species pairs
# There are only 734 rows in distance matrix because some of the species got lumped.

edist <- cophenetic(erictree[[1]])
epairs <- list()
for (i in 1:nrow(edist)) {
  x <- edist[i, edist[i, ] > 0]
  sisters <- names(x)[which(x == min(x))]
  epairs[[i]] <- sisters
  names(epairs)[[i]] <- dimnames(edist)[[1]][i]
}


# erictree <- read.nexus('~/verts/BirdzillaEricson2.tre')
# 
# # Load one of ten randomly selected trees (random numbers generated on 20 Jan)
# treeids <- c(507, 281, 191, 3, 920, 124, 990, 44, 654, 22)
# t1 <- erictree[[treeids[1]]]
# 
# # Find which species on our species list are not given as tips of the tree.
# 
# tlabel1 <- t1$tip.label
# tlabel1 <- gsub('_', ' ', tlabel1)
# 
# phymatch <- scinames$n %in% tlabel1 
# phymatchidx <- rep(NA,length(tlabel1))
# 
# for (i in 1:length(tlabel1)) {
# 	phymatchidx[i] <- c(which(scinames$n == tlabel1[i]), NA)[1]
# }
# 
# spp_not_in_tree <- scinames$n[!phymatch]