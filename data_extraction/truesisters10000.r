
library(ape)
allbirds <- read.tree('~/verts/trees/AllBirdsEricson1.tre')

for (i in 2:10) allbirds <- c(allbirds, read.tree(paste0('~/verts/trees/BirdzillaEricson',i,'.tre')))

# Use all the trees at this point.
# Select 50 trees at random to use.
#set.seed(6880535)
#idx <- sample(1:1000, 50)
#allbirds <- allbirds[idx]

# Calculate a cophenetic distance matrix.
ericsondist <- lapply(allbirds, cophenetic.phylo)



# Must sort the distance matrices so that they are all the same order.
roworder <- dimnames(ericsondist[[1]])[[1]]
ericsondist <- lapply(ericsondist, function(x) x[roworder, roworder])

# Convert list of matrices to an array
library(abind)
ericsondist <- abind(ericsondist, along = 3)

# Mean and upper and lower quantile of branch lengths across all trees
ericsondistmean <- apply(ericsondist, 1:2, mean)
ericsondistlower <- apply(ericsondist, 1:2, quantile, probs = 0.025)
ericsondistupper <- apply(ericsondist, 1:2, quantile, probs = 0.975)
#ericsondistmean <- apply(simplify2array(ericsondist), 1:2, mean) #This overflows! :-O

# A true sister pair is defined as taxa A and B such that the nearest neighbor of taxon A is taxon B, and the nearest neighbor of taxon B is taxon A and only taxon A.

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

sisters_df <- identify_sisters(ericsondistmean)
sisters_df_lower <- identify_sisters(ericsondistlower)
sisters_df_upper <- identify_sisters(ericsondistupper)

write.csv(sisters_df, file = '~/verts/truesisters12jun.csv', row.names = FALSE)
write.csv(sisters_df_lower, file = '~/verts/truesisters12jun_lower.csv', row.names = FALSE)
write.csv(sisters_df_upper, file = '~/verts/truesisters12jun_upper.csv', row.names = FALSE)
