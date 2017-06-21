# Code to get lists of sisters from each of the 10000 matrices, in parallel.
# Split into 250 child jobs, then loop through 40 of them per job.

n_tasks <- 250 # 40 matrices per task.
n_mats <- 10000

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
idx <- round(seq(0,n_mats,length.out=n_tasks + 1))
idxmin <- idx[task]+1
idxmax <- idx[task+1]

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

for (i in idxmin:idxmax) {
	dist_i <- read.csv(paste0('/mnt/ls15/scratch/users/qdr/distmats/sortedmat', i, '.csv'))
	sisters_df <- identify_sisters(dist_i)
	write.csv(sisters_df, file = paste0('~/verts/sisterlists/truesisters', i, '.csv'), row.names = FALSE)
}
