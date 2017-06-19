# Sort distance matrices. Later we will have to load them again.
# Save each one as a separate csv.
# Done in parallel.

library(ape)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

fileno <- 70:72
distno <- c(99,100,100)
matno <- c(7099, 7100, 7200)

#load('/mnt/ls15/scratch/users/qdr/distmats/distmat1.r')
#roworder <- dimnames(ericsondist[[1]])[[1]]
#save(roworder, file = '~/verts/code/roworder.r')

load(paste0('/mnt/ls15/scratch/users/qdr/distmats/distmat',fileno,'.r'))
load('~/verts/code/roworder.r')
write.csv(ericsondist[[distno]][roworder, roworder], file = paste0('/mnt/ls15/scratch/users/qdr/distmats/sortedmat',matno,'.csv'), row.names = TRUE)
