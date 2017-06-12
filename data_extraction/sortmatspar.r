# Sort distance matrices. Later we will have to load them again.
# Save each one as a separate csv.
# Done in parallel.

library(ape)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

#load('/mnt/ls15/scratch/users/qdr/distmats/distmat1.r')
#roworder <- dimnames(ericsondist[[1]])[[1]]
#save(roworder, file = '~/verts/code/roworder.r')

load(paste0('/mnt/ls15/scratch/users/qdr/distmats/distmat',task,'.r'))
load('~/verts/code/roworder.r')
for (j in 1:100) {
	write.csv(ericsondist[[j]][roworder, roworder], file = paste0('/mnt/ls15/scratch/users/qdr/distmats/sortedmat',j + (task-1)*100,'.csv'), row.names = TRUE)
}
