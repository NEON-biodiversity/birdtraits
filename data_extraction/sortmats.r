# Sort distance matrices. Later we will have to load them again.
# Save each one as a separate csv.

library(ape)

load('/mnt/ls15/scratch/users/qdr/distmats/distmat1.r')
roworder <- dimnames(ericsondist[[1]])[[1]]
ericsondist <- lapply(ericsondist, function(x) x[roworder, roworder])

matno <- 0
for (i in 1:100) {
	load(paste0('/mnt/ls15/scratch/users/qdr/distmats/distmat',i,'.r'))
	ericsondist <- lapply(ericsondist, function(x) x[roworder, roworder])
	for (j in 1:100) {
		matno <- matno + 1
		write.csv(ericsondist[[j]], file = paste0('/mnt/ls15/scratch/users/qdr/distmats/sortedmat',matno,'.csv'), row.names = TRUE)
	}
}