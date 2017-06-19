# nested for loops to calculate summary statistics from 10000 distance matrices, loading only one row of each matrix at a time.
# Done in parallel
# Version 2: load ~40 rows at a time instead of 1 row at a time in 40 iterations of a loop.

n_tasks <- 250 # approx. 40 rows per task.
n_rows <- 9994
n_mats <- 10000

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
idx <- round(seq(0,n_rows,length.out=n_tasks + 1))
idxmin <- idx[task]+1
idxmax <- idx[task+1]
n_rows_task <- length(idxmin:idxmax)

ericsondist <- list()

pb <- txtProgressBar(0, n_mats, style = 3)

	for (j in 1:n_mats) {
		pb <- setTxtProgressBar(pb, j)
		# Load ~40 rows of each matrix in succession.
		ericsondist[[j]] <- read.csv(paste0('/mnt/ls15/scratch/users/qdr/distmats/sortedmat', j, '.csv'), header = FALSE, skip = idxmin, nrows = n_rows_task)
	} 
	
close(pb)
	
	# bind ericson mats into an array
	ericsondist <- abind::abind(ericsondist, along = 3)
	
	# Calculate summary stats for each of the ~40 rows.
	ericsondistmean <- apply(ericsondist, 1:2, mean)
	ericsondistlower <- apply(ericsondist, 1:2, quantile, probs = 0.025)
	ericsondistupper <- apply(ericsondist, 1:2, quantile, probs = 0.975)
	
# load the header names and assign names to the summary matrices.
# header <- read.csv('/mnt/ls15/scratch/users/qdr/distmats/sortedmat1.csv', header=TRUE, nrows=1)
# dnames <- list(dimnames(header)[[1]], dimnames(header)[[1]])
# dimnames(ericsondistmean) <- dnames
# dimnames(ericsondistlower) <- dnames
# dimnames(ericsondistupper) <- dnames

save(ericsondistmean, ericsondistlower, ericsondistupper, file = paste0('~/verts/matslices/distmeans',task,'.r')