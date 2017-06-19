# nested for loops to calculate summary statistics from 10000 distance matrices, loading only one row of each matrix at a time.
# Done in parallel

n_tasks <- 250 # approx. 40 rows per task.
n_rows <- 9994
n_mats <- 10000

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
idx <- round(seq(0,n_rows,length.out=n_tasks + 1))
idxmin <- idx[task]+1
idxmax <- idx[task+1]
n_rows_task <- length(idxmin:idxmax)

ericsondistmean <- matrix(0, nrow=n_rows_task, ncol=n_rows)
ericsondistlower <- matrix(0, nrow=n_rows_task, ncol=n_rows)
ericsondistupper <- matrix(0, nrow=n_rows_task, ncol=n_rows)

pb <- txtProgressBar(0, n_rows_task, style = 3)

for (i in 1:n_rows_task) {
	rowidx_i <- idxmin + i - 1
	setTxtProgressBar(pb, i)
	
	each_row_i <- matrix(0, nrow=n_mats, ncol=n_rows)
	for (j in 1:10000) {
		# Load row i of each matrix in succession.
		each_row_i[j,] <- as.numeric(read.csv(paste0('/mnt/ls15/scratch/users/qdr/distmats/sortedmat', j, '.csv'), header = FALSE, skip = rowidx_i, nrows = 1))
	} 
	# Calculate summary stats for each row.
	ericsondistmean[rowidx_i, ] <- apply(each_row_i, 2, mean)
	ericsondistlower[rowidx_i, ] <- apply(each_row_i, 2, quantile, probs = 0.025)
	ericsondistupper[rowidx_i, ] <- apply(each_row_i, 2, quantile, probs = 0.975)
}

close(pb)

# load the header names and assign names to the summary matrices.
# header <- read.csv('/mnt/ls15/scratch/users/qdr/distmats/sortedmat1.csv', header=TRUE, nrows=1)
# dnames <- list(dimnames(header)[[1]], dimnames(header)[[1]])
# dimnames(ericsondistmean) <- dnames
# dimnames(ericsondistlower) <- dnames
# dimnames(ericsondistupper) <- dnames

save(ericsondistmean, ericsondistlower, ericsondistupper, file = paste0('~/verts/matslices/distmeans',task,'.r')