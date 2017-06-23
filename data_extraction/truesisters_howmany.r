# code to see how many of each of the lists the pairs are in.

all_lists <- list()
pb <- txtProgressBar(0,10000,style=3)

for (i in 1:10000) {
	all_lists[[i]] <- read.csv(paste0('~/verts/sisterlists/truesisters', i, '.csv'), stringsAsFactors = FALSE)
	setTxtProgressBar(pb,i)
}

close(pb)

all_lists <- do.call('rbind', all_lists)


library(dplyr)
sister_summary <- all_lists %>%
	mutate(bothnames = paste(sister1, sister2, sep = '_')) %>%
	group_by(bothnames, sister1, sister2) %>%
	summarize(proportion = n()/10000,
			  meandist = mean(dist),
			  mediandist = median(dist),
			  dist025 = quantile(dist, probs = 0.025),
			  dist975 = quantile(dist, probs = 0.975)) %>%
	ungroup %>%		  
	select(-bothnames)
	
write.csv(sister_summary, file = '~/verts/sister_summary.csv', row.names = FALSE)
	