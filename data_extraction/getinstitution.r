# Get institution IDs
vn <- read.csv('vertnet_latest_birds.csv', stringsAsFactors = FALSE)
vnids <- vn[,c('organismid','institutionid','institutioncode')]
write.csv(vnids, '~/verts/institution_ids.csv', row.names = FALSE)
