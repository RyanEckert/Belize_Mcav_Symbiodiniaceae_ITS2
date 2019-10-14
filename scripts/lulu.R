#!/usr/bin/env Rscript

## Load libraries
library(lulu); packageVersion("lulu")

#-------------------------------------------------------------------------------
## Set up data and run LULU curation
# Let's take a look at how many ASVs we get by clustering at 84 and 97% similarity

alldat <- read.csv("ASVtable.txt",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
matchlist <- read.table("match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
curated_result_84 <- lulu(alldat, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)
curated_result_97 <- lulu(alldat, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 97, minimum_relative_cooccurence = 0.95)

#------------------------------------------------------------------------------
## Examine curated reads and save results for later
# Check out the curated results tables for the 97 and 84%

curated_result_84$curated_otus
curated_result_97$curated_otus

curated_result_84$curated_table
curated_result_97$curated_table

curated_result_84$otu_map
curated_result_97$otu_map

# write.csv(curated_result_97$otu_map, "lulu_otu_map_97.csv")

# Compare retained and discarded ASVs at 84 and 97% similarity in this case we retain 18 and 24 ASVs respectively.
# Print to file to view later
sink("curationResults.txt", append=F, split=T)
cat('There were',curated_result_84$curated_count,'ASVs retained, and',curated_result_84$discarded_count,'ASVs discarded when clustered at 84% similarity\n')
cat('There were',curated_result_97$curated_count,'ASVs retained, and',curated_result_97$discarded_count,'ASVs discarded when clustered at 97% similarity\n')
cat('\n')
cat('END TRANSMISSION')
sink()

# Save outputs from lulu to pick up analyses later if needed.
saveRDS(curated_result_84, file="curated_result_84.rds")
saveRDS(curated_result_97, file="curated_result_97.rds")

# If you need to read in the currated results, uncomment below and run
#curated_result_84 = readRDS("curated_result_84.rds")
#curated_result_97 = readRDS("curated_result_97.rds")

# Clear the R environment
rm(list=ls())
