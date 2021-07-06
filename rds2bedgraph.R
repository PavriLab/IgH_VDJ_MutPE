#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('/Volumes/groups/pavri/Kimon/ursi/B18_MutPE/round7/results/invitro/tracks', '/Volumes/groups/pavri/Kimon/ursi/B18_MutPE/round7/process/invitro/pileup/134329_B18_d4_1.aln.1-1.point.stats.RDS', '/Volumes/groups/pavri/Kimon/ursi/B18_MutPE/round7/process/invitro/pileup/134329_B18_d4_1.aln.2-2.point.stats.RDS')
outpref <- args[1]
files <- args[2:length(args)]


###########
# bedGraph
###########

for (sf in files) {
	# sf <- files[1]
	
	rds <- readRDS(sf)
	outfile <- file.path(outpref, sub('stats.RDS', 'bedGraph', basename(sf)))
	# Create bedGraph of mutations
	# Use the normal mapped coordinates for the bedGraph, not the shifted ones used for peak labelling in the plots.
	bG <- unique(rds$posdata[(mutated), .(chr, pos-1, pos, aggrfreq)])
	if (dim(bG)[1] != 0) {
		cat(paste0("track type=bedGraph name=", sub('.stats.RDS', '', basename(sf)), "\n"), file=outfile)
		fwrite(bG, file=outfile, col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE, append=TRUE)
	} else {
		message("No mutations to put in the bedGraph.")
	}
}


