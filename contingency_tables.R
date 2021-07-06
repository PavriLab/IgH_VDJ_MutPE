#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('~/test.substract.stats.RDS')

for (sf in args) {
	# sf <- args[1]    # .stats.RDS post-processed mpileup summary
	
	rds <- readRDS(sf)
	
	
	# Deconvolve mutation type
	rds$posdata[grepl('>', type) & !grepl('N>', type), c('was', 'became') := list(substr(type, start=1, stop=1), substr(type, start=3, stop=3))]
	
	####################
	# Overall mutations
	####################
	
	contingency1 <- data.frame('mut'=c('to_A', 'T', 'G', 'C'),
														 from_A=c(0,
														 				 sum(rds$posdata[was == 'A' & became == 'T' & freq >= 0, count]),
														 				 sum(rds$posdata[was == 'A' & became == 'G' & freq >= 0, count]),
														 				 sum(rds$posdata[was == 'A' & became == 'C' & freq >= 0, count])),
														 
														 T=c(sum(rds$posdata[was == 'T' & became == 'A' & freq >= 0, count]),
														 		0,
														 		sum(rds$posdata[was == 'T' & became == 'G' & freq >= 0, count]),
														 		sum(rds$posdata[was == 'T' & became == 'C' & freq >= 0, count])),
														 
														 G=c(sum(rds$posdata[was == 'G' & became == 'A' & freq >= 0, count]),
														 		sum(rds$posdata[was == 'G' & became == 'T' & freq >= 0, count]),
														 		0,
														 		sum(rds$posdata[was == 'G' & became == 'C' & freq >= 0, count])),
														 
														 C=c(sum(rds$posdata[was == 'C' & became == 'A' & freq >= 0, count]),
														 		sum(rds$posdata[was == 'C' & became == 'T' & freq >= 0, count]),
														 		sum(rds$posdata[was == 'C' & became == 'G' & freq >= 0, count]),
														 		0) )
	# Fix number, order and colour of categories for point mutations and combined files.
	if(!any('mutlev' == names(rds))){
		rds$mutlev <- c("C>G", "C>A", "C>T", "C>N",
								"G>C", "G>A", "G>T", "G>N", 
								"A>T", "A>C", "A>G", "A>N", 
								"T>A", "T>C", "T>G", "T>N",
								"N>C", "N>G", "N>A", "N>T")
	}
	if(!any('mutcol' == names(rds))){
		rds$mutcol <- c("C>G"='black', "C>A"='black', "C>T"='black', "C>N"='black',
								"G>C"='black', "G>A"='black', "G>T"='black', "G>N"='black',
								"A>T"='grey30', "A>C"='grey30', "A>G"='grey30', "A>N"='grey30',
								"T>A"='grey30', "T>C"='grey30', "T>G"='grey30', "T>N"='grey30',
								"N>C"='grey60', "N>G"='grey60', "N>A"='grey60', "N>T"='grey60')
	}
	
	
	
	##########################
	# Potential AIDER targets
	##########################
	
	# contingency2 <- data.frame(hot_CG = sum(rds$posdata[sweetspot & !sourspot, count]),
	# 													 cold_CG = sum(rds$posdata[sourspot & !sweetspot, count]),
	# 													 ambiguous_CG = sum(rds$posdata[sourspot & sweetspot, count]),
	# 													 other_CG = sum(rds$posdata[!sweetspot & !sourspot & muttype=='G:C', count]),
	# 													 AT = sum(rds$posdata[muttype == 'A:T', count]) )
	contingency2 <- data.frame(hot_CG = sum(rds$posdata[sweetspot & aggrfreq >= 0, count]),
														 other_CG = sum(rds$posdata[(!sweetspot)  & aggrfreq >= 0 & muttype=='G:C', count]),
														 AT = sum(rds$posdata[muttype == 'A:T' & aggrfreq >= 0, count]) )
	
	# Fix number, order and colour of categories for point mutations and combined files.
	if(!any('catlev' == names(rds))){
		rds$catlev <- c("C:G in WRCH/DGYW", "C:G in WRCH/DGYW and SYC/GRS", "C:G in SYC/GRS", "other C:G", "A:T")
	}
	if(!any('catcol' == names(rds))){
		rds$catcol <- c("C:G in WRCH/DGYW"="#DD0000", 
										"C:G in WRCH/DGYW and SYC/GRS"="black", 
										"C:G in SYC/GRS"="black", 
										"other C:G"="black", 
										"A:T"="grey60")
	}
	
	
	#######
	# Done
	#######
	
	# Write out contingency tables
	mtxf <- sub('.RDS', '.mtx.txt', sf, fixed=TRUE)
	
	cat(c(basename(sf), "\n\n"), file=mtxf, append=FALSE)
	suppressWarnings( write.table(contingency2, file=mtxf, sep="\t", row.names=FALSE, quote=FALSE, append=TRUE) ) # appending with headers triggers warning
	cat("\n", file=mtxf, append=TRUE)
	suppressWarnings( write.table(contingency1, file=mtxf, sep="\t", row.names=FALSE, quote=FALSE, append=TRUE) )
	cat("\n", file=mtxf, append=TRUE)
	
	
	# Update the R object
	rds[['mutcntg']] <- contingency1
	rds[['patcntg']] <- contingency2
	
	saveRDS(rds, file=sf)
}

