#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)
# args <- c('NULL', '~/test.substract', '~/test.substract.stats.RDS')
# args <- c('NULL', '~/test.substract', '/Volumes/groups/pavri/Kimon/rushad/mutpe/stats/86746_B18_AIDER_r2.aln.point.stats.RDS')

if (args[1] == 'NULL') {
	hardymax <- NULL
} else {
	hardymax <- as.numeric(args[1])  # override ymax across all inputs
}
outpref <- args[2]              # file /path/prefix. Ending will be appended automatically.
files <- args[3:length(args)]   # one or more preprocess .RDS with contingency table plots and position plots already prepared.


################
# No-A:T layout
################

pdf(file=paste0(outpref, '.noAT_layout.pdf'), height=3, width=8)

for (sf in files) {
	# sf <- files[1]
	rds <- readRDS(sf)
	
	
	####################
	# Order and colours
	####################
	
	# Reuse order and colour coding from contingency table plots if present.
	if(! any('catlev' == names(rds))) {
		rds$catlev <- c("C:G in WRCH/DGYW", "C:G in WRCH/DGYW and SYC/GRS", "C:G in SYC/GRS", "other C:G", "A:T")
	}
	# Override predefined colours
	rds$catcol <- c("C:G in WRCH/DGYW"="#DD0000", 
									"C:G in WRCH/DGYW and SYC/GRS"="black", 
									"C:G in SYC/GRS"="black", 
									"other C:G"="black", 
									"A:T"="transparent")
	
	rds$posdata[, mutcat := factor(mutcat, ordered=TRUE, levels=rds$catlev)]
	
	
	################
	# Position plot - cropped, without coverage, hiding A:T mutations
	################
	
	# Determine Y range, instead of extending annotations beyond the plot area and causing PDF/Illustrator ploblems.
	ymin <- 0
	if (!is.null(hardymax)){
		ymax <- hardymax
	} else {
		ymax <- max(rds$posdata[, aggrfreq])
	}
	
	# Extend Y range to make room for annotations
	annotheight <- (ymax - ymin) * 0.1        # Height of each annotation track
	nannot <- 2                               # Number of annotation tracks
	yext <- ymin - (nannot * annotheight)     # extended to include space for the annotation tracks
	ymax <- 1.1 * ymax                        # Make a bit of room above
	
	# Base plot
	pp <- ggplot(data.frame(x=c(rds$xmin, rds$xmax), y=c(yext, ymax)), 
							 aes(x, y)) +       # token data to create area
		scale_x_continuous(expand=c(0,1), breaks=seq(rds$xmin, rds$xmax, 50)) +      # xmin already incorporates the -1 so as to be 0 when the first position is 1
		labs(x="Position") +
		theme(legend.position = 'bottom',
					legend.title = element_blank(),
					legend.text = element_text(size=8),
					panel.grid.major.x = element_blank(),
					panel.grid.minor.x = element_blank(),
					axis.text = element_text(size=8),
					axis.title=element_text(size=9),
					panel.background = element_rect(fill="white"),
					panel.border = element_blank())
	# Add annotations
	pplain <- pp
	if (nrow(rds$cdr) > 0) {
		for (i in 1:nrow(rds$cdr)) {
			pplain <- pplain +
				geom_rect(xmin=rds$cdr[i, plotxleft], xmax=rds$cdr[i, plotxright], ymin=yext, ymax=ymax, fill='grey90', alpha=0.3)
		}
	}
	if (nrow(rds$bed[origin=='warm', ]) > 0) {
		pplain <- pplain +
			geom_rect(data=rds$bed[origin=='warm', ],
								aes(xmin=plotxleft, xmax=plotxright, ymin=ymin - (0.8 * annotheight), ymax=ymin - (0.2 * annotheight)), inherit.aes=FALSE, fill='black')
	}
	if (nrow(rds$bed[origin=='hot', ]) > 0) {
		pplain <- pplain +
			geom_rect(data=rds$bed[origin=='hot', ],
								aes(xmin=plotxleft, xmax=plotxright, ymin=ymin - (0.8 * annotheight), ymax=ymin - (0.2 * annotheight)), inherit.aes=FALSE, fill='red')
	}
	if (nrow(rds$bed[origin=='cold', ]) > 0) {
		pplain <- pplain +
			geom_rect(data=rds$bed[origin=='cold', ],
								aes(xmin=plotxleft, xmax=plotxright, ymin=ymin - (1.8 * annotheight), ymax=ymin - (1.2 * annotheight)), inherit.aes=FALSE, fill='black')
	}
	# Crop data nicely (no hidden objects out of frame, no overlaps with the annotation). DO NOT SAVE modified rds back onto itself!
	rds$posdata[aggrfreq > ymax, aggrfreq := ymax]
	rds$posdata[aggrfreq < ymin, aggrfreq := ymin]
	# Add data
	pnoat <- pplain +
		geom_bar(data=unique(rds$posdata[, .(plotpos, aggrfreq, mutcat)]),    # Crop negative frequancies in the substractions.
						 aes(x=plotpos, y=aggrfreq, fill=mutcat, colour=mutcat), stat='identity', position='identity', size=0, width=1) +
		scale_fill_manual(values=rds$catcol) +
		scale_colour_manual(values=rds$catcol) +
		coord_cartesian(xlim=c(rds$xmin, rds$xmax), ylim=c(yext, ymax)) +
		labs(x='', y='Mutation Frequency', caption=basename(sf)) + theme(legend.position='none')
	
	
	#######
	# Done
	#######
	
	print( pnoat )
}
dev.off()

