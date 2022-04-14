#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)


args <- commandArgs(trailingOnly = TRUE)
# args <- c('0.97', 'auto', 'auto', '~/test.substract', '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round6/process/HDR1/pileup/101095_B18_HDR1_NP_e1r.aln.point.stats.RDS')

flagtop <- as.numeric(args[1])  # 1-a top percentile of frequencies for which to add position markers
outpref <- args[4]              # file /path/prefix. Ending will be appended automatically.
files <- args[5:length(args)]   # one or more preprocess .RDS with contingency table plots and position plots already prepared.


if (args[2] == 'NULL') {           # ymin auto-determined for each input
	hardymin <- NULL
} else if (args[3] == 'auto') {    # ymin shared across all inputs, but determined automatically such that no peaks are cropped.
	hardymin <- 0
	for (sf in files){
		rds <- readRDS(sf)
		hardymin <- min(hardymin, min(rds$posdata[, aggrfreq]))  # 0 or negative, never positive
	}
} else {
	hardymin <- as.numeric(args[3])  # ymin shared across all inputs and set to a manually chosen value. Probably to enable cropping negative values without causing problems with the annotations.
}

if (args[3] == 'NULL') {					 # ymax auto-adjusted for each input
	hardymax <- NULL
} else if (args[3] == 'auto') {    # ymax shared across all inputs, but determined automatically such that no peaks are cropped.
	hardymax <- 0
	for (sf in files){
		rds <- readRDS(sf)
		hardymax <- max(hardymax, max(rds$posdata[, aggrfreq]))  # 0 or positive, never negative
	}
} else {
	hardymax <- as.numeric(args[3])  # ymax shared across all inputs and set to a manually chosen value (allows cropping of peaks)
}


##############
# Ursi layout 
##############

pdf(file=paste0(outpref, '.labels_layout.pdf'), height=3, width=8)

for (sf in files) {
	# sf <- files[1]
	message(sf)

	rds <- readRDS(sf)
	
	####################
	# Order and colours
	####################
	
	# Reuse order and colour coding from contingency table plots if present.
	if(! any('catlev' == names(rds))) {
		rds$catlev <- c("C:G in WRCH/DGYW", 
		#"C:G i'n WRCH/DGYW and SYC/GRS", 
		#"C:G i'n SYC/GRS",
		"other C:G", 
		"A:T")
	}
	if(! any('catcol' == names(rds))) {
		rds$catcol <- c("C:G in WRCH/DGYW"="#DD0000", 
						#"C:G in WRCH/DGYW and SYC/GRS"="black", 
						#"C:G in SYC/GRS"="black", 
						"other C:G"="black", 
						"A:T"="grey50")
	}
	
	rds$posdata[, mutcat := factor(mutcat, ordered=TRUE, levels=rds$catlev)]
	
	
	################
	# Position plot - uncropped, with coverage
	################
	
	# Determine Y range, instead of extending annotations beyond the plot area and causing PDF/Illustrator ploblems.
	if (!is.null(hardymin)) {
		ymin <- hardymin
	} else {
		ymin <- min(0, max(rds$posdata[, aggrfreq]))  # 0 if no negative values are present.
	}
	if (!is.null(hardymax)){
		ymax <- hardymax
	} else {
		ymax <- max(rds$posdata[, aggrfreq])
	}
	
	
	# Extend Y range to make room for annotations
	annotheight <- (ymax - ymin) * 0.1        # Scale the height of the annotation to be visually constant regardless of Y range
	nannot <- 2                               # Number of annotation tracks
	yext <- ymin - (nannot * annotheight)     # extended to include space for the annotation tracks
	ymax <- 1.1 * ymax                        # Make a bit of room above
	
	# Determine top frequencies to flag
	maxdepth <- rds$posdata[, max(depth, na.rm=TRUE)]
	topthresh <- quantile(rds$posdata[, aggrfreq], flagtop)
	rds$posdata[, istop := aggrfreq >= topthresh]
	rds$posdata[, scaledepth := (depth / maxdepth) * ymax]
	stopifnot(all.equal(max(rds$posdata[, scaledepth]), ymax))
	# Base plot
	pp <- ggplot(data.frame(x=c(rds$xmin, rds$xmax), y=c(yext, ymax)), 
							 aes(x, y)) +       # token data to create area
		scale_x_continuous(expand=c(0,1), breaks=seq(rds$xmin, rds$xmax, 50)) +      # xmin already incorporates the -1 so as to be 0 when the first position is 1labs(x="Position") +
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
	pflag <- pplain +
		geom_line(data=unique(rds$posdata[, .(plotpos, scaledepth)]),
							aes(x=plotpos, y=scaledepth), inherit.aes=FALSE, colour='steelblue3', size=rel(0.5)) +
		geom_bar(data=unique(rds$posdata[, .(plotpos, aggrfreq, mutcat)]),    # Crop negative frequancies in the substractions.
						 aes(x=plotpos, y=aggrfreq, fill=mutcat, colour=mutcat), stat='identity', position='identity', size=0, width=1) +
		geom_label_repel(data=unique(rds$posdata[(istop), .(plotpos, aggrfreq, mutcat)]),
										 aes(x=plotpos, y=aggrfreq, label=plotpos, colour=mutcat), inherit.aes=FALSE, min.segment.length=0, show.legend=FALSE, size=rel(1.8)) +
		scale_y_continuous(sec.axis = sec_axis(trans= ~ . * maxdepth / ymax, name = "Coverage")) +
		scale_fill_manual(values=rds$catcol) +
		scale_colour_manual(values=rds$catcol) +
		coord_cartesian(xlim=c(0, rds$length), ylim=c(yext, ymax)) +
		labs(x= '', y='Mutation Frequency', caption=basename(sf)) + theme(legend.position='none')
	
	
	#######
	# Done
	#######
	
	print( pflag )
}
dev.off()

