#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)


args <- commandArgs(trailingOnly = TRUE)
# args <- c('NULL', '~/test.substract', '~/test.substract.stats.RDS')

if (args[1] == 'NULL') {
  hardymax <- NULL
} else {
  hardymax <- as.numeric(args[1])  # override ymax across all inputs
}
outpref <- args[2]              # file /path/prefix. Ending will be appended automatically.
files <- args[3:length(args)]   # one or more preprocess .RDS with contingency table plots and position plots already prepared.


#####################
# Publication layout
#####################

pdf(file=paste0(outpref, '.freqs_layout.pdf'), height=6, width=12)

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
  if(! any('catcol' == names(rds))) {
    rds$catcol <- c("C:G in WRCH/DGYW"="#DD0000", 
                    "C:G in WRCH/DGYW and SYC/GRS"="black", 
                    "C:G in SYC/GRS"="black", 
                    "other C:G"="black", 
                    "A:T"="grey60")
  }
  
  rds$posdata[, mutcat := factor(mutcat, ordered=TRUE, levels=rds$catlev)]
  
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
  
  
  #######################
  # Summary frequencies #
  #######################
  
  # Reshape overall mutation frequencies
  cont1 <- melt(as.data.table(rds$mutcntg), id.vars='mut', variable.name='ation', value.name='Frequency')
  setnames(cont1, c('to', 'from', 'Frequency'))
  cont1[, to := sub('to_', '', to, fixed=TRUE)]
  cont1[, from := sub('from_', '', from, fixed=TRUE)]
  cont1[, mutation := paste(from, to, sep='>')]
  cont1 <- cont1[!is.na(Frequency) & Frequency > 0,]
  cont1[, mutation := ordered(mutation, levels=rds$mutlev)]
  cont1[, Frequency := Frequency/sum(Frequency)]
  # Plot
  pcont1 <- ggplot(cont1, aes(x=mutation, y=Frequency*100, fill=mutation)) +
    geom_bar(stat='identity') +
    geom_hline(yintercept=0) +
    scale_y_continuous(expand=c(0,0)) +  # breaks=c(0,1,5,10,20,30,40,50,75,100) # for sqrt scaled only
    scale_fill_manual(values=rds$mutcol) +
    labs(x='', y='% of mutations') +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
          legend.position='none',
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=9) )
  
  # Reshape AIDER target frequencies
  cont2 <- as.data.table(t(rds$patcntg))
  names(cont2) <- 'Frequency'
  cont2[, Category := names(rds$patcntg)]
  cont2[Category=='hot_CG', Category := 'C:G in WRCH/DGYW']
  cont2[Category=='cold_CG', Category := 'C:G in SYC/GRS']
  cont2[Category=='other_CG', Category := 'other C:G']
  cont2[Category=='ambiguous_CG', Category := 'C:G in WRCH/DGYW and SYC/GRS']
  cont2[Category=='AT', Category := 'A:T'] 
  cont2[, Category := ordered(Category, levels=rds$catlev)]
  cont2[, Frequency := Frequency/sum(Frequency)]
  # Need some helper fields to help place bar percentages in the stacked bars.
  setorder(cont2, -Category)
  cont2[, lead := shift(Frequency,1)]
  cont2[is.na(lead), lead := 0]
  cont2[, cumfreq := cumsum(Frequency)]
  cont2[, texty := cumfreq - (Frequency / 2)]
  # Plot
  pcont2 <- ggplot(cont2, aes(x='', y=Frequency * 100, fill=Category) ) +
    geom_bar(stat='identity', position='stack',) +
    geom_text(aes(y=texty*100, label=round(Frequency*100,1), colour=Category), size=rel(2.5)) +
    scale_y_continuous(expand=c(0, 0)) +
    scale_fill_manual(values=rds$catcol) +
    scale_colour_manual(values=c('C:G in WRCH/DGYW'='white', 'other C:G'='white', 'A:T'='black', 'C:G in SYC/GRS'='black', 'C:G in WRCH/DGYW and SYC/GRS'='black')) +
    labs(x='', y='% of mutations') +
    theme_bw() +
    theme(legend.position='right',
          axis.ticks.length.x = unit(0, 'cm'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text=element_text(size=9),
          axis.title=element_text(size=9),
          legend.text=element_text(size=7),
          legend.title=element_blank())
  
  #######################################
  # Positions sorted by mutation frequency
  #######################################
  
  # Decreasing
  clump <- unique(rds$posdata[, .(pos, aggrfreq, mutcat)])
  setorder(clump, -aggrfreq, na.last = TRUE)
  clump[, pos := factor(pos, ordered=TRUE, levels=pos)]
  
  pclump <- ggplot(clump, aes(x=pos, y=aggrfreq, fill=mutcat)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=rds$catcol) +
    labs(x='Reordered positions', y='Mutation Frequency') +
    theme_minimal() +
    theme(legend.position='bottom',
          legend.title=element_blank(),
          legend.text = element_text(size=8),
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_text(size=9),
          axis.ticks.length.x = unit(0, 'cm'))
  
  
  ################
  # Position plot - cropped to positives only, no coverage
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
  ptrack <- pplain +
    geom_bar(data=unique(rds$posdata[, .(plotpos, aggrfreq, mutcat)]),    # Crop negative frequancies in the substractions.
             aes(x=plotpos, y=aggrfreq, fill=mutcat, colour=mutcat), stat='identity', position='identity', size=0, width=1) +
    scale_fill_manual(values=rds$catcol) +
    scale_colour_manual(values=rds$catcol) +
    coord_cartesian(xlim=c(0, rds$length), ylim=c(yext, ymax)) +
    labs(x='Position', y='Mutation Frequency')
  
  #######
  # Done
  #######
  
  print( 
      (ptrack + theme(legend.position='none')) / 
        ((pclump + theme(legend.position=c(0.8, 0.8))) | 
           ((pcont2 + theme(legend.position='none')) |
              pcont1 + labs(caption=basename(sf))) )
  )
}
dev.off()




