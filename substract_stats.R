#!/usr/bin/env Rscript

library(data.table)
library(matrixStats)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process_vdj/in-vivo/wt/86765_B18_GCB4.aln.point.stats', '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process_vdj/in-vivo/wt/86759_B18_AIDKO_GCB2.aln.point.stats', '~/test.substract.stats')

PDA <- fread(args[1], colClasses=c(chr="character", pos="integer", count="numeric", depth="numeric"))   # left side of substraction
PDB <- fread(args[2], colClasses=c(chr="character", pos="integer", count="numeric", depth="numeric"))   # right side of substraction
out <- args[3]          # destination file

PDC <- merge(PDA, PDB, by=c('chr', 'pos', 'type'), suffixes=c('.a', '.b'), all=TRUE)
PDC[is.na(count.a), count.a := 0]
PDC[is.na(count.b), count.b := 0]
# The depth is the same for all outcomes at each given position. There should at least be a depth entry for the reference base at each position. So really there are at most two values: a (possibly repeated) depth, and possibly NA introduced by the merger. Using max(na.rm=TRUE) to propagate the real depth to the NAs.
suppressWarnings( PDC[, depth.a := max(depth.a, na.rm=TRUE), by=.(chr, pos)] ) # Generates warning if all values are NA (ie either table is missing a bit at an end)
suppressWarnings( PDC[, depth.b := max(depth.b, na.rm=TRUE), by=.(chr, pos)] )
PDC[!is.finite(depth.a), depth.a := 0]    # Overhanging ends
PDC[!is.finite(depth.b), depth.b := 0]

# Scale counts down to the sample with the lowest coverage.
# Because of the base quality filter in the pileup, each position has different covereage from which the count is derived.
# Normalizing to the sample-wide read count would be unrepresentative for positions with local coverage dips that may not be reflected in the other sample.
PDC[, goal := rowMins(as.matrix(PDC[, .(depth.a, depth.b)]))]
PDC[, count.a.scaled := (count.a / depth.a) * goal]
PDC[, count.b.scaled := (count.b / depth.b) * goal]
PDC[, depth.a.scaled := (depth.a / depth.a) * goal]
PDC[, depth.b.scaled := (depth.b / depth.b) * goal]
# Divisions by 0 can occur when either depth is 0. There can be no meaningful substracted mutation frequency for those.
PDC[!is.finite(count.a.scaled), count.a.scaled := 0]
PDC[!is.finite(count.b.scaled), count.b.scaled := 0]
PDC[!is.finite(depth.a.scaled), depth.a.scaled := 0]
PDC[!is.finite(depth.b.scaled), depth.b.scaled := 0]


# Set count and depth for the substracted track. A minus B
PDC[, count := floor(count.a.scaled - count.b.scaled)]
PDC[, depth := floor(rowMins(as.matrix(PDC[, .(depth.a.scaled, depth.b.scaled)])))] # lowest coverage between the two: Only as strong as the weakest link.

PDC <- PDC[, .(chr, pos, type, count, depth)]


fwrite(PDC, file=out, quote=FALSE, sep="\t")
