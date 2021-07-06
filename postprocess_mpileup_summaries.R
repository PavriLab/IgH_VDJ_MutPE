#!/usr/bin/env Rscript

library(data.table)


args <- commandArgs(trailingOnly = TRUE)
# args <- c('~/test.substract.stats', 'B18', 'na:0:0', '0.3')
# args <- c('/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process_vdj/in-vitro/wt/86746_B18_AIDER_r2.aln.point.stats', 'B18', '-:0:0', '0.3')

sf <- args[1]                       # .stats input file: seq \t pos \t type \t count \t depth
vdj <- args[2]                      # B18, B18_1a, B18_1b, CH12, CH12_1, CH12_1a, CH12_1b, Ramos
offset <- args[3]                   # Correct for deletions in references. "REF:START:LENGTH,REF:START:LENGTH" ie. "HDR2:280:6"
                                    # For no corrections needed you can use "-:0:0". REF will be grepl'ed, so avoid REF names with dashes. 
                                    # "na:0:0" will be auto-converted to "-:0:0". This is to enable nesting bash commands (ie. for cluster submission), as the dash gets misinterpreted for a flag and needs to use up quotations.
allelecutoff <- as.numeric(args[4]) # Ignore mutations more frequent that this (ie. clonal mutations).


# Deconvolve the corrections string.
if (offset == "na:0:0")
	offset <- "-:0:0"
offset <- strsplit(offset, ':', fixed=TRUE)[[1]]

###########
# Presets
###########

# View ranges
if (vdj == 'B18' || vdj == 'B18_1') {
	showfrom <- 2022
	showto <- 2381
} else if (vdj == 'B18_1a') {
	showfrom <- 2022
	showto <- 2209
} else if (vdj == 'B18_1b') {
	showfrom <- 2195
	showto <- 2381
} else if (vdj == 'CH12' || vdj == 'CH12_1') {
	showfrom <- 1694
	showto <- 2059
} else if (vdj == 'CH12_1a') {
	showfrom <- 1694
	showto <- 1904
} else if (vdj == 'CH12_1b') {
	showfrom <- 1874
	showto <- 2077
} else if (vdj == 'Ramos1') {
	showfrom <- 3684
	showto <- 3929
} else if (vdj == 'Ramos2') {
	showfrom <- 3927
	showto <- 4147
} else if (vdj == 'Ramos3') {
	showfrom <- 4148
	showto <- 4375
} else if (vdj == 'Ramos') {
	showfrom <- 3684
	showto <- 4375
} else {
	stop(paste('Unknown view window', vdj, '.'))
}

# # CDRs and hotspots (dev/debug paths)
# if (vdj == 'B18' || vdj == 'B18_1' || vdj == 'B18_1a' || vdj == 'B18_1b') {
# 	cdrs <- data.table(xleft=c(2112, 2169, 2300) - 0.3,
# 										 xright=c(2126, 2216, 2353) + 0.3 )
# 	bedfw <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_WRCH.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# 	bedrv <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_DGYW.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# 	bedfwn <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_SYC.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# 	bedrvn <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_GRS.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# } else if (vdj == 'CH12' || vdj == 'CH12_1' || vdj == 'CH12_1a' || vdj == 'CH12_1b') {
# 	cdrs <- data.table(xleft=c(1784, 1841, 1988) - 0.3,
# 										 xright=c(1798, 1891, 2026) + 0.3 )
# 	bedfw <- fread('/Volumes/groups/pavri/Kimon/ursi/cH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_WRCH.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# 	bedrv <- fread('/Volumes/groups/pavri/Kimon/ursi/cH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_DGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# 	bedfwn <- fread('/Volumes/groups/pavri/Kimon/ursi/cH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_SYC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# 	bedrvn <- fread('/Volumes/groups/pavri/Kimon/ursi/cH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_GRS.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# } else if (grepl('Ramos', vdj)) {
# 	cdrs <- data.table(xleft=c(3899, 3980, 4107) - 0.3,
# 										 xright=c(3943, 4021, 4181) + 0.3 )
# 
# 	bedfw <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_WRCH.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# 	bedrv <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_DGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# 	bedfwn <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_SYC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# 	bedrvn <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_GRS.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
# } else {
# 	stop('No BEDfile defined for the viewing window', vdj, '.')
# }

## CDRs and hotspots
if (vdj == 'B18' || vdj == 'B18_1' || vdj == 'B18_1a' || vdj == 'B18_1b') {
	cdrs <- data.table(xleft=c(2112, 2169, 2300) - 0.3,
										 xright=c(2126, 2216, 2353) + 0.3 )
	bedfw <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_WRCH.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	bedrv <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_DGYW.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	bedfwn <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_SYC.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	bedrvn <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_GRS.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
} else if (vdj == 'CH12' || vdj == 'CH12_1' || vdj == 'CH12_1a' || vdj == 'CH12_1b') {
	cdrs <- data.table(xleft=c(1784, 1841, 1988) - 0.3,
										 xright=c(1798, 1891, 2026) + 0.3 )
	bedfw <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_WRCH.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	bedrv <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_DGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	bedfwn <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_SYC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	bedrvn <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_GRS.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	
} else if (grepl('Ramos', vdj)) {
	cdrs <- data.table(xleft=c(3899, 3980, 4107) - 0.3,
										 xright=c(3943, 4021, 4181) + 0.3 )
	bedfw <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_WRCH.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	bedrv <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_DGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	bedfwn <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_SYC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
	bedrvn <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_GRS.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
} else {
	stop('No BEDfile defined for the viewing window', vdj, '.')
}

# Clean up BED
names(bedfw) <- c('chr', 'xleft', 'xright', 'seq')
names(bedrv) <- c('chr', 'xleft', 'xright', 'seq')
names(bedfwn) <- c('chr', 'xleft', 'xright', 'seq')
names(bedrvn) <- c('chr', 'xleft', 'xright', 'seq')
# Isolate core coordinate of the pattern
bedfw[, x := xleft + 3]    # wrCh in 1-based coordinate from 0 bases coordinate range
bedrv[, x := xleft + 2]    # dGyw in 1-based coordinate
bedfwn[, x := NA_real_]
bedrvn[, x := NA_real_]
# Combine
bedfw[, fw := TRUE]
bedrv[, fw := FALSE]
bedfwn[, fw := TRUE]
bedrvn[, fw := FALSE]
bedfw[, origin := 'hot']
bedrv[, origin := 'hot']
bedfwn[, origin := 'cold']
bedrvn[, origin := 'cold']
bed <- rbind(bedfw, bedrv, bedfwn, bedrvn) # remove overlaps (same hot C/G base)
setorder(bed, xleft)

# Convert coordinate ranges from 0-based right-open to 1-based right-closed
bed[, xleft := xleft + 1] 
# Pad segment ends to more closely line up with bar edges rather than centres
bed[, xleft := xleft - 0.3]
bed[, xright := xright + 0.3]
# Distinguish between strict and loose pattern.
# The seq field is already rev'comp'ed for reverse strand matches, by the way the BEDfiles were made.
hotspots <- bed[origin == 'hot' & grepl('[Aa][Gg][Cc][Tt]', seq, perl=TRUE), ] [, origin := 'hot']
warmspots <- bed[origin == 'hot' & !grepl('[Aa][Gg][Cc][Tt]', seq, perl=TRUE), ] [, origin := 'warm']
coldspots <- bed[origin=='cold', ]
# Update designations.
bed <- rbind(hotspots, warmspots, coldspots)
# Prepare annotations. Lose duplicates from palindromes.
hotspots <- unique(bed[origin == 'hot', .(xleft, xright)])
warmspots <- unique(bed[origin == 'warm', .(xleft, xright)])
coldspots <- unique(bed[origin=='cold', .(xleft, xright)])

# Shift annotation coordinates so x axis starts at 1.
hotspots <- hotspots - showfrom + 1
warmspots <- warmspots - showfrom + 1
coldspots <- coldspots - showfrom + 1

# Create relative positions
bed[, plotx := x - showfrom + 1]
bed[, plotxleft := xleft - showfrom + 1]
bed[, plotxright := xright - showfrom + 1]
cdrs[, plotxleft := xleft - showfrom + 1]
cdrs[, plotxright := xright - showfrom + 1]


########
# Input
########

# Nothing to do if file is empty
info = file.info(sf)
if (info$size == 0) {
	message(paste("No content found in ", sf))
	q(status=0, save="no")
}

# Input
posdata <- fread(sf)

# Nothing to do if file is empty
if (dim(posdata)[1] == 0) {
	message(paste("No content found in ", sf))
	q(status=0, save="no")
}

# Unify headers between my new names and Tobi's older names
if (all(names(posdata) == c('amplicon', 'pos', 'type', 'freq', 'total')))
	names(posdata) <- c('chr', 'pos', 'type', 'count', 'depth')
stopifnot(all(names(posdata) == c('chr', 'pos', 'type', 'count', 'depth')))


##############
# Coordinates
##############

# Correct the coordinates to account for reference deletions
# Create the shifted position. Retain the original position.
posdata[, newpos := pos]
sel <- ( grepl(offset[1], posdata$chr) &
				 posdata$pos >= as.integer(offset[2]) )
posdata[(sel), newpos := pos + as.integer(offset[3])]

# Apply the viewing window
posdata <- posdata[newpos >= showfrom & newpos <= showto, ]

# Fill in gaps in the coordinates. Prevents plot lines from jumping across gaps.
setorder(posdata, chr, newpos)
posdata <- merge(posdata, data.table(newpos=showfrom:showto), by="newpos", all=TRUE)
posdata[is.na(chr), chr := posdata[!is.na(chr), chr][1]]
posdata[is.na(count), count := 0]
posdata[is.na(depth), depth := 0]

# Create relative positions
posdata[, plotpos := pos - showfrom + 1]


##############
# Mutations
##############

# Relabel reference entries and remove any indel entries
posdata[grepl('^[NACGT]$', type), type := 'ref']
posdata[, mutated := (type != 'ref')]
posdata <- posdata[! grepl('\\+|-', type, perl=TRUE), ]


##############
# Frequencies
##############

# Calculate frequencies 
posdata[(mutated), freq := count / depth]
posdata[!is.finite(freq), freq := 0]      # Missing count, or lack of coverage

# Remove mutations with high frequencies, as likely clonal SNPs.
n <- nrow(posdata[(mutated),])
posdata <- posdata[!mutated | (mutated & freq <= allelecutoff), ]
diffn <- n - nrow(posdata[(mutated),])
if (diffn > 0)
	message(paste0("Mutations removed as likely clonal SNPs (>=", allelecutoff, "): ", diffn, " / ", n))

# Aggregate frequencies by position
posdata[(mutated), aggrcount := sum(count), by=c('chr', 'newpos')]    # Aggrecated mutation count (all mutation types per position)
posdata[is.na(aggrcount), aggrcount := 0]
posdata[, aggrfreq := aggrcount / depth]


###########
# Patterns
###########

# Categorize mutations

# by source base type
posdata[(mutated), muttype := 'A:T']
posdata[mutated & (grepl('C>', type, fixed=TRUE) | grepl('G>', type, fixed=TRUE)), muttype := 'G:C']
# by presence on a warm or hot spot
posdata[(mutated), sweetspot := FALSE]
posdata[mutated & grepl('C>', type, fixed=TRUE) & pos %in% bed[fw & origin %in% c('hot', 'warm'), x], sweetspot := TRUE]
posdata[mutated & grepl('G>', type, fixed=TRUE) & pos %in% bed[(!fw) & origin %in% c('hot', 'warm'), x], sweetspot := TRUE]
posdata[(mutated), sourspot := FALSE]
posdata[mutated & grepl('C>', type, fixed=TRUE), 
				sourspot := vapply(posdata[mutated & grepl('C>', type, fixed=TRUE), pos], 
													 function(p) {
													 	any(bed[fw & origin == 'cold', p > xleft & p < xright]) # in position in the cold pattern
													 }, logical(1)) ] # Remember, the x have been padded with a decimal value
posdata[mutated & grepl('G>', type, fixed=TRUE), 
				sourspot := vapply(posdata[mutated & grepl('G>', type, fixed=TRUE), pos], 
													 function(p) {
													 	any(bed[(!fw) & origin == 'cold', p > xleft & p < xright])  # in position in the cold pattern
													 }, logical(1)) ] # Remember, the x have been padded with a decimal value
posdata[, mutcat := muttype]
posdata[muttype == 'G:C', mutcat := 'other C:G']
posdata[(sweetspot), mutcat := 'C:G in WRCH/DGYW']
# posdata[sourspot & !sweetspot, mutcat := 'C:G in SYC/GRS']      # not claimed by a hot/warm spot
# posdata[sourspot & sweetspot, mutcat := 'C:G in WRCH/DGYW and SYC/GRS']      # overlaps


#######
# Done
#######

rds <- list(vdj=vdj, 
			offset=offset, 
			showfrom=showfrom,
			length=showto-showfrom + 1,
			xmin=posdata[1, plotpos] - 1,           # In case filtering loses positions later, this gives me the X dimensions of the area.
			xmax=posdata[nrow(posdata), plotpos],   # In case filtering loses positions later, this gives me the X dimensions of the area.
			allelecutoff=allelecutoff,
			posdata=posdata[(mutated), ], 
			bed=bed, 
			cdr=cdrs)

# Save as R list object, so that I have all the related info in one place for later plotting, without spamming with separate files.
saveRDS(rds, file=paste0(sf, '.RDS'))
