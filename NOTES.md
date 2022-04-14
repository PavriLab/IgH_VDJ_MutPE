# Notes on the pipeline

## Genomic indices

The VDJs are added as an extra "chromosome" to the stock reference genome. The same sequences are used as for procap/proseq, but the indices are mad ewith bowtie2 here.

## Presets

`-V` controls preset values.

* In `mutPE_workflow.sh` this controls what adaptor sequences are searched and trimmed.

* In `postprocess_mpileup_summaries.R` it controls the coordinates of the VDJ amplicon relative to the VDJ "chromosome". 
It also defines the CDR coordinates and which hotspot annotation file is used.

* In `combine_amplicons.R` it also defines the amplicon coordinates relative to the "chromosome", for the purpose of concatenating the consecutive amplicons of older runs.


## Hotspot annotation

This is created by representing the hotspot/coldspot motives as regex (in fwd and in reverse-complemented form) and scanning the VDJ "chromosome" for pattern matches, then formatting the coordinates as a bedfile.
`sequtilities.py` is not part of the pipeline but implements the `regex2bed` functionality for this.
