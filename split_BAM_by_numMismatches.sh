#!/usr/bin/env sh

# split_BAM_by_numMismatches.sh <path_to_script> <path/file> <MIN#> <MAX#>

samtools view -h $2 | perl ${1}/split_SAM_by_numMismatches.pl $3 $4 | samtools sort -O BAM -o ${2/.bam/.$3-$4.bam}
