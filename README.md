# MutPE

This pipeline aligns paired-end reads from amplicon sequencing to the respective amplicon reference and plots mismatch frequencies per position.

The pipeline is structured as a Bash script that calls other scripts and submits jobs to the SLURM scheduler. Within the script, different sets of tasks are grouped as modules that can be enebled or disabled to facilitate partial reruns.

Certain file paths and property values are hard-coded in the pipeline and/or its scripts and can be selected for a predefined set of VDJ contexts. New sets of predefined values will need to be created and added to the code in the relevant places for use with new VDJ contexts.

# Dependencies

Dependencies are assumed to be loaded in advance in the execution environment, there is no management of dependencies included in the pipeline.

* HPC with SLURM scheduler
* python 3.6.7, with pandas, pysam, biopython, levenshtein
* R 3.5.1, with data.table, ggplot2, ggrepel, patchwork
* fastqc 0.11.8
* multiqc 1.7
* cutadapt 2.5
* trimmomatic 0.39
* flash 1.2.11 (the read-pair merging tool, not Adobe Flash!)
* bowtie2 2.3.5
* samtools 1.9
