
Instructions for MutPE-seq analyses

### 1. System requirements

The pipelines are designed to run on a HPC cluster with Linux operating system and SLURM workload manager: https://slurm.schedmd.com

### 2. Installation guide

The pipelines software dependencies can be installed with conda:

1. install conda with the bioconda channel: https://bioconda.github.io/

2. create a new env with the required packages:

`conda create -n igh_pipe bc umi_tools r-data.table fastqc multiqc cutadapt trimmomatic bowtie bedtools samtools ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph pysam pandas intervaltree pybedtools ucsc-bedsort bowtie2 flash r-patchwork python-levenshtein biopython r-matrixstats r-ggrepel`

4. clone the pipeline repository locally: https://github.com/PavriLab/IgH_VDJ_MutPE

### 3. Instructions for use

Overview of the pipeline: [MutPESeq_Pipeline.pptx](MutPESeq_Pipeline.pptx)

Before running the pipeline, it is necessary to define which reference will be used to align the reads. The repositories contain all the references used for the analyses in the paper on the `reference_and_annotation` folder.

The pipeline is run with the bash script `mutPE_workflow.sh` contained within the code repository. Assuming all input fastq are stored in `input_data/example/<x>_1.fastq.gz`, `input_data/example/<x>_2.fastq.gz`, where `<x>` is the sample name, for example:
```bash
REFLEN=361
REF=B1-8hi
REF_INDEX_PREFIX=data/B1-8hi
bash mutPE_workflow.sh -d "$(pwd)" -D input_data -H “$(pwd)” -r out -A 0.3 -R $REFLEN -V $REF -b example -i $REF_INDEX_PREFIX -S data/subs/example.txt -P -1 -2 -3 -4 -5
```

`data/subs/example.txt` contains the list of paired samples for background subtraction.

`$REF` is the name of the reference.

`$REF_INDEX_PREFIX` is the prefix of the reference bowtie2 index.

The output plots will be in the `out/example` directory.

See also: [NOTES.md](NOTES.md)