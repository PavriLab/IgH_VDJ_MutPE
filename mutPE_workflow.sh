#!/usr/bin/env sh

echo "MutPE Workflow"
# Process the samples of one Run, from .bam/.fq

## Parameters ##
function usage() {
    echo "Usage (look in source for more parameters):"
    echo "      $0 -d base_Dir -b Batch_dir -i bowtie2_Index [-D Data_dir -p Processing_dir -r Results_dir -a Aux_dir] [-o Offset]"
		exit 1
}
# Defaults
user='kimon.fr'
data='data'
process='process'
results='results'
aux='aux'
offsets='na:0:0'
fastasuffix='fa'
XDIR='/groups/pavri/Kimon/mutpe'
refLen=361                      # fragment size
mergelen=30                     # +- nt from refLen for merge filtering
flash=1                         # merge read pairs
unpaired=0                      # ignore pairs and align single-end
ovlpmm=0.10						# FLASh : Maximum allowed mismatch rate when merging overlapping mates
maxOver=250                     # FLASh : max allowed overlap of mates
minOver=80                      # FLASh : min allowed overlap of mates
adapterr=0.1                    # cutadapt : mismatches in adapter identification
minL1=200                       # cutadapt : minimum allowed read1 length post-trim
minL2=100                       # cutadapt : minimum allowed read2 length post-trim
# clip1=-5                        # cutadapt : trim from read1
# clip2=-150                      # cutadapt : trim from read2
orient="--fr"                   # bowtie2 : paired-end orientation ("" / "--fr" / "--rf" / "--ff")
# minfrag=370                     # bowtie2 : paired-end max fragment size
# maxfrag=350                     # bowtie2 : paired-end min fragment size
minalq=13						# samtools : Minimum allowed read mapping quality for contributing to the pileup
minbasq=30						# samtools : Minimum allowed base quality for contributing to mutation stats
                                          # Reads are trimmed for base quality at 25 (hard-coded), so that's kind of a floor for $minbasq
dopre=1
doal=1
dopile=1
doviz=1
dopost=1
allowoutties=0
vdj="B18"                     # Default
allelecutoff=1                # No cutoff for allele frequencies

# Parse options.
while getopts 'a:A:b:d:D:f:F:H:i:l:L:m:M:o:Op:Pr:R:sS:UV:W:x:X:z12345' flag; do
  case "${flag}" in
    a) aux="${OPTARG}" ;;     		# Dir where the bowtie index is located, relative to base
    A) allelecutoff="${OPTARG}";; # Ignore mutations with frequencies higher than this, as likely clonal SNPs.
    b) run="${OPTARG}" ;;         # Batch subdir. If defined, it will be appended to $data, $process and $results
    # c) clip2="${OPTARG}" ;;
    # C) clip1="${OPTARG}" ;;
  	d) base="${OPTARG}" ;;        # Base dir
    D) data="${OPTARG}" ;;        # BAM Data dir relative to base
    f) mergelen="${OPTARG}";;
    F) fastasuffix="${OPTARG}";;  # Reference fasta suffix. (assumes the prefix is the same as the bowtie2 index prefix)
    H) XDIR="${OPTARG}" ;;        # Path to the mutpe scripts directory (ie. the local repo clone).
    i) bowtie2idx="${OPTARG}" ;;  # Bowtie2 index prefix
    l) minL2="${OPTARG}" ;;
    L) minL1="${OPTARG}" ;;
    m) minOver="${OPTARG}" ;;
    M) maxOver="${OPTARG}" ;;
    o) offsets="${OPTARG}" ;;     # Reference deletions: "REF:START:LENGTH" ie. "HDR2:2301:6,HDR4:2301:6"
    O) allowoutties=1;;           # also allow read pairs that overlap on their 5' ends instead of their 3' ends.
    p) process="${OPTARG}" ;;     # Dir in which intermediate steps output files are created, relative to base
    P) flash=0;;                  # paired-end instead of merged pairs
    r) results="${OPTARG}" ;;     # Dir for the final files, relative to base
    R) refLen="${OPTARG}";;
    s) issra=1 ;;                 # Input is .sra format
    S) subs="${OPTARG}" ;;       # file specifying the substraction patterns (from_prefix \t minus_prefix)
    U) unpaired=1;;                 # single-end. Only applicable if also -P
    V) vdj="${OPTARG}" ;;          # Configure adapter trimming and coordinates range for B18, CH12, Ramos1, Ramos2, Ramos3
    W) user="${OPTARG}" ;;          # CBE user name as shown when monitoring jobs with squeue
    x) ovlpmm="${OPTARG}" ;;
    X) adapterr="${OPTARG}" ;;
    z) renamefq=1 ;;               # Files are uncompressed .fq instead of compressed .fastq.gz
    1) dopre=0;;                   # skip all steps before alignment, up to and including trimming
    2) doal=0;;                    # Skip alignment, including read pair merge and filtering
    3) dopile=0;;                  # Skip stratification and pileup
    4) dopost=0;;                  # Skip multiqc and collective readcounts
    5) doviz=0;;                   # Skip visualisations
    *) usage ;;
  esac
done

if [ ! -z "$run" ]; then
    data="${data}/${run}"
    process="${process}/${run}"
    results="${results}/${run}"
fi

# Configure adapters for trimming
if [ "$vdj" == 'B18' ] || [ "$vdj" == 'CH12' ]; then
    # First-round PCR primer 3' ends (also trimming the ends of the amplicon)
    five='TCTCCACAGGTGTCCACTCC'     # cutadapt : read1 5'
    three='GTGAGTCCTTACAACCTCTC'    # cutadapt : read1 3', revcomp of $FIVE
    FIVE='GAGAGGTTGTAAGGACTCAC'     # cutadapt : read2 5'
    THREE='GGAGTGGACACCTGTGGAGA'    # cutadapt : read2 3', revcomp of $five
    # echo "Mouse VDJ adapters"
elif [ "$vdj" == 'B18_1' ]; then    # for round 1, where the amplicons were different
    # Two amplicons in the same sequencing pool/barcode, so two primers to trim
    five='CTTTCTCTCCACAGGTGTCCACTCC -g GGATTGATCCTAATAGTGGTGGTAC'     # cutadapt : read1 5'
    three='GTTCAAGAGCAAGGCCACACTGAC -a GGTGAGTCCTTACAACCTCTCTCTT'    # cutadapt : read1 3', revcomp of $FIVE
    FIVE='GTCAGTGTGGCCTTGCTCTTGAAC -G AAGAGAGAGGTTGTAAGGACTCACC'     # cutadapt : read2 5'
    THREE='GGAGTGGACACCTGTGGAGAGAAAG -A GTACCACCACTATTAGGATCAATCC'    # cutadapt : read2 3', revcomp of $five
elif [ "$vdj" == 'CH12_1' ]; then    # for round 1, where the amplicons were different
    # Two amplicons in the same sequencing pool/barcode, so two primers to trim
    five='CTTTCTCTCCACAGGTGTCCACTCC -g ATCCTAGCAATGGTGGTACTAACTAC'     # cutadapt : read1 5'
    three='GTTCAAGAGCAAGGCCACACTGAC -a GGTGAGTCCTTACAACCTCTCTCTT'    # cutadapt : read1 3', revcomp of $FIVE
    FIVE='GTCAGTGTGGCCTTGCTCTTGAAC -G AAGAGAGAGGTTGTAAGGACTCACC'     # cutadapt : read2 5'
    THREE='GGAGTGGACACCTGTGGAGAGAAAG -A GTAGTTAGTACCACCATTGCTAGGAT'    # cutadapt : read2 3', revcomp of $five
elif [ "$vdj" == 'Ramos1' ] || [ "$vdj" == 'Ramos2' ] || [ "$vdj" == 'Ramos3' ] || [ "$vdj" == 'Ramos' ]; then
    # # Middle of first-round PCR primer (preserve the amplicon-specific 3' ends)
    five='TCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
      FIVE='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
    three='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    THREE='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGA'
    # # Just the shared middle part
    # five='GCTCTTCCGATCT'
    # FIVE='GCTCTTCCGATCT'
    # three='AGATCGGAAGAGC'
    # THREE='AGATCGGAAGAGC'
    # echo "Human VDJ adapters"
else
    echo "$vdj is not among the pre-defined configurations for adapter trimming."
    exit 1
fi


wait_for_jobs(){
  echo "waiting"
  sleep 60  # seconds, give time to the scheduler to put up the task
  sleeptime=120  # ask every 2 mins, for the first 10 mins
  n=1
  while true; do
    if [ $(squeue | grep "$user" | grep -c $1) -eq 0 ]; then
      break
    else
      echo "sleep another" $((sleeptime / 60)) "minutes..."
      sleep $sleeptime
    fi
    n=$((n + 1))
    if [ $n -eq 5 ]; then
      sleeptime=300  # if still running after 10 mins, ask every 5 mins
    fi
    if [ $n -eq 10 ]; then
      sleeptime=600  # if still running after 30 mins, ask every 10 mins
    fi
  done
}


prefix="$(echo $base | tr '/' '_')_$(echo $data | tr '/' '_')"

echo ""
if [[ ! -d "${base}/${process}" ]]; then
	mkdir -p "${base}/${process}"
  echo "Created destination: ${base}/${process}"
fi
if [[ ! -d "${base}/${results}" ]]; then
	mkdir -p "${base}/${results}"
  echo "Created destination: ${base}/${results}"
fi
mkdir -p ./mutpe_logs

# Read count tracking
mkdir -p ${base}/${process}/readcounts
metrics="${base}/${process}/readcounts/metrics.txt"


if [[ "$dopre" -eq 1 ]]; then
    if [[ $renamefq ]]; then
        echo ""
        echo "Renaming & compressing ${data}/*.fq to *.fastq.gz"
    	${XDIR}/fileutilities.py T ${base}/${data} --dir fq | ${XDIR}/fileutilities.py P --loop mv {abs} {dir}/{bas}.fastq \; srun gzip {dir}/{bas}.fastq
    fi
    # if [[ $issra ]]; then
    #     echo ""
    #     echo "Converting ${data}/*.sra to *.fastq.gz"
    # 	${XDIR}/fileutilities.py T ${base}/${data} --dir .sra | ${XDIR}/fileutilities.py P --loop srun ${XDIR}/sra2fastq.sh ${base}/${data} ${base}/${data}/{val}
    # fi


    if [[ ! $renamefq ]] && [[ ! $issra ]]; then
        echo ""
        echo "Extracting ${data}/*.bam to *_1/2.fastq.gz"
    	${XDIR}/fileutilities.py T ${base}/${data} --dir .bam | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_fastq.out ,-e ./mutpe_logs/slurm_{ali}_fastq.err ,--time=0-00:10:00 ,-J bam2fq ,--wrap "'samtools fastq -c 1 -1 ${base}/${data}/{ali}_1.fastq.gz -2 ${base}/${data}/{ali}_2.fastq.gz {abs}'"
    	wait_for_jobs bam2fq
    fi

    echo ""
    echo "Running FastQC for *_1/2.fastq.gz"
    mkdir -p ${base}/${process}/fastqc_raw
    ${XDIR}/fileutilities.py T ${base}/${data}/*_1.fastq.gz ${base}/${data}/*_2.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fastqc ,--wrap "'fastqc -o ${base}/${process}/fastqc_raw {abs}'"
    wait_for_jobs fastqc


    ${XDIR}/fileutilities.py T ${base}/${data}/*fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fgzcnt ,--wrap "'${XDIR}/countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs fgzcnt


    echo ""
    echo "Trimming *_1/2.fastq.gz"
    mkdir -p ${base}/${process}/trim
    ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_trimmomatic.out ,-e ./mutpe_logs/slurm_{ali}_trimmomatic.err ,--time=0-00:10:00 ,-J trimomma ,-c 4 ,--wrap "'trimmomatic PE -threads 4 -phred33 -summary ${base}/${process}/trim/{ali}.log ${base}/${data}/{ali}_1.fastq.gz ${base}/${data}/{ali}_2.fastq.gz ${base}/${process}/trim/{ali}_1.trim3.fastq /dev/null ${base}/${process}/trim/{ali}_2.trim3.fastq /dev/null SLIDINGWINDOW:25:25 TRAILING:25 MINLEN:$(( minL1 < minL2 ? minL1 : minL2 ))'"
    wait_for_jobs trimmoma
    ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_cutadapt.out ,-e ./mutpe_logs/slurm_{ali}_cutadapt.err ,--time=0-00:10:00 ,-J cutadapt ,-c 4 ,--wrap "'cutadapt --cores 4 -q 25 -g $five -a $three -G $FIVE -A $THREE -e $adapterr -m ${minL1}:${minL2} -o ${base}/${process}/trim/{ali}_1.trim53.fastq.gz -p ${base}/${process}/trim/{ali}_2.trim53.fastq.gz ${base}/${process}/trim/{ali}_1.trim3.fastq ${base}/${process}/trim/{ali}_2.trim3.fastq > ${base}/${process}/trim/{ali}.log'"
    wait_for_jobs cutadapt


    ${XDIR}/fileutilities.py T ${base}/${process}/trim/*trim3.fastq --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fqcnt ,--wrap "'${XDIR}/countfq.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs fqcnt
    ${XDIR}/fileutilities.py T ${base}/${process}/trim/*trim53.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,-J fgzcnt ,--wrap "'${XDIR}/countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs fgzcnt

    rm ${base}/${process}/trim/*.trim3.fastq


    echo ""
    echo "Running FastQC for *_1/2.trim53.fastq.gz"
    mkdir -p ${base}/${process}/fastqc_posttrim
    ${XDIR}/fileutilities.py T ${base}/${process}/trim/*trim53.fastq.gz --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fastqc ,--wrap "'fastqc -o ${base}/${process}/fastqc_posttrim {abs}'"
    wait_for_jobs fastqc
fi

if [[ "$doal" -eq 1 ]]; then
    if [[ "$flash" -eq 1 ]]; then
        echo ""
        echo "Merging read pairs into ${process}/*.trim53.extendedFrags.fastq.gz"
        mkdir -p ${base}/${process}/flash
    	if [[ "$allowoutties" -eq 0 ]]; then
        	${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_flash.out ,-e ./mutpe_logs/slurm_{ali}_flash.err ,--time=0-00:15:00 ,-J flash ,--wrap "'flash -t 1 -M $maxOver -m $minOver -x $ovlpmm -d ${base}/${process}/flash -o {ali} -z ${base}/${process}/trim/{ali}_1.trim53.fastq.gz ${base}/${process}/trim/{ali}_2.trim53.fastq.gz 2>&1 > ${base}/${process}/flash/{ali}.log'"
        else
            ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_flash.out ,-e ./mutpe_logs/slurm_{ali}_flash.err ,--time=0-00:15:00 ,-J flash ,--wrap "'flash -O -t 1 -M $maxOver -m $minOver -x $ovlpmm -d ${base}/${process}/flash -o {ali} -z ${base}/${process}/trim/{ali}_1.trim53.fastq.gz ${base}/${process}/trim/{ali}_2.trim53.fastq.gz 2>&1 > ${base}/${process}/flash/{ali}.log'"
        fi
        wait_for_jobs flash
    	# I have no use for the uncombined and they get caught up in *_1/2.fastq.gz
    	rm ${base}/${process}/flash/*notCombined_1.fastq.gz ${base}/${process}/flash/*notCombined_2.fastq.gz


        ${XDIR}/fileutilities.py T ${base}/${process}/flash --dir 'extendedFrags.fastq.gz$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fgzcnt ,--wrap "'${XDIR}/countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
        wait_for_jobs fgzcnt


        echo ""
        echo "Running FastQC for *.extendedFrags.fastq.gz"
        mkdir -p ${base}/${process}/fastqc_postmerge
        ${XDIR}/fileutilities.py T ${base}/${process}/flash --dir 'extendedFrags.fastq.gz$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00  ,-J fastqc ,--wrap "'fastqc -o ${base}/${process}/fastqc_postmerge {abs}'"
        wait_for_jobs fastqc

        echo ""
        echo "Filtering out unacceptable merged lengths from *.extendedFrags.fastq.gz"
        ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_flashcutadapt.out ,-e ./mutpe_logs/slurm_{ali}_flashcutadapt.err ,--time=0-00:10:00 ,-J cutadapt ,-c 4 ,--wrap "'cutadapt --cores 4 -m $((refLen - mergelen)) -M $((refLen + mergelen)) -o ${base}/${process}/flash/{ali}.extendedFrags.fltr.fastq.gz ${base}/${process}/flash/{ali}.extendedFrags.fastq.gz > ${base}/${process}/trim/{ali}_mrg.log'"
        wait_for_jobs cutadapt


        ${XDIR}/fileutilities.py T ${base}/${process}/flash --dir 'extendedFrags.fastq.gz$' | fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J fgzcnt ,--wrap "'${XDIR}/countfqgz.sh {val} {abs} ${metrics}_{bas}.txt'"
        wait_for_jobs fgzcnt


        echo ""
        echo "Aligning *.extendedFrags.fltr.fastq.gz single-end to ${aux}/${bowtie2idx}"
    	mkdir -p ${base}/${process}/bowtie2
    	${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_bowtie2.out ,-e ./mutpe_logs/slurm_{ali}_bowtie2.err ,--time=0-01:00:00 ,-J bowtie2 ,-c 4 ,--wrap "'bowtie2 --no-unal -x ${base}/${aux}/${bowtie2idx} -U ${base}/${process}/flash/{ali}.extendedFrags.fltr.fastq.gz -p 2 --very-sensitive-local 2> ${base}/${process}/bowtie2/{ali}.trim53.log | samtools sort -n -@ 2 -O BAM -o ${base}/${process}/bowtie2/{ali}.aln.bam'"
    	wait_for_jobs bowtie2
    else
        echo ""
        mkdir -p ${base}/${process}/bowtie2
        if [[ $unpaired -eq 0 ]]; then
            # Paired-end
            echo "Aligning *.trim53.fastq.gz paired-end to ${aux}/${bowtie2idx}"
        	${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_bowtie2.out ,-e ./mutpe_logs/slurm_{ali}_bowtie2.err ,-J bowtie2 ,-c 4 ,--time=0-01:00:00 ,--wrap "'bowtie2 -p 2 --very-sensitive-local --no-unal --no-discordant --no-contain --no-mixed -x ${base}/${aux}/${bowtie2idx} -1 ${base}/${process}/trim/{ali}_1.trim53.fastq.gz -2 ${base}/${process}/trim/{ali}_2.trim53.fastq.gz 2> ${base}/${process}/bowtie2/{ali}.trim53.log | samtools sort -n -@ 2 -O BAM -o ${base}/${process}/bowtie2/{ali}.aln.bam'"

            # This seems to remove the non-overlapping part of read2 instead of the overlap between read1 and read2. Fail.
            # echo ""
            # echo "Clipping pair overlaps in *.alnovlp.bam"
            # ${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J clipOver ,-c 4 ,--wrap "'bam clipOverlap --readName --stats --in ${base}/${process}/{ali}.alnovlp.bam --out ${base}/${process}/{ali}.aln.bam'"
            # wait_for_jobs clipOver
            # rm ${base}/${process}/*.alnovlp.bam
        else
            # OR Single-end
            echo "Aligning *.trim53.fastq.gz single-end to ${aux}/${bowtie2idx}"
        	${XDIR}/fileutilities.py T ${base}/${data} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_bowtie2.out ,-e ./mutpe_logs/slurm_{ali}_bowtie2.err ,-J bowtie2 ,-c 4 ,--time=0-01:00:00 ,--wrap "'bowtie2 -p 2 --very-sensitive-local --no-unal -x ${base}/${aux}/${bowtie2idx} -U ${base}/${process}/trim/{ali}_1.trim53.fastq.gz,${base}/${process}/trim/{ali}_2.trim53.fastq.gz 2> ${base}/${process}/bowtie2/{ali}.trim53.log | samtools sort -n -@ 2 -O BAM -o ${base}/${process}/bowtie2/{ali}.aln.bam'"
        fi
    	wait_for_jobs bowtie2
    fi


    ${XDIR}/fileutilities.py T ${base}/${process}/bowtie2/*.aln.bam --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-00:10:00 ,-J bamcnt ,--wrap "'${XDIR}/countbam.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs bamcnt
fi


if [[ "$dopile" -eq 1 ]]; then
    echo ""
    echo "Stratifying *.aln.bam by number of mismatches"
    mkdir -p ${base}/${process}/strata
    # Unstratified
    ln -s $(realpath ${base}/${process}/bowtie2/*aln.bam) ${base}/${process}/strata/

    # # REMEMBER to UPDATE the call to PLOTMETRICS.R , if changing the strata
    # # REMEMBER to UPDATE the SUBSTRACTIONS command below, if changing the strata
    # # REMEMBER to UPDATE the PLOT LAYOUT commands below, if changing the strata

    # # Original Yeap et al stratification.
    # ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J stratBAM ,--time=0-0:15:00 ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 1 2'"
    # ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J stratBAM ,--time=0-0:15:00 ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 3 10'"
    # ${XDIR}/fileutilities.py T ${base}/${process} --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J stratBAM ,--time=0-0:15:00 ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh {abs} 11 30'"

    # # More detailed stratification.
    # # Fewer reads have many mutations, so higher strata still need to be pooled.
    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_strat1.out ,-e ./mutpe_logs/slurm_{ali}_strat1.err ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh ${XDIR} {abs} 1 1'"
    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_strat2.out ,-e ./mutpe_logs/slurm_{ali}_strat2.err ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh ${XDIR} {abs} 2 2'"
    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_strat3.out ,-e ./mutpe_logs/slurm_{ali}_strat3.err ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh ${XDIR} {abs} 3 3'"
    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_strat4.out ,-e ./mutpe_logs/slurm_{ali}_strat4.err ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh ${XDIR} {abs} 4 4'"
    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_strat5.out ,-e ./mutpe_logs/slurm_{ali}_strat5.err ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh ${XDIR} {abs} 5 15'"
    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_strat16.out ,-e ./mutpe_logs/slurm_{ali}_strat16.err ,--time=0-00:15:00 ,-J stratBAM ,--wrap "'${XDIR}/split_BAM_by_numMismatches.sh ${XDIR} {abs} 16 30'"

    wait_for_jobs stratBAM


    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln(\.\d+-\d+)?.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,--time=0-0:15:00 ,-J bamcnt ,--wrap "'${XDIR}/countbam.sh {val} {abs} ${metrics}_{bas}.txt'"
    wait_for_jobs bamcnt


    echo ""
    echo "Piling up *\d+-\d+.bam onto ${aux}/${bowtie2idx}.${fastasuffix}"
    mkdir -p ${base}/${process}/pileup
    # Apply map and base quality thresholds to reduce noise.
    # Do not penalize alignment qualities for high mismatches.
    # Do not allow unpaired reads (for unmerged use-case).
    if [ ! -e "${base}/${aux}/${bowtie2idx}.${fastasuffix}" ]; then
        # This has caught me out several times already, so check explicitly
        echo "${base}/${aux}/${bowtie2idx}.${fastasuffix} not found! Did you forget to link it there?"
        exit 1
    fi
    # mpileup wants them sorted, it doesn't know it's an amplicon. But I don't need both files, and I don't want the .sorted suffix in everything downstream.
    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln(\.\d+-\d+)?.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,-o ./mutpe_logs/slurm_{ali}_sort.out ,-e ./mutpe_logs/slurm_{ali}_sort.err ,--time=0-00:15:00 ,-J samsort ,--wrap "'samtools sort {abs} > {dir}/{bas}.sorted.bam'"
    wait_for_jobs samsort
    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln(\.\d+-\d+)?.sorted.bam$' | perl -e 'while(<>){~s/.sorted$//;print}' | ${XDIR}/fileutilities.py P --loop mv {abs} {dir}/{ali}.bam
    ${XDIR}/fileutilities.py T ${base}/${process}/strata --dir 'aln(\.\d+-\d+)?.bam$' | ${XDIR}/fileutilities.py P --loop sbatch ,--mem=10G ,-o ./mutpe_logs/slurm_{ali}_mpileup.out ,-e ./mutpe_logs/slurm_{ali}_mpileup.err ,--time=0-00:15:00 ,-J sampile ,--wrap "'samtools mpileup -f ${base}/${aux}/${bowtie2idx}.${fastasuffix} -d 5000000 -q $minalq -Q $minbasq {abs} > ${base}/${process}/pileup/{bas}.pileup'"
    wait_for_jobs sampile

    echo ""
    echo "Counting *\d+-\d+.pileup"
    ${XDIR}/fileutilities.py T ${base}/${process}/pileup --dir '\.pileup$' | ${XDIR}/fileutilities.py P --loop sbatch ,--time=0-00:10:00 ,-o /dev/null ,-e /dev/null ,-J summpil ,--wrap "'${XDIR}/summarize_mpileup.py -p {abs} -m {dir}/{bas}.point.stats -i {dir}/{bas}.indel.stats'"
    wait_for_jobs summpil
fi

if [[ "$dopost" -eq 1 ]]; then
    echo ""
    echo "Compiling readcounts from $metrics"
    # Combine and plot all metrics files
    printf '%s\t%s\n' "file" "count" > $metrics
    cat ${metrics}_* >> $metrics
    #rm ${metrics}_*
    head -n 1 $metrics > ${base}/${process}/tmp
    tail -n +2 $metrics | sort >> ${base}/${process}/tmp
    mv ${base}/${process}/tmp $metrics
    # UPDATE strata HERE
    sbatch --time=0-0:05:00 -o /dev/null -e /dev/null ${XDIR}/plotmetrics.R $metrics ${base}/${results}/${run/\//_}_readcounds.pdf "1-1" "2-2" "3-3" "4-4" "5-15" "16-30"

    echo ""
    echo "Collecting MultiQC for ${run}"
    mkdir -p ${base}/${process}/multiqc_pre
    mkdir -p ${base}/${process}/multiqc_posttrim
    sbatch -o /dev/null -e /dev/null multiqc -f -o ${base}/${process}/multiqc_pre ${base}/${process}/fastqc_raw
    if [[ $flash -eq 1 ]]; then
        mkdir -p ${base}/${process}/multiqc_postmerge
        sbatch --time=0-0:05:00 -o /dev/null -e /dev/null multiqc -f -o ${base}/${process}/multiqc_posttrim ${base}/${process}/fastqc_posttrim
        sbatch --time=0-0:05:00 -o /dev/null -e /dev/null multiqc -f -o ${base}/${process}/multiqc_postmerge ${base}/${process}/fastqc_postmerge ${base}/${process}/bowtie2
    else
        sbatch --time=0-0:05:00 -o /dev/null -e /dev/null multiqc -f -o ${base}/${process}/multiqc_posttrim ${base}/${process}/fastqc_posttrim ${base}/${process}/bowtie2
    fi
fi


if [[ "$doviz" -eq 1 ]]; then
    echo ""
    # UPDATE strata HERE
    # In retrospect using '-' for the null offset was a bad idea. With both quotation symbols already used to prevent from being seen as a flag, I can't nest it further for sbatch.
    # echo "Postprocessing ${process} mutation frequencies"
    sbatch --time=0-0:05:00 --mem=1G -J stats2rds -o /dev/null -e /dev/null ${XDIR}/stats2rds.sh ${XDIR} ${base}/${process}/pileup $vdj $offsets $allelecutoff
    wait_for_jobs stats2rds
    sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/pileup/*aln.point.stats.RDS
    sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/pileup/*aln.1-1.point.stats.RDS
    sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/pileup/*aln.2-2.point.stats.RDS
    sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/pileup/*aln.3-3.point.stats.RDS
    sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/pileup/*aln.4-4.point.stats.RDS
    sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/pileup/*aln.5-15.point.stats.RDS
    sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/pileup/*aln.16-30.point.stats.RDS
    wait_for_jobs rds2cntg

    echo "Plotting ${process} into ${results}"
    # Plot by stratum, with shared Y scale throughout each stratum for easier comparisons.
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/$(basename $process).all.point.Y ${base}/${process}/pileup/*aln.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/$(basename $process).1-1.point.Y ${base}/${process}/pileup/*aln.1-1.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/$(basename $process).2-2.point.Y ${base}/${process}/pileup/*aln.2-2.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/$(basename $process).3-3.point.Y ${base}/${process}/pileup/*aln.3-3.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/$(basename $process).4-4.point.Y ${base}/${process}/pileup/*aln.4-4.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/$(basename $process).5-15.point.Y ${base}/${process}/pileup/*aln.5-15.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/$(basename $process).16-30.point.Y ${base}/${process}/pileup/*aln.16-30.point.stats.RDS
    # Same as before, but with individual Y scales
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/$(basename $process).all.point ${base}/${process}/pileup/*aln.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/$(basename $process).1-1.point ${base}/${process}/pileup/*aln.1-1.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/$(basename $process).2-2.point ${base}/${process}/pileup/*aln.2-2.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/$(basename $process).3-3.point ${base}/${process}/pileup/*aln.3-3.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/$(basename $process).4-4.point ${base}/${process}/pileup/*aln.4-4.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/$(basename $process).5-15.point ${base}/${process}/pileup/*aln.5-15.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/$(basename $process).16-30.point ${base}/${process}/pileup/*aln.16-30.point.stats.RDS
    # Plot by stratum, with individual Y scales and extra stats, but without labels and coverage
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/$(basename $process).all.point ${base}/${process}/pileup/*aln.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/$(basename $process).1-1.point ${base}/${process}/pileup/*aln.1-1.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/$(basename $process).2-2.point ${base}/${process}/pileup/*aln.2-2.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/$(basename $process).3-3.point ${base}/${process}/pileup/*aln.3-3.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/$(basename $process).4-4.point ${base}/${process}/pileup/*aln.4-4.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/$(basename $process).5-15.point ${base}/${process}/pileup/*aln.5-15.point.stats.RDS
    sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/$(basename $process).16-30.point ${base}/${process}/pileup/*aln.16-30.point.stats.RDS
    wait_for_jobs mutpeviz

    echo "BedGraphs for ${results}"
    mkdir -p ${base}/${results}/tracks
    sbatch --time=0-0:01:00 -J mutbed -o /dev/null -e /dev/null ${XDIR}/rds2bedgraph.R ${base}/${results}/tracks ${base}/${process}/pileup/*RDS

    if [ ! -z "$subs" ] ; then
        echo ""
        echo "Creating substractions for ${process}"
        mkdir -p ${base}/${process}/substractions
        mkdir -p ${base}/${results}/substractions
        # UPDATE strata HERE
        # Looping cannot take just '' or '.' target values, they get auto-substituted for current directory, which breaks the unstratified.
        # So I have to keep some more of the filename, even though it's repetitive.
        sbatch --time=0-0:02:00 --mem=1G -J mutpesub -o /dev/null -e /dev/null ${XDIR}/fileutilities.py L $subs --loop ${XDIR}/substract_stats.R {abs}.aln.point.stats {ali}.aln.point.stats ${base}/${process}/substractions/\$\(basename {abs}\)_\$\(basename {ali}\).sub.point.stats
        sbatch --time=0-0:02:00 --mem=1G -J mutpesub -o /dev/null -e /dev/null ${XDIR}/fileutilities.py L $subs --loop ${XDIR}/substract_stats.R {abs}.aln.1-1.point.stats {ali}.aln.1-1.point.stats ${base}/${process}/substractions/\$\(basename {abs}\)_\$\(basename {ali}\).sub.1-1.point.stats
        sbatch --time=0-0:02:00 --mem=1G -J mutpesub -o /dev/null -e /dev/null ${XDIR}/fileutilities.py L $subs --loop ${XDIR}/substract_stats.R {abs}.aln.2-2.point.stats {ali}.aln.2-2.point.stats ${base}/${process}/substractions/\$\(basename {abs}\)_\$\(basename {ali}\).sub.2-2.point.stats
        sbatch --time=0-0:02:00 --mem=1G -J mutpesub -o /dev/null -e /dev/null ${XDIR}/fileutilities.py L $subs --loop ${XDIR}/substract_stats.R {abs}.aln.3-3.point.stats {ali}.aln.3-3.point.stats ${base}/${process}/substractions/\$\(basename {abs}\)_\$\(basename {ali}\).sub.3-3.point.stats
        sbatch --time=0-0:02:00 --mem=1G -J mutpesub -o /dev/null -e /dev/null ${XDIR}/fileutilities.py L $subs --loop ${XDIR}/substract_stats.R {abs}.aln.4-4.point.stats {ali}.aln.4-4.point.stats ${base}/${process}/substractions/\$\(basename {abs}\)_\$\(basename {ali}\).sub.4-4.point.stats
        sbatch --time=0-0:02:00 --mem=1G -J mutpesub -o /dev/null -e /dev/null ${XDIR}/fileutilities.py L $subs --loop ${XDIR}/substract_stats.R {abs}.aln.5-15.point.stats {ali}.aln.5-15.point.stats ${base}/${process}/substractions/\$\(basename {abs}\)_\$\(basename {ali}\).sub.5-15.point.stats
        sbatch --time=0-0:02:00 --mem=1G -J mutpesub -o /dev/null -e /dev/null ${XDIR}/fileutilities.py L $subs --loop ${XDIR}/substract_stats.R {abs}.aln.16-30.point.stats {ali}.aln.16-30.point.stats ${base}/${process}/substractions/\$\(basename {abs}\)_\$\(basename {ali}\).sub.16-30.point.stats
        wait_for_jobs mutpesub

        sbatch --time=0-0:05:00 --mem=1G -J stats2rds -o /dev/null -e /dev/null ${XDIR}/stats2rds.sh ${XDIR} ${base}/${process}/substractions $vdj $offsets $allelecutoff
        wait_for_jobs stats2rds
        sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/substractions/*sub.point.stats.RDS
        sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/substractions/*sub.1-1.point.stats.RDS
        sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/substractions/*sub.2-2.point.stats.RDS
        sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/substractions/*sub.3-3.point.stats.RDS
        sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/substractions/*sub.4-4.point.stats.RDS
        sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/substractions/*sub.5-15.point.stats.RDS
        sbatch --time=0-0:01:00 --mem=1G -J rds2cntg -o /dev/null -e /dev/null ${XDIR}/contingency_tables.R ${base}/${process}/substractions/*sub.16-30.point.stats.RDS
        wait_for_jobs rds2cntg


        echo "Plotting ${process} into ${results}"
        # Plot by stratum, with shared Y scale throughout each stratum for easier comparisons.
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/substractions/$(basename $process).sub.all.point.Y ${base}/${process}/substractions/*sub.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/substractions/$(basename $process).sub.1-1.point.Y ${base}/${process}/substractions/*sub.1-1.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/substractions/$(basename $process).sub.2-2.point.Y ${base}/${process}/substractions/*sub.2-2.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/substractions/$(basename $process).sub.3-3.point.Y ${base}/${process}/substractions/*sub.3-3.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/substractions/$(basename $process).sub.4-4.point.Y ${base}/${process}/substractions/*sub.4-4.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/substractions/$(basename $process).sub.5-15.point.Y ${base}/${process}/substractions/*sub.5-15.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 auto auto ${base}/${results}/substractions/$(basename $process).sub.16-30.point.Y ${base}/${process}/substractions/*sub.16-30.point.stats.RDS
        # Same as before, but with individual Y scales
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/substractions/$(basename $process).sub.all.point ${base}/${process}/substractions/*sub.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/substractions/$(basename $process).sub.1-1.point ${base}/${process}/substractions/*sub.1-1.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/substractions/$(basename $process).sub.2-2.point ${base}/${process}/substractions/*sub.2-2.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/substractions/$(basename $process).sub.3-3.point ${base}/${process}/substractions/*sub.3-3.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/substractions/$(basename $process).sub.4-4.point ${base}/${process}/substractions/*sub.4-4.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/substractions/$(basename $process).sub.5-15.point ${base}/${process}/substractions/*sub.5-15.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_labels_layout.R 0.97 'NULL' 'NULL' ${base}/${results}/substractions/$(basename $process).sub.16-30.point ${base}/${process}/substractions/*sub.16-30.point.stats.RDS
        # Plot by stratum, with individual Y scales and extra stats, but without labels and coverage
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/substractions/$(basename $process).sub.all.point ${base}/${process}/substractions/*sub.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/substractions/$(basename $process).sub.1-1.point ${base}/${process}/substractions/*sub.1-1.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/substractions/$(basename $process).sub.2-2.point ${base}/${process}/substractions/*sub.2-2.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/substractions/$(basename $process).sub.3-3.point ${base}/${process}/substractions/*sub.3-3.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/substractions/$(basename $process).sub.4-4.point ${base}/${process}/substractions/*sub.4-4.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/substractions/$(basename $process).sub.5-15.point ${base}/${process}/substractions/*sub.5-15.point.stats.RDS
        sbatch --time=0-0:02:00 -J mutpeviz -o /dev/null -e /dev/null ${XDIR}/plot_freqs_layout.R 'NULL' ${base}/${results}/substractions/$(basename $process).sub.16-30.point ${base}/${process}/substractions/*sub.16-30.point.stats.RDS
        wait_for_jobs mutpeviz


        echo "BedGraphs for ${results}"
        mkdir -p ${base}/${results}/substractions/tracks
        sbatch --time=0-0:05:00 -J mutbed -o /dev/null -e /dev/null ${XDIR}/rds2bedgraph.R ${base}/${results}/substractions/tracks ${base}/${process}/substractions/*RDS
    fi
fi

# echo ""
# echo "Cleaning up file-pollution"
# Trash (most) intermediate files
# if [[ "$dopre" -eq 1 ]]; then
#     rm ${base}/${process}/*trim53.fastq.gz
# fi
# if [[ "$flash" -eq 1 ]] && [[ "$doaln" -eq 1 ]]; then
#     rm ${base}/${process}/*extendedFrags*.fastq.gz
#     if [ "$allowoutties" -eq 1];  then
#         rm ${base}/${process}/*histogram.outie
#         rm ${base}/${process}/*histogram.innie
#         rm ${base}/${process}/*hist.innie
#         rm ${base}/${process}/*hist.outie
#     else
#         rm ${base}/${process}/*histogram
#         rm ${base}/${process}/*hist
#     fi
# fi
# if [[ "$dopile" -eq 1 ]]; then
#     ${XDIR}/fileutilities.py T ${base}/${process} --dir '\d+\-\d+.bam' | ${XDIR}/fileutilities.py P --loop rm {abs}
#     rm ${base}/${process}/*pileup
# fi


echo ""
echo "Finished ${run}"
