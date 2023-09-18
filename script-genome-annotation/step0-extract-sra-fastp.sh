#!/bin/bash
#SBATCH --job-name Cflex-rnaseq-step1-zhu
#SBATCH --nodelist node2
#SBATCH --nodes 1
#SBATCH --ntasks 8
#SBATCH -o ./logs/Cflex_rnaseq-step1.out
#SBATCH -e ./logs/Cflex_rnaseq-step1.err

source activate RNA_Pipeline

if [ ! -d "./unfiltered_fastq" ]; then
	mkdir -p ./unfiltered_fastq
fi
if [ ! -d "./filtered_fastq" ]; then
	mkdir -p ./filtered_fastq
fi

fastq-dump --split-3 --gzip -outdir ./unfiltered_fastq/ SRR10829620.sra &
fastq-dump --split-3 --gzip -outdir ./unfiltered_fastq/ SRR2970609.sra &
wait

fastp \
	-i ./unfiltered_fastq/SRR10829620_1.fastq.gz \
	-I ./unfiltered_fastq/SRR10829620_2.fastq.gz \
	-o ./filtered_fastq/SRR10829620_1.fastq.gz \
	-O ./filtered_fastq/SRR10829620_2.fastq.gz \
	-j ./filtered_fastq/SRR10829620_Report.json \
	-h ./filtered_fastq/SRR10829620_Report.html &
fastp \
	-i ./unfiltered_fastq/SRR2970609_1.fastq.gz \
	-I ./unfiltered_fastq/SRR2970609_2.fastq.gz \
	-o ./filtered_fastq/SRR2970609_1.fastq.gz \
	-O ./filtered_fastq/SRR2970609_2.fastq.gz \
	-j ./filtered_fastq/SRR2970609_Report.json \
	-h ./filtered_fastq/SRR2970609_Report.html &
wait
