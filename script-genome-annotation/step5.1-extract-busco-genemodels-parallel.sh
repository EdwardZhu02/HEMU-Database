#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -o ./logs/extract-busco-genemodels.out
#SBATCH -e ./logs/extract-busco-genemodels.log 

source activate maker

if [ ! -d "./augustus/round1" ]; then
	mkdir -p ./augustus/round1
fi

cd ./augustus/round1

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' \
	../../snap/round1/Microstegium_rnd1.all.gff | \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
bedtools getfasta \
	-fi ../../microstegium_genome.fa \
	-bed - \
	-fo Microstegium_rnd1.all.maker.transcripts1000.fasta