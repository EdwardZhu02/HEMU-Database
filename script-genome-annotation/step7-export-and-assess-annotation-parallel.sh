#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o ./logs/assess-annotation.out
#SBATCH -e ./logs/assess-annotation.log

source activate maker

# Create the final merged gff file.
# The “-n” option would produce a gff file without genome sequences.
gff3_merge -s -d \
	./Chrysopogon_rnd2.maker.output/Chrysopogon_rnd2_master_datastore_index.log \
	>./maker-round2-annot-gff/Chrysopogon_rnd2.all.gff

gff3_merge -s -n -d \
	./Chrysopogon_rnd2.maker.output/Chrysopogon_rnd2_master_datastore_index.log \
	>./maker-round2-annot-gff/Chrysopogon_rnd2.noseq.gff

# 6.1 Evaluation of Gene Models

# Count the number of gene models and the gene lengths after each round.
cat ./maker-round2-annot-gff/Themeda_rnd2.noseq.gff | \
	awk '{ if ($3 == "gene") print $0 }' | \
	awk '{ sum += ($5 - $4) } END { print NR, sum / NR }' \
	> ./maker-round2-annot-gff/genemodel_num_and_avg_genelen.txt

# Visualize the AED distribution.
perl step7norun-visualize-AED-distribution.pl \
	-b 0.025 \
	./maker-round2-annot-gff/Chrysopogon_rnd2.noseq.gff \
	> ./maker-round2-annot-gff/aed_distribution.txt