#!/bin/bash

# Environment req: /home/test3/miniconda3/envs/maker

# Check log and snap folder
if [ ! -d "./logs" ]; then
	mkdir ./logs
fi
if [ ! -d "./snap/round1" ]; then
	mkdir -p ./snap/round1
fi

cd ./snap/round1

echo 'Started model extraction.'

gff3_merge -d ../../Microstegium_rnd1.maker.output/Microstegium_rnd1_master_datastore_index.log
# Export 'confident' gene models from MAKER
# Use models with AED>0.25 and Aa_count>=50, helps get rid of junky models.
maker2zff -x 0.25 -l 50 Microstegium_rnd1.all.gff
#rename 's/genome/Microstegium_rnd1.zff.length50_aed0.25/g' *
echo 'Completed exporting gene models.'

# Stat gathering and validation
fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
fathom genome.ann genome.dna -validate > validate.log 2>&1

# Collect training sequences and annotations, plus 1000 surrounding bp for training
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

# Create training params
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
echo 'Completed creating training params.'

# HMM assembly
hmm-assembler.pl Microstegium_rnd1.zff.length50_aed0.25 ./params > Microstegium_rnd1.zff.length50_aed0.25.hmm
echo 'Completed assembling hidden markov model.'