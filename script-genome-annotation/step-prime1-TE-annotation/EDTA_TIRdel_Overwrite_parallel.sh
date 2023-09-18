#!/bin/bash
#SBATCH --job-name EDTA-parallel-zhu
#SBATCH --nodelist master
#SBATCH --nodes 1
#SBATCH --ntasks 8
#SBATCH -o ./EDTA_parallel.out
#SBATCH -e ./EDTA_parallel.err

source activate EDTA

project_name=GCA_002018215.1

mkdir ./GCA_002018215.1.fa.mod.EDTA.raw/rsync_empty
rsync --delete-before -d ./GCA_002018215.1.fa.mod.EDTA.raw/rsync_empty/ ./GCA_002018215.1.fa.mod.EDTA.raw/TIR/
rm -r ./GCA_002018215.1.fa.mod.EDTA.raw/TIR
rm -r ./GCA_002018215.1.fa.mod.EDTA.raw/rsync_empty

perl /data5/Andropogoneae_LTR_Project/AUX_SOFTWARES/EDTA/EDTA.pl \
  --genome ./GCA_002018215.1.fa \
  --step all \
  --overwrite 1 \
  --sensitive 0 \
  --anno 1 \
  --threads 8 \
