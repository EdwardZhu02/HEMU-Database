#!/bin/bash

#project_path=$(cd `dirname $0`; pwd)
#project_name="${project_path##*/}"

perl /data5/Andropogoneae_LTR_Project/AUX_SOFTWARES/EDTA/EDTA.pl \
  --genome ./Zluxu_genome.fa \
  --step all \
  --overwrite 1 \
  --sensitive 0 \
  --anno 1 \
  --threads 12 \

# Do not use RepeatModeler to identify remaining TEs. May cost less time.
# Do not pick up previous results