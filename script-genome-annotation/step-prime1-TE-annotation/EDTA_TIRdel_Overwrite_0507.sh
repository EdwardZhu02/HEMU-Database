#!/bin/bash

project_path=$(cd `dirname $0`; pwd)
project_name="${project_path##*/}"

mkdir ./${project_name}.fa.mod.EDTA.raw/rsync_empty
rsync --delete-before -d ./${project_name}.fa.mod.EDTA.raw/rsync_empty/ ./${project_name}.fa.mod.EDTA.raw/TIR/
rm -r ./${project_name}.fa.mod.EDTA.raw/TIR
rm -r ./${project_name}.fa.mod.EDTA.raw/rsync_empty

perl /data5/Andropogoneae_LTR_Project/AUX_SOFTWARES/EDTA/EDTA.pl \
  --genome ./${project_name}.fa \
  --step all \
  --overwrite 1 \
  --sensitive 0 \
  --anno 1 \
  --threads 8 \

# Do not use RepeatModeler to identify remaining TEs. May cost less time.
# Do not pick up previous results