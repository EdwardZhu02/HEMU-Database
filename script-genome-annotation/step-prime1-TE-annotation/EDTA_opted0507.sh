#!/bin/bash

project_path=$(cd `dirname $0`; pwd)
project_name="${project_path##*/}"

for i in $(ls *.fna);
do
  cat ${i}
  # Output a blank line after each fna
  echo
done > ./${project_name}.original.fa

#cat ./${project_name}.fa | grep ">"
# rm ./*.fna
# rm ./*.jsonl # Delete Sequence Report
wait

perl /data5/Andropogoneae_LTR_Project/AUX_SOFTWARES/translate_fasta_headers_perl/translate_fasta_headers.pl \
  --out=./${project_name}.fa \
  ./${project_name}.original.fa # Input for long-seq name FA

wait

mkdir ./${project_name}_FA_sequence_id_translation
mv ./${project_name}.fa.translation.tab ./${project_name}_FA_sequence_id_translation/
rm ./${project_name}.original.fa
wait

perl /data5/Andropogoneae_LTR_Project/AUX_SOFTWARES/EDTA/EDTA.pl \
  --genome ./${project_name}.fa \
  --step all \
  --overwrite 1 \
  --sensitive 0 \
  --anno 1 \
  --threads 8 \

# Do not use RepeatModeler to identify remaining TEs. May cost less time.
# Do not pick up previous results