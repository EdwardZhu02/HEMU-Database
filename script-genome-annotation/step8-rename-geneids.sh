#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o ./logs/rename_geneids.out
#SBATCH -e ./logs/rename_geneids.log

source activate maker

# create naming table
maker_map_ids --prefix Hdipl_ --justify 6 ./maker-round2-annot-gff/Hyparrhenia_rnd2.all.gff > ./maker-final-curation/Hyparrhenia_rnd2.all.maker.name.map

# replace names in GFF files
map_gff_ids ./maker-final-curation/Hyparrhenia_rnd2.all.maker.name.map ./maker-final-curation/Hyparrhenia_rnd2.all.gff
map_gff_ids ./maker-final-curation/Hyparrhenia_rnd2.all.maker.name.map ./maker-final-curation/Hyparrhenia_rnd2.noseq.gff

# replace names in FASTA headers
map_fasta_ids ./maker-final-curation/Hyparrhenia_rnd2.all.maker.name.map ./maker-final-curation/Hyparrhenia_rnd2.all.maker.transcripts.fasta
map_fasta_ids ./maker-final-curation/Hyparrhenia_rnd2.all.maker.name.map ./maker-final-curation/Hyparrhenia_rnd2.all.maker.proteins.fasta