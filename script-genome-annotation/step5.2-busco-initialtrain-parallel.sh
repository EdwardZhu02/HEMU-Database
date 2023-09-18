#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o ./logs/busco-initialtrain.out
#SBATCH -e ./logs/busco-initialtrain.log 

source activate buscov5

if [ ! -d "./augustus/round1" ]; then
	echo "Please run step5.1 first."
	exit
fi
cd ./augustus/round1

busco \
	-i Microstegium_rnd1.all.maker.transcripts1000.fasta \
	-o Microstegium_rnd1_maker \
	--lineage_dataset /public1/home/sc80041/zhuyuzhi/databases/busco_db/embryophyta_odb10/ \
	--mode genome \
	--evalue 1e-03 \
	--long --augustus \
	--cpu 40 --offline