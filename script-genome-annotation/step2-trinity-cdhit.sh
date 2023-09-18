#!/bin/bash

# Environment req: /data21/zhuyuzhi/condaenvs/trinity

# Clean existing tmp files
#if [ -d "./trinity_result" ]; then
#	rm -rf ./trinity_result
#fi

# Check log folder
if [ ! -d "./logs" ]; then
	mkdir ./logs
fi

cd ./repRNAseq/
Trinity --seqType fq --no_version_check \
	--left left_all.fq \
	--right right_all.fq \
	--CPU 40 \
	--max_memory 160G \
	--min_contig_length 300 \
	--output ../trinity_result

cd ..
if [ -f "./trinity_result/Trinity.fasta" ]; then
	cd-hit-est \
		-i ./trinity_result/Trinity.fasta \
		-o ./cdhit_rmredun_Trinity.fasta \
		-c 0.98a -d 0 \
		-T 20 -M 40960
else
	echo "Trinity run failed, check specs."
fi
