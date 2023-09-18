#!/bin/bash

# Environment req: /data21/zhuyuzhi/condaenvs/trinity

# Check log folder
if [ ! -d "./logs" ]; then
	mkdir ./logs
fi
mkdir ./repRNAseq

cd /data21/wongzj/Themeda/Themeda_triandra
for accession in ` ls *.fastq.gz | cut -d"_" -f1 | uniq `
do  
	cd /data21/zhuyuzhi/maker_test_Themeda/repRNAseq
	echo "===Read Filtering: ${accession}==="

	ln -s /data21/wongzj/Themeda/Themeda_triandra/${accession}_1.fastq.gz ./${accession}_1.fastq.gz
	ln -s /data21/wongzj/Themeda/Themeda_triandra/${accession}_2.fastq.gz ./${accession}_2.fastq.gz

	fastp \
		-i ./${accession}_1.fastq.gz -I ./${accession}_2.fastq.gz \
		-o ./${accession}_1_Filtered.fastq.gz -O ./${accession}_2_Filtered.fastq.gz \
		-j ./${accession}_Report.json \
		-h ./${accession}_Report.html &
done
wait

cd /data21/wongzj/Themeda/Themeda_triandra
for accession in ` ls *.fastq.gz | cut -d"_" -f1 | uniq `
do
	cd /data21/zhuyuzhi/maker_test_Themeda/repRNAseq
	gunzip ./${accession}_1_Filtered.fastq.gz &
	gunzip ./${accession}_2_Filtered.fastq.gz &
done
wait

cd /data21/zhuyuzhi/maker_test_Themeda/repRNAseq
echo "===Merge Left and Right Reads==="
cat ./*_1_Filtered.fastq > left_all.fq &
cat ./*_2_Filtered.fastq > right_all.fq &
wait

echo "===Delete Temp Reads==="
rm -rf ./*_Filtered.fastq
wait
