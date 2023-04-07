#!/bin/bash

PROJECT_DIR='/mnt/data/callers_proj/'
REFERENCE='GRCh37.primary_assembly.genome.fa'

for j in NOTRIM FASTP
do
	for BAM in ${PROJECT_DIR}/preprocessing_comp/bams/HG00*_LOWEXOME_${j}*.dedup.bam
	do
		i=$( basename $BAM )
		docker run -v ${PROJECT_DIR}:/data google/deepvariant:latest /opt/deepvariant/bin/run_deepvariant --model_type=WES --ref=/data/raw_data/$REFERENCE --reads=/data/preprocessing_comp/bams/$i --regions=/data/raw_data/gencode.v19.revamp.bed --output_vcf=/data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_DV_STANDART.vcf --output_gvcf=/data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.DV.g.vcf --num_shards=3 &
	done
	wait
done
