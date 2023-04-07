#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'


for j in NOTRIM
do
	for i in *${j}.dedup.bam
	do
		docker run -v /mnt/data/:/mnt/data google/deepvariant:latest /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/mnt/data/callers_proj/raw_data/$REFERENCE --reads=/mnt/data/callers_proj/preprocessing_comp/bams/WGS/$i --regions=/mnt/data/callers_proj/raw_data/gencode.v19.revamp.bed --output_vcf=/mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_DV_STANDART.vcf --output_gvcf=/mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.DV.g.vcf --num_shards=3 &
	done
	wait
done
