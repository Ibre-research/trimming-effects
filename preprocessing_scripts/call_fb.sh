#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'
CALLINGREGIONS='gencode.v19.revamp.bed'

i=$1

# FreeBayes
freebayes --standard-filters -t /mnt/data/callers_proj/raw_data/$CALLINGREGIONS -f /mnt/data/callers_proj/raw_data/$REFERENCE $i > /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.raw.vcf

grep -P '^#' /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.raw.vcf  > header
grep -vP '^#' /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.raw.vcf | sort -k1,1 -k2,2n | cat header - > /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.resort.vcf

docker run -v /mnt/data/callers_proj/:/data broadinstitute/gatk:latest gatk LeftAlignAndTrimVariants \
	-V /data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.resort.vcf \
	-R /data/raw_data/$REFERENCE --split-multi-allelics \
	-O /data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.split.vcf

docker run -v /mnt/data/callers_proj/:/data broadinstitute/gatk:latest gatk VariantFiltration \
	-V /data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.FB.split.vcf -filter "RPR < 1" --filter-name "RPR1" \
	-filter "RPL < 1" --filter-name "RPL1" -filter "SAF < 1" --filter-name "SAF1" \
	-filter "SAR < 1" --filter-name "SAR1" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "QUAL / AO < 10.0" \
	--filter-name "QUAQLbyAO10" -O /data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_FB_STANDART.vcf


