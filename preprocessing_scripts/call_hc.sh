#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'
CALLINGREGIONS='gencode.v19.revamp.bed'
DBSNP='hg19_v0_Homo_sapiens_assembly19.dbsnp138.hg19.vcf.gz'
MILLS_INDELS='hg19_v0_Mills_and_1000G_gold_standard.indels.b37.hg19.vcf.gz'
GENOMES='hg19_v0_Homo_sapiens_assembly19.known_indels_20120518.hg19.vcf.gz'
KG_HICONF_INDELS='hg19_v0_1000G_phase1.snps.high_confidence.b37.hg19.vcf.gz'
HAPMAP='hg19_v0_hapmap_3.3.b37.hg19.vcf.gz'
KG_SITES='b37_1000G_phase3_v4_20130502.sites.hg19.vcf.gz'

i=$1

docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" BaseRecalibrator -R /mnt/data/callers_proj/raw_data/$REFERENCE -I /mnt/data/callers_proj/preprocessing_comp/bams/$i --known-sites /mnt/data/callers_proj/raw_data/bundle/$DBSNP --known-sites /mnt/data/callers_proj/raw_data/bundle/$MILLS_INDELS --known-sites /mnt/data/callers_proj/raw_data/bundle/$GENOMES -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.recal.table -L /mnt/data/callers_proj/raw_data/$CALLINGREGIONS
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" ApplyBQSR -R /mnt/data/callers_proj/raw_data/$REFERENCE -I /mnt/data/callers_proj/preprocessing_comp/bams/$i -bqsr /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.recal.table -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.recal.bam -L /mnt/data/callers_proj/raw_data/$CALLINGREGIONS
	
# Raw HaplotypeCaller
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g"  HaplotypeCaller -R /mnt/data/callers_proj/raw_data/$REFERENCE -I /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.recal.bam -bamout /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.bam -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.RAW.vcf -L /mnt/data/callers_proj/raw_data/$CALLINGREGIONS
        
# CNN scoring and filtration
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" CNNScoreVariants -R /mnt/data/callers_proj/raw_data/$REFERENCE -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.RAW.vcf -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN1.vcf -tensor-type reference -L /mnt/data/callers_proj/raw_data/$CALLINGREGIONS -ip 20
#docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" CNNScoreVariants -I /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.bam -R /mnt/data/callers_proj/raw_data/$REFERENCE -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN1.vcf -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN2.vcf -tensor-type read_tensor -L /mnt/data/callers_proj/raw_data/$CALLINGREGIONS -ip 20

docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" FilterVariantTranches -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN1.vcf --output /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_HC_1DCNN.vcf \
  --info-key CNN_1D --snp-tranche 99.9 --indel-tranche 99.0 \
  --resource /mnt/data/callers_proj/raw_data/bundle/$MILLS_INDELS \
  --resource /mnt/data/callers_proj/raw_data/bundle/$KG_HICONF_INDELS \
  --resource /mnt/data/callers_proj/raw_data/bundle/$HAPMAP 
#docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" FilterVariantTranches -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.CNN2.vcf --output /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_HC_2DCNN.vcf \
#  --resource /mnt/data/callers_proj/raw_data/bundle/$MILLS_INDELS \
#  --resource /mnt/data/callers_proj/raw_data/bundle/$KG_HICONF_INDELS \
#  --resource /mnt/data/callers_proj/raw_data/bundle/$HAPMAP \
#  --info-key CNN_2D --snp-tranche 99.9 --indel-tranche 99.5
	
# Hard filters
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" SelectVariants -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.RAW.vcf -select-type SNP -select-type MIXED -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.snps.vcf
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" VariantFiltration -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.snps.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" \
  --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.snps.flt.vcf
  	
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" SelectVariants -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.RAW.vcf -select-type INDEL -select-type MIXED -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.indels.vcf
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" VariantFiltration -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.indels.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.indels.flt.vcf
	
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" MergeVcfs -I /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.snps.flt.vcf -I /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.indels.flt.vcf -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp.vcf
	
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" CalculateGenotypePosteriors -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp.vcf -supporting /mnt/data/callers_proj/raw_data/bundle/$KG_SITES \
  -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp2.vcf
	
docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ broadinstitute/gatk:latest gatk --java-options "-Xmx8g" VariantFiltration -V /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp2.vcf --genotype-filter-expression "GQ < 20" --genotype-filter-name "GQ20" -O /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp3.vcf

bcftools filter -s GQ20 -e "FORMAT/FT[*]!=''" /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}.HC.hardfilter_tmp3.vcf > /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_HC_HARDFILTER.vcf
