#!/usr/bin/env bash

#-rw-r--r-- 1 barbitoff barbitoff 139M Sep 20 16:51 /mnt/data/callers_proj/giab/giab42/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
#-rw-r--r-- 1 barbitoff barbitoff 169M Dec  7  2020 /mnt/data/callers_proj/giab/giab42/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
#-rw-r--r-- 1 barbitoff barbitoff 159M Dec  7  2020 /mnt/data/callers_proj/giab/giab42/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
#-rw-r--r-- 1 barbitoff barbitoff 161M Dec  7  2020 /mnt/data/callers_proj/giab/giab42/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
#-rw-r--r-- 1 barbitoff barbitoff 152M Sep 20 16:50 /mnt/data/callers_proj/giab/giab42/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
#-rw-r--r-- 1 barbitoff barbitoff 118M Sep 20 16:49 /mnt/data/callers_proj/giab/giab42/HG006_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
#-rw-r--r-- 1 barbitoff barbitoff 119M Sep 20 16:48 /mnt/data/callers_proj/giab/giab42/HG007_GRCh37_1_22_v4.2.1_benchmark.vcf.gz


#MAIN_DIR=$1

function make_comparison()
{
  PARAMS=("$@") 
  CURR_VCF="${PARAMS[0]}"

  TMP=$(basename "${CURR_VCF}")

  SAMPLE=$(cut -d '_' -f 1 <<< "${TMP}")
  TYPE=$(cut -d '_' -f 2 <<< "${TMP}")
  CALLER=$(cut -d '_' -f 4 <<< "${TMP}")
  ALIGNER=$(cut -d '_' -f 3 <<< "${TMP}")

  echo $SAMPLE $TYPE $ALIGNER $CALLER 

  if [[ "${SAMPLE}" == "HG001" ]]; then
      ref_bed=/mnt/data/callers_proj/giab/giab42/master.bed
      ref_vcf=/mnt/data/callers_proj/giab/giab42/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/mnt/data/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/mnt/data/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG002" ]]; then
      ref_bed=/mnt/data/callers_proj/giab/giab42/master.bed
      ref_vcf=/mnt/data/callers_proj/giab/giab42/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/mnt/data/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/mnt/data/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG003" ]]; then
      ref_bed=/mnt/data/callers_proj/giab/giab42/master.bed
      ref_vcf=/mnt/data/callers_proj/giab/giab42/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/mnt/data/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/mnt/data/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG004" ]]; then
      ref_bed=/mnt/data/callers_proj/giab/giab42/master.bed
      ref_vcf=/mnt/data/callers_proj/giab/giab42/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/mnt/data/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/mnt/data/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG005" ]]; then
      ref_bed=/mnt/data/callers_proj/giab/giab42/master.bed
      ref_vcf=/mnt/data/callers_proj/giab/giab42/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/mnt/data/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/mnt/data/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG006" ]]; then
      ref_bed=/mnt/data/callers_proj/giab/giab42/master.bed
      ref_vcf=/mnt/data/callers_proj/giab/giab42/HG006_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/mnt/data/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/mnt/data/callers_proj/bed/exome_sureselect.bed
      fi
  fi
  if [[ "${SAMPLE}" == "HG007" ]]; then
      ref_bed=/mnt/data/callers_proj/giab/giab42/master.bed
      ref_vcf=/mnt/data/callers_proj/giab/giab42/HG007_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
      if [[ "${TYPE}" == "GENOME" ]]; then
        bed=/mnt/data/callers_proj/raw_data/gencode.v19.revamp.bed
      else
        bed=/mnt/data/callers_proj/bed/exome_sureselect.bed
      fi
  fi

  HAPPY=/mnt/data/callers_proj/happy/bin/hap.py
  REF=/mnt/data/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa
  REF_RTG=/mnt/data/callers_proj/GRCh37.primay_assembly.genome.rtg/
  RTG_SH=/mnt/data/callers_proj/rtg-tools-3.12/rtg
  EVAL_DIR=/mnt/data/callers_proj/preprocessing_comp/eval_v42/${TMP%%.vcf.gz}_eval_data/
#  rm -rf $EVAL_DIR
  mkdir $EVAL_DIR

  python $HAPPY $ref_vcf $CURR_VCF -r $REF -f $ref_bed --threads 4 --engine vcfeval -T $bed -o "${EVAL_DIR}/report" --engine-vcfeval-template $REF_RTG --engine-vcfeval-path $RTG_SH --stratification /mnt/data/callers_proj/GRCh37-giab-stratifications/v2.0-GRCh37-stratifications_plusCustom.tsv
  echo $SAMPLE $TYPE $ALIGNER $CALLER 
  echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
}
export -f make_comparison

/usr/bin/parallel -j 7 make_comparison ::: /mnt/data/callers_proj/preprocessing_comp/bgzipped/*FB_STAN*.vcf.gz
wait

#for i in $( tree /mnt/data/callers_proj/eval_stratified | grep -B 1 report.runinf | grep -oP 'HG.*' | sed 's/_eval_data/.vcf.gz/' | perl -pe 's|^|/mnt/data/callers_proj/bgzipped/|' ) ; do make_comparison $i ; done


