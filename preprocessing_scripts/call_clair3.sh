#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'
CALLINGREGIONS='gencode.v19.revamp.bed.gz'
PREFIX='/mnt/data/callers_proj/'

i=$1

mkdir ${PREFIX}/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_clair3/

#wget http://www.bio8.cs.hku.hk/clair_models/illumina/12345.tar
#tar -xf 12345.tar

singularity run --bind ${PREFIX}:/data ${PREFIX}/preprocessing_comp/clair3.simg  \
	/opt/bin/run_clair3.sh --bam_fn=/data/preprocessing_comp/bams/$i \
	--ref_fn=/data/raw_data/$REFERENCE --bed_fn=/data/raw_data/$CALLINGREGIONS \
	--threads=12 --platform="ilmn" --model_path="/opt/models/ilmn" \
	--output=/data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_clair3/


