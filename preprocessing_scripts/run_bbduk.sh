#!/bin/bash


RAW_FQ_DIR='/media/array/callers_proj/raw_data/fastqs'
OUR_FQ_DIR='/media/array/callers_proj/preprocessing_comp/fastqs'
OUR_LOGS_DIR='/media/array/callers_proj/preprocessing_comp/logs'
ADAPTERS='/home/barbitoff/anaconda3/envs/preproc/share/trimmomatic/adapters/TruSeq3-PE.fa'


for FQ in ${RAW_FQ_DIR}/HG001*R1.fastq.gz
do
	TMPTAG=${FQ%%.R1.fastq.gz}
	TAG=${TMPTAG##*fastqs/}
	SMP=$( echo $TAG | grep -oP '.*OME' )
	echo "Trimming $TAG"

	while [ $( ps -AF | grep trimm | wc -l ) -ge 12 ] ; do sleep 1 ; done
	
	bbduk.sh in1=$FQ in2=${TMPTAG}.R2.fastq.gz out1=${OUR_FQ_DIR}/${TAG}_BBDUK.R1.fastq.gz out2=${OUR_FQ_DIR}/${TAG}_BBDUK.R2.fastq.gz ref=$ADAPTERS ktrim=r k=23 mink=11 hdist=1 tpe tbo &> ${OUR_LOGS_DIR}/${TAG}.BBDUK.log &

done
