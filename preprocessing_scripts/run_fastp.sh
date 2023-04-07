#!/bin/bash


RAW_FQ_DIR='/mnt/data/callers_proj/preprocessing_comp/fastqs'
OUR_FQ_DIR='/mnt/data/callers_proj/preprocessing_comp/fastqs'
OUR_LOGS_DIR='/mnt/data/callers_proj/preprocessing_comp/logs'

for FQ in ${RAW_FQ_DIR}/HG00*_LOWEXOME_NOTR*R1.fastq.gz
do
	TMPTAG=${FQ%%.R1.fastq.gz}
	TAG=${TMPTAG##*fastqs/}
	SMP=$( echo $TAG | grep -oP '.*OME' )
	echo "Trimming $TAG"

	while [ $( ps -AF | grep fastp | wc -l ) -ge 12 ] ; do sleep 1 ; done
	
	/mnt/data/callers_proj/preprocessing_comp/software/fastp -Q -L -i $FQ -I ${TMPTAG}.R2.fastq.gz -o ${OUR_FQ_DIR}/${TAG}_FASTP.R1.fastq.gz -O ${OUR_FQ_DIR}/${TAG}_FASTP.R2.fastq.gz &> ${OUR_LOGS_DIR}/${TAG}.FASTP.log &

done
