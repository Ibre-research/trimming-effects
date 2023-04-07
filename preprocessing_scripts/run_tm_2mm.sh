#!/bin/bash


RAW_FQ_DIR='/media/array/callers_proj/raw_data/fastqs'
OUR_FQ_DIR='/media/array/callers_proj/preprocessing_comp/fastqs'
OUR_LOGS_DIR='/media/array/callers_proj/preprocessing_comp/logs'
TM_DIR='/home/barbitoff/anaconda3/envs/preproc/share/trimmomatic'

for FQ in ${RAW_FQ_DIR}/HG001*R1.fastq.gz
do
	TMPTAG=${FQ%%.R1.fastq.gz}
	TAG=${TMPTAG##*fastqs/}
	SMP=$( echo $TAG | grep -oP '.*OME' )
	echo "Aligning $TAG"

	while [ $( ps -AF | grep trimm | wc -l ) -ge 12 ] ; do sleep 1 ; done
	
	java -jar ${TM_DIR}/trimmomatic.jar PE $FQ ${TMPTAG}.R2.fastq.gz ${OUR_FQ_DIR}/${TAG}.paired.R1.fastq.gz ${OUR_FQ_DIR}/${TAG}.unpaired.R1.fastq.gz ${OUR_FQ_DIR}/${TAG}.paired.R3.fastq.gz ${OUR_FQ_DIR}/${TAG}.unpaired.R2.fastq.gz ILLUMINACLIP:${TM_DIR}/adapters/TruSeq3-PE.fa:2:30:10:2:True &> ${OUR_LOGS_DIR}/${TAG}.TM.log &

done
