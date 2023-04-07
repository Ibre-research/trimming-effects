#!/bin/bash

OUT_LOGS_DIR="/mnt/data/callers_proj/preprocessing_comp/logs"

for i in `seq 1 7`
do
	for FQ in $( find fastqs/HG00${i}_GENOME_NOTRIM/* | grep -P '_R1.*fastq.gz$' )
	do
		TMPTAG=${FQ%%.fastq.gz}
		TAG=${TMPTAG##*/}
		IN_FQ_ROOT=${FQ%/*}
		OUT_FQ_ROOT=${IN_FQ_ROOT/_NOTRIM/_FASTP}
		SMP=$( echo $FQ | grep -oP '[^\/]+OME' )
		echo "Trimming $TAG for $SMP, R1 $FQ, R2 ${FQ/_R1_/_R2_}; storing in $OUT_FQ_ROOT"

		while [ $( ps -AF | grep fastp | wc -l ) -ge 12 ] ; do sleep 1 ; done
		/mnt/data/callers_proj/preprocessing_comp/software/fastp -Q -L -i $FQ -I ${FQ/_R1_/_R2_} -o ${OUT_FQ_ROOT}/${TAG}_FASTP.R1.fastq.gz -O ${OUT_FQ_ROOT}/${TAG}_FASTP.R2.fastq.gz &> ${OUT_LOGS_DIR}/${TAG}.FASTP.log &
	done
done
