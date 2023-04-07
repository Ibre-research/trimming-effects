#!/bin/bash

for i in HG00{6,7}*FASTP*.dedup.bam
do
#	sleep 2
#	while [ $( docker ps | grep 'clair' | wc -l ) -ge 10 ] ; do sleep 1 ; done
	/mnt/data/callers_proj/preprocessing_comp/call_clair3.sh $i
done
