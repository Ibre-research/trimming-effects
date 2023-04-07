#!/bin/bash

for i in *LOWEXOME_{FASTP,NOTRIM}*.dedup.bam
do
	sleep 2
	while [ $( ps -AF | grep 'call_hc.sh' | wc -l ) -ge 12 ] ; do sleep 1 ; done
	/mnt/data/callers_proj/preprocessing_comp/call_hc.sh $i &
done
