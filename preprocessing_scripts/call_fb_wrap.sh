#!/bin/bash

#sleep 18000

for i in *{NOTRIM,FASTP}*.dedup.bam
do
	sleep 2
	while [ $( ps -AF | grep 'call_fb.sh' | wc -l ) -ge 8 ] ; do sleep 1 ; done
	/mnt/data/callers_proj/preprocessing_comp/call_fb.sh $i &
done
