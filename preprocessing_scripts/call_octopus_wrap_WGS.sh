#!/bin/bash

for i in *.dedup.bam
do
	sleep 1
	while [ $( docker ps | grep octopus | wc -l ) -ge 6 ] ; do sleep 1 ; done
	/mnt/data/callers_proj/preprocessing_comp/call_octopus_WGS.sh $i &
done
