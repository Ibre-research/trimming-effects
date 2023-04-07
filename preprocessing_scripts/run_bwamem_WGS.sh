#!/bin/bash

REFERENCE='/mnt/data/callers_proj/raw_data/GRCh37.primary_assembly.genome.fa'
DATENOW=$( date )
echo "Alignment started at: $DATENOW"

for FQDIR in fastqs/HG00*_GENOME_FASTP
do
	cd $FQDIR
	TAG=${FQDIR##fastqs/}
	echo "Aligning $TAG"
	for R1 in $( find * | grep -P 'R1.fastq.gz' )
	do
		R2=${R1/.R1.fastq.gz/.R2.fastq.gz}
		NODIR=$( basename $R1 )
		RGID=${NODIR%%.R1*}
		echo $R1 $R2 $RGID
		bwa mem -M -t 24 -R "@RG\tID:$RGID\tSM:S${TAG}\tLB:1\tPL:illumina" \
			$REFERENCE ../../${FQDIR}/$R1 ../../${FQDIR}/$R2 \
			2> ../../logs/${TAG}.bwa.log | samtools view -bS - > ../../bams/WGS/${TAG}_${NODIR%%.R1.fastq.gz}_BWA.bam
	done
	cd ../../
done

DATENOW=$( date )
echo "Alignment finished at: $DATENOW"

