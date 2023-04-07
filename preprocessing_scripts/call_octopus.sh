#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'

i=$1

docker run -v /mnt/data/callers_proj/:/data dancooke/octopus octopus -R /data/raw_data/$REFERENCE -I /data/preprocessing_comp/bams/$i -t /data/raw_data/gencode.v19.revamp.octopus.list -o /data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_OCTOPUS_STANDART.vcf --threads 4

docker run -v /mnt/data/callers_proj/:/data dancooke/octopus octopus -R /data/raw_data/$REFERENCE -I /data/preprocessing_comp/bams/$i --filter-vcf /data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_OCTOPUS_STANDART.vcf --forest /opt/octopus/resources/forests/germline.v0.7.4.forest -o /data/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_OCTOPUS_FOREST.vcf --threads 4
