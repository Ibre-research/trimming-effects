#!/bin/bash

REFERENCE='GRCh37.primary_assembly.genome.fa'

i=$1

#docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ dancooke/octopus octopus -R /mnt/data/callers_proj/raw_data/$REFERENCE -I /mnt/data/callers_proj/preprocessing_comp/bams/WGS/$i -t /mnt/data/callers_proj/raw_data/gencode.v19.revamp.octopus.list -o /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_OCTOPUS_STANDART.vcf --threads 4

docker run -v /mnt/data/callers_proj/:/mnt/data/callers_proj/ dancooke/octopus octopus -R /mnt/data/callers_proj/raw_data/$REFERENCE -I /mnt/data/callers_proj/preprocessing_comp/bams/WGS/$i --filter-vcf /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_OCTOPUS_STANDART.vcf --forest /opt/octopus/resources/forests/germline.v0.7.4.forest -o /mnt/data/callers_proj/preprocessing_comp/${i%%.dedup.bam}/${i%%.dedup.bam}_OCTOPUS_FOREST.vcf --threads 4
