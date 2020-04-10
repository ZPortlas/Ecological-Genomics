#!/bin/bash

bismark --bowtie2 --multicore 1 \
    --genome /data/project_data/epigenetics/reference_genome \
    --output_dir /data/project_data/epigenetics/aligned_output \
    -1 /data/project_data/epigenetics/trimmed_fastq/AH_F25_3_1.fq.gz \
    -2 /data/project_data/epigenetics/trimmed_fastq/AH_F25_3_2.fq.gz \
    --rg_tag --rg_id AH_F25_3 --rg_sample AH_F25_3 --gzip --local --maxins 1000