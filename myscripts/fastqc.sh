#!/bin/bash

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/DG*fastq.gz

do

fastqc ${file} -o ~/Ecological-Genomics/myresults/fastqc

done
