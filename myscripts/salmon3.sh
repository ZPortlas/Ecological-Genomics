#!/bin/bash

cd /data/project_data/RS_RNASeq/fastq/cleanreads/

for file in NOR*C*.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_cds_index -l A -r ${file} --validateMappings --seqBias -o /data/project_data/RS_RNASeq/salmon/allmapping/${file}

done

for file in NOR*D*.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_cds_index -l A -r ${file} --validateMappings --seqBias -o /data/project_data/RS_RNASeq/salmon/allmapping/${file}

done