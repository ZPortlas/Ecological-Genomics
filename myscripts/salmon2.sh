#!/bin/bash

cd /data/project_data/RS_RNASeq/fastq/cleanreads/

for file in NOR*C*.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC_index -l A -r ${file} --validateMappings -seqBias -o /data/project_data/RS_RNASeq/salmon/HCmapping/${file}

done

for file in NOR*D*.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC_index -l A -r ${file} --validateMappings -seqBias -o /data/project_data/RS_RNASeq/salmon/HCmapping/${file}

done