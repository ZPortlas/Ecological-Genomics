#!/bin/bash

cd /data/project_data/RS_RNASeq/fastqc/cleanreads/

for file in NOR*C*.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC27_index -l A -r ${file} --validateMappings -o /data/project_data/RS_RNASeq/salmon/cleanedreads/${file}

done

for file in NOR*D*.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC27_index -l A -r ${file} --validateMappings -o /data/project_data/RS_RNASeq/salmon/cleanedreads/${file}

done