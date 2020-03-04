#!/bin/bash

for file in /data/project_data/RS_RNASeq/fastqc/cleanreads/NOR*C*.cl.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC27_index -l A -r ${file} --validateMappings -o /data/project_data/RS_RNASeq/salmon/cleanedreads

done

for file in /data/project_data/RS_RNASeq/fastqc/cleanreads/NOR*D*.cl.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC27_index -l A -r ${read} --validateMappings -o /data/project_data/RS_RNASeq/salmon/cleanedreads

done