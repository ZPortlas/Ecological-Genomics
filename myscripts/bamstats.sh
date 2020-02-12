#!/bin/bash

# set repo

myrepo="/users/z/p/zportlas/Ecological-Genomics"

mypop="DG"

output="/data/project_data/RS_ExomeSeq/mapping"

echo "Num.reads R1 R2 Paired MateMapped Singletons MateMappedDiffChr" >${myrepo}/myresults/${mypop}.flagstats.txt

for file in ${output}/BWA/${mypop}*sorted.rmdup.bam

do
f=${file/sorted.rmdup.bam/}
name=`basename ${f}`
echo ${name} >> ${myrepo}/myresults/${mypop}.names.txt
samtools flagstat ${file} | awk 'NR>=6&&NR<=12 {print $1}' | column -x
done >> ${myrepo}/myresults/${mypop}.flagstats.txt

# calculate depth of coverage from our bam files

for file in ${output}/BWA/${mypop}*sorted.rmdup.bam

do
samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${myrepo}/myresults/${mypop}.coverage.txt

