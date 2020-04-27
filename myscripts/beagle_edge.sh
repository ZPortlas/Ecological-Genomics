#!/bin.bash

cd /data/project_data/GroupProjects/GWAS_env/bam/edge/

ANGSD -GL 2 -out edge -nThreads 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -bam edge_bam.list

# Estimate covariance matrix and individual admixture proportions
# python pcangsd.py -beagle edge.beagle.gz -admix -selection -o edge_PC -threads 1