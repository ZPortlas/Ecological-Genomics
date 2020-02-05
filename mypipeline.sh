#!/bin/bash

# We'll use this as a wrapper to run our different mapping scripts

# Path to my repo:

myrepo="/users/s/r/zportlas/Ecological_Genomics/Spring_2020"

# echo ${myrepo}

# My population:

mypop="DG"

# Directory to cleaned and paired reads:

input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"

# Directory in common space where we'll store mapping output for everyone to use

output="/data/project_data/RS_ExomeSeq/mapping"

# Run mapping.sh

source ./mapping.sh

# Run the post-processing steps

source ./process_bam.sh
