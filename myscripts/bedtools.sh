#!/bin/bash

/data/popgen/bedtools2/bin/bedtools closest -a ~/Ecological-Genomics/AA_AH_diffmeth.bed \
      -b /data/project_data/epigenetics/GCA_900241095.1_Aton1.0_genomic.fa_annotation_table.bed \
      -D b | \
      awk '!($10=="-")' > AA_AH_hits.bed
	  
	  
	  
/data/popgen/bedtools2/bin/bedtools closest -a ~/Ecological-Genomics/AA_HH_diffmeth.bed \
      -b /data/project_data/epigenetics/GCA_900241095.1_Aton1.0_genomic.fa_annotation_table.bed \
      -D b | \
      awk '!($10=="-")' > AA_HH_hits.bed
	  