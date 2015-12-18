#!/bin/bash

#########################
#
#This script aims to align all unmapped samples to Viral/Bacterial genomes.

# version: V-0.0.0

# 1. Download Viral/Bacterial genomes from NCBI(ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/) to local server

# 2. Use Bowtie/SNAP to build index of reference genomes.

# 3. Use Bowtie/SNAP to map samples to indexes.
#
# This script use SNAP as default aligner. And perform the part 2.
#########################

module use /apps/modulefiles/test
module load snap/1.0.20

viruses_root_dir=/data/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses

for viruses_dir in $viruses_root_dir/*

do 
	if test -d $viruses_dir
		then
		viruses_name=${viruses_dir##*/}
			if [[ -f $viruses_dir/${viruses_name}".fna" ]]; then
				bsub -q "short" -J ${viruses_name}":Indexing" -oo $viruses_dir/snap_index_jobout snap-aligner index  $viruses_dir/${viruses_name}".fna" $viruses_dir/snap_index/
			fi
	fi


done
