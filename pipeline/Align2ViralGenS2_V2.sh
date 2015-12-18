#!/bin/bash

#########################
# In preprocess_pipeline_core.sh, we filtered reads of LC, phiX, hg38 and rRNA
# This script aims to align LEFT unmapped samples to Viral/Bacterial genomes.

# version: V-0.0.1

# 1. Download Viral/Bacterial genomes from NCBI(ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/) to local server

# 2. Use SNAP to build index of reference genomes.

# 3. Use SNAP to map samples to indexes.
#
# This script use SNAP as default aligner. And perform the part 3.
#########################
module use /apps/modulefiles/test
module load snap/1.0.20

usage()
{
	echo "usage: `basename $0` [sample_name] [queue_name]"
}
#######################
#Input
if [[ $# -ne 2 ]]; then
	usage
	exit 1
fi
input_sample_name=$1
queue_name=$2
input_sample_dir=/PHShome/tw786/neurogen/Tao/kraken_output/$input_sample_name
if [[ ! -d $input_sample_dir ]]; then
	echo $input_sample_dir "doesn't exist, exit!"
	exit 1
fi

output_dir=/data/neurogen/Tao/snap_run_output
[[ -d $output_dir ]] || mkdir $output_dir

output_sample_dir=$output_dir/$input_sample_name
[[ -d $output_sample_dir ]] || mkdir $output_sample_dir

input_sample=$input_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap_unmapped.fastq

if [[ ! -f $input_sample ]]; then
	echo ${input_sample} "doesn't exist! Exit"
	exit 1
fi
# sample_file_kb=`du -k "$input_sample" | cut -f1`
# if [[ $sample_file_kb -le 1000000 ]]; then
# 	queue_name=normal
# else
# 	queue_name=normal
# fi


viruses_root_dir=/data/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses
for viruses_dir in $viruses_root_dir/*
do
	if [[ -d $viruses_dir/snap_index ]]; then
		viruses_name=${viruses_dir##*/}
		result_dir=$output_sample_dir/$viruses_name
		[[ -d $result_dir ]] || mkdir $result_dir
		if  test -f $result_dir/snap_output.sam ; then
			# echo $viruses_name finished >> $output_sample_dir/jobs_summary.txt
			continue
						
		else
			bsub -q "$queue_name" -J ${sample_name}"-mapping-"${viruses_name} -oo $result_dir/snap_mapping_jobout snap-aligner single  \
			$viruses_dir/snap_index -mrl 20 -fastq $input_sample -o $result_dir/snap_output.sam -F a  
			if [[ $? -eq 0 ]]; then
				echo successfully submit $viruses_name >> $result_dir/jobs_summary.txt
			fi
			

		fi
	fi	
	# echo $viruses_dir
	# break	
done





