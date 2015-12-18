#!/bin/bash

#########################
#
#This script aims to align all unmapped samples to Viral/Bacterial genomes.

# version: V-0.0.0

# 1. Download Viral/Bacterial genomes from NCBI(ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/) to local server

# 2. Use Bowtie/SNAP to build index of reference genomes.

# 3. Use Bowtie/SNAP to map samples to indexes.
#
# This script use SNAP as default aligner. And perform the part 3.
#########################

module use /apps/modulefiles/test
module load snap/1.0.20

viruses_root_dir=/data/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses

samples_root_dir=/data/neurogen/rnaseq_PD/run_output #Check!!

output_dir=/data/neurogen/Tao/run_output
[[ -d $output_dir ]] || mkdir $output_dir

for sample_dir in $samples_root_dir/*
do
	#echo $sample_dir
	if [[ -d $sample_dir ]]; then
		sample_name=${sample_dir##*/}
		#echo $sample_name
		if [[ $sample_name == "ILB"* ]]; then
			if [[ ${sample_name:(-6):1} -ne "1" ]]; then
				[[ -d $output_dir/$sample_name ]] || mkdir $output_dir/$sample_name
				if [[ -f $output_dir/$sample_name/sample_finish_mark ]]; then
					continue
				fi

				bam_file_path=$sample_dir/unmapped.bam
				bam_file_kb=`du -k "$bam_file_path" | cut -f1`
				if [[ $bam_file_kb -le 1000000 ]]; then
					queue_name=short
				else
					queue_name=big
				fi

				for viruses_dir in $viruses_root_dir/*
				do
					if [[ -d $viruses_dir/snap_index ]]; then
						viruses_name=${viruses_dir##*/}
						result_dir=$output_dir/$sample_name/$viruses_name
						[[ -d $result_dir ]] || mkdir $result_dir
						if  test -f $result_dir/snap_output.sam ; then
							echo $viruses_name finished >> $output_dir/$sample_name/jobs_summary.txt
							continue
						
						else
							bsub -q "$queue_name" -J ${sample_name}"-mapping-"${viruses_name} -oo $result_dir/snap_mapping_jobout snap-aligner single \
							$viruses_dir/snap_index -bam $sample_dir/unmapped.bam -o $result_dir/snap_output.sam -F a > $result_dir/stdout.txt
							echo submit $viruses_name >> $output_dir/$sample_name/jobs_summary.txt

						fi
					fi
					
				done
				
				echo "Sample finished!" > $output_dir/$sample_name/sample_finish_mark 
			fi
		
		fi

	fi
	
done
