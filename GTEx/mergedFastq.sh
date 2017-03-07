#!/bin/bash

unmapped_bam_dir=$1
output_root_dir=$2
body_site=$3
if [[ ! -d $unmapped_bam_dir ]]; then
	echo $unmapped_bam_dir "doesn't exist"
	exit 1
fi
[[ -d $output_root_dir ]] || mkdir $output_root_dir
job_out_dir=$output_root_dir/jobs_out
[[ -d $job_out_dir ]] || mkdir $job_out_dir
error_out_dir=$output_root_dir/errors_out
[[ -d $error_out_dir ]] || mkdir $error_out_dir
for unmapped_bam_file_dir in $unmapped_bam_dir/*; do
	if [[ $unmapped_bam_file_dir != *.unmapped.bam ]]; then
		continue
	fi
	unmapped_bam_file_name=${unmapped_bam_file_dir##*/}
	srr_name=${unmapped_bam_file_name%%.unmapped.bam}
	srr_id=${srr_name##SRR}
	output_sample_dir=$output_root_dir/$srr_name
	if [[ ! -d $output_sample_dir ]]; then
		mkdir $output_sample_dir
		bsub -q normal -N -oo $job_out_dir/${srr_name}.jobout -eo $error_out_dir/${srr_name}.jobout \
		-J "$body_site[$srr_id]" ./mergedFastq_core.sh $output_sample_dir $unmapped_bam_file_dir $srr_name
	fi
	
	#break
done