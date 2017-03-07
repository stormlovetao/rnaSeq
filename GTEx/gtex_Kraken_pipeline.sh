#!/bin/bash
if [[ $# -ne 1 ]]; then
	echo $0 body_site
	exit 1
fi
body_site=$1
input_sample_dir=/PHShome/tw786/neurogen/external_download/gtex/dbGap/2016march/${body_site}_unmapped
output_dir=/PHShome/tw786/neurogen/Tao/GTEx_output/$body_site

if [[ ! -d $input_sample_dir ]]; then
	echo $input_sample_dir "doesn't exist!"
	exit 1
fi

[[ -d $output_dir ]] || mkdir $output_dir

for srr_unmapped_bam_path in $input_sample_dir/*; do
	
	srr_name=${srr_unmapped_bam_path##*/}
	if [[ $srr_name != SRR* ]]; then
		continue
	fi
	srr_name=${srr_name%%.unmapped*}
	if [[ ! -d $output_dir/$srr_name ]]; then
		bsub -q big-multi -n 4 -M 30000 -R 'rusage[mem=30000]' ./gtex_Kraken_pipeline_core.sh $srr_name $input_sample_dir $output_dir
	fi
	break
done
	

