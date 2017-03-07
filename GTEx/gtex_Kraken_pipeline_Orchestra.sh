#!/bin/bash

input_sample_dir=/home/tw83/twang/GTEx/dbgap/Amygdala/dbGaP-6101/bam
output_dir=/home/tw83/twang/GTEx_output
n=1
for srr_folder_path in $input_sample_dir/*; do
	if [[ $n -gt 5 ]]; then
		break
	fi
	let n=n+1
	srr_name=${srr_folder_path##*/}
	if [[ $srr_name != SRR* ]]; then
		continue
	fi
	if [[ ! -d $output_dir/$srr_name ]]; then
		bsub -q mcore -W 24:00 -n 8 -M 90000 -R 'rusage[mem=90000]' ./gtex_Kraken_pipeline_core_Orchestra.sh $srr_name
	fi
done
	
