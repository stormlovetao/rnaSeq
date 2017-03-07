#!/bin/bash

raw_input_dir=/data/neurogen/rnaseq_PD/run_output/
kraken_output_dir=/data/neurogen/Tao/BRAINCODE_output
for file in $raw_input_dir/*;do
	folder_name=${file##*/}
	if [[ ! -f $file/unmapped.bam ]]; then
		continue
	fi
	if [[ ! -d $kraken_output_dir/$folder_name ]]; then
		bsub -q big -M 120000 -R 'rusage[mem=120000]' ./BC_pipeline_forsample.sh $folder_name
	elif [[ ! -f $kraken_output_dir/$folder_name/kraken_output ]]; then
		bsub -q big -M 120000 -R 'rusage[mem=120000]' ./BC_pipeline_forsample.sh $folder_name
	fi
done
