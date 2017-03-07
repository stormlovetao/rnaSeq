#!/bin/bash

kraken_out_dir=/data/neurogen/Tao/kraken_output



for sample_dir in $kraken_out_dir/*
do
	if [[ -d $sample_dir ]]; then
		sample_name=${sample_dir##*/}
		if [[ ! -f $sample_dir/kraken_output.report ]]; then
				echo $sample_name 'does not have kraken_output.report'
		fi
		if [[ -f $sample_dir/kraken_output.report ]] && [[ -f $sample_dir/unmapped.fastq ]]; then

			bsub -q normal sh /PHShome/tw786/MyOwnScript/rmTmpfiles.sh $sample_dir
			
		else
			
			continue
			
		fi

	fi

done

