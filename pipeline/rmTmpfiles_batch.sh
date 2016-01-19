#!/bin/bash

kraken_out_dir=/data/neurogen/Tao/kraken_output



for sample_dir in $kraken_out_dir/*
do
	if [[ -d $sample_dir ]]; then
		sample_name=${sample_dir##*/}
		
		if [[ -f $sample_dir/kraken_output.report ]]; then

			bsub -q normal sh /PHShome/tw786/MyOwnScript/rmTmpfiles.sh $sample_dir
			
		else
			echo $sample_name 'does not have kraken_output.report'
			continue
		fi

	fi

done

