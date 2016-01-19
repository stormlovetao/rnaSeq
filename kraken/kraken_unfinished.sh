
kraken_out_dir=/data/neurogen/Tao/kraken_output
for sample_dir in $kraken_out_dir/*
do
	if [[ -d $sample_dir ]]; then
		sample_name=${sample_dir##*/}
		if [[ -f $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.fastq ]]; then
			if [[ ! -f $sample_dir/kraken_output.report ]]; then
				bsub -q big-multi -n 8 -M 128000 -R 'rusage[mem=128000]' sh /PHShome/tw786/MyOwnScript/kraken_pipeline_core.sh $sample_name
			fi
		fi
	fi
done