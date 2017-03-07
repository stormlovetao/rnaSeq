#!/bin/bash
# bsub -q big-multi -n 8 -M 40000 -R 'rusage[mem=40000]' 
prog=$0
BCRootdir=/PHShome/tw786/neurogen/Tao/BRAINCODE_output
prog_dir=/PHShome/tw786/MyOwnScript
if [[ ! -d $BCRootdir ]]; then
	echo "Not found:" $BCRootdir
fi
function error_exit
{
	echo "$1" 1>&2
	exit 1
}

for sample_path in $BCRootdir/*
do
	if [[ ! -d $sample_path ]]; then
		continue
	fi

	sample_name=${sample_path##*/}

	if [[ ! -f $sample_path/kraken_output ]]; then
		continue
	fi
	if [[ -f $sample_path/kraken_output_M_S.fa.gz ]]; then
		continue
	fi
	echo start $sample_path
	bsub -q big -M 30000 -R 'rusage[mem=40000]' sh $prog_dir/BC_rerun_blast.sh $sample_path
done
