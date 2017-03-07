#!/bin/bash
prog_dir=/home/tw83/Myscript
GTEx_out_dir=/home/tw83/twang/GTEx_output
tissue=$1
tissue_path=$GTEx_out_dir/$tissue
function error_exit
{
	echo "$1" 1>&2
	exit 1
}
if [[ ! -d $tissue_path ]]; then
	error_exit "Not found: " tissue_path
fi

for srr_path in $tissue_path/*
do
	if [[ ! -d $srr_path ]]; then
		continue
	fi
	srr_name=${srr_path##*/}

	if [[ ! -f $srr_path/unmapped_qc.fastq.gz ]]; then
		continue
	fi
	bsub -q short -W 1:00 -o $srr_path/vendor_filter_run.log sh $prog_dir/GTEx_filter_vendorfail.sh $tissue $srr_name
done
