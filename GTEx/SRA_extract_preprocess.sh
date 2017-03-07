#!/bin/bash

#From SRA to bam
if [[ $# -ne 1 ]]; then
	echo $0 body_site
	exit 1
fi
body_site=$1
download_dir=/home/tw83/twang/GTEx/dbgap/$body_site

if [[ ! -d $download_dir ]]; then
	echo $download_dir "is not exist, exit!"
	exit 1
fi
SRA_dir=/home/tw83/twang/GTEx/dbgap/$body_site/dbGaP-6101/sra
if [[ ! -d $SRA_dir ]]; then
	echo $SRA_dir "is not exist, exit!"
	exit 1
fi
work_dir=/home/tw83/twang/GTEx/dbgap/$body_site/dbGaP-6101
out_dir=/home/tw83/twang/GTEx/dbgap/$body_site/dbGaP-6101/bam
[[ -d $out_dir ]] || mkdir $out_dir
for sra_file in $SRA_dir/*;do
	sra_name=${sra_file##*/}
	sra_prefix=${sra_name%%.sra} 
	if [[ $sra_file == *.sra ]]; then
		if [[ ! -d $out_dir/$sra_prefix ]]; then
			mkdir $out_dir/$sra_prefix
			bsub -q short -W 12:00 -M 8000 -R 'rusage[mem=8000]' ./SRA_extract_preprocess_core.sh $sra_file $out_dir/$sra_prefix $sra_prefix $work_dir
		fi
	fi
	break
done