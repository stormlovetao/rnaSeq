#!/bin/bash

if [[ $# -ne 2 ]]; then
	echo "check input"
	exit 1
fi

body_site=$1
srr_name=$2
sra_prefix=$srr_name
#echo $body_site $srr_name
work_dir=/home/tw83/twang/GTEx/dbgap/allLeft/dbGaP-6101
success_dir=/home/tw83/twang/GTEx/dbgap/allLeft/dbGaP-6101/success
[[ -d $success_dir ]] || mkdir $success_dir
cd $work_dir
[[ -d $body_site ]] || mkdir $body_site
sra_file=$work_dir/sra/${srr_name}.sra
if [[ ! -f $sra_file ]]; then
	echo $sra_file 'not exists, exit!'
	exit 1
fi
out_dir=$work_dir/$body_site
#tmp_dir=/n/scratch2/tw83
tmp_dir=$work_dir/tmp
[[ -d $tmp_dir ]] || mkdir $tmp_dir
sam-dump -r -u $sra_file | samtools view -b -S - > $tmp_dir/${sra_prefix}.bam
if [[ $? -ne 0 ]]; then
	echo "ERROR exists in sam-dump commond"
	exit 1
fi

samtools view -f 4 $tmp_dir/${sra_prefix}.bam | samtools view -b -S - > $out_dir/${sra_prefix}.unmapped.bam

if [[ $? -ne 0 ]]; then
	echo "ERROR exists in samtools view -f 4..."
else
	touch $success_dir/$srr_name 
	[[ -f $tmp_dir/${sra_prefix}.bam ]] && rm $tmp_dir/${sra_prefix}.bam
fi

