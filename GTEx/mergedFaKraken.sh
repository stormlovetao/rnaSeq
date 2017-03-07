#!/bin/bash


work_dir=$1
body_site=$2

if [[ ! -d $work_dir ]]; then
	echo $work_dir does not exists
	exit 1
fi

cd $work_dir

##############################################
# merge fa.gz
if [[ ! -f ${body_site}_merged_unmapped_qc.fa.gz ]]; then
	echo "${body_site}_merged_unmapped_qc.fa.gz doesn't exist"
	exit 1
	
fi



##############################################
# Kraken
##############################################
# #######################
# #Step: 
kraken_standard_db=/PHShome/tw786/neurogen/Tao/kraken_standard_db
kraken_bin=kraken
if [[ ! -f ${body_site}_merged_kraken_output ]]; then
	$kraken_bin  --threads 8                                     \
	--db  $kraken_standard_db                                     \
	--fasta-input ${body_site}_merged_unmapped_qc.fa.gz            \
	--gzip-compressed												\
	> ${body_site}_merged_kraken_output
fi

if [[ $? -ne 0 ]]; then
	echo "ERROR exists in kraken"
fi
