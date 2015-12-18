#!/bin/bash

# This scrip should be run on a server with more than 128 GB mem, and 16 or more cpus
# bsub -q big-multi -n 16 -M 128000 -R 'rusage[mem=128000]' ./kraken_pipeline_core.sh sample_name
#######################
#module load
module load snap/1.0.20
module load tophat/2.1.0
module load bowtie2
usage()
{
	echo "usage: `basename $0` [sample_name]"
}
error_exit()
{
	echo $1
	exit 1
}
#######################
#Input
if [[ $# -ne 1 ]]; then
	usage
	exit 1
fi
input_sample_name=$1

#Output
kraken_output_dir=/data/neurogen/Tao/kraken_output
[ -d $kraken_output_dir ] || mkdir $kraken_output_dir

output_sample_dir=$kraken_output_dir/$input_sample_name
[ -d $output_sample_dir ] || mkdir $output_sample_dir

#######################

# Check if input fastq file exists!

if [[ ! -e $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap_unmapped.fastq ]]; then
	error_exit "Your input fastq file does not exist!"
fi



# #######################
# #Step: 
if [[ ! -f $output_sample_dir/kraken_output ]]; then
	$HOME/bin/kraken   --threads 16                                                       \
	--db  /PHShome/tw786/neurogen/Tao/kraken_standard_db                                   \
	--fastq-input $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap_unmapped.fastq  \
	--only-classified-output > $output_sample_dir/kraken_output

fi

if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in kraken"
fi

if [[ ! -f $output_sample_dir/kraken_output.translate ]]; then
	$HOME/bin/kraken-translate                                                                 \
	--db /PHShome/tw786/neurogen/Tao/kraken_standard_db                                         \
	$output_sample_dir/kraken_output > $output_sample_dir/kraken_output.translate

fi


if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in kraken-translate"
fi
# $HOME/bin/kraken-translate   --mpa-format                                                      \
# --db /PHShome/tw786/neurogen/Tao/kraken_standard_db                                             \
# $output_sample_dir/kraken_output > $output_sample_dir/kraken_output.translate.mpa

if [[ ! -f $output_sample_dir/kraken_output.report ]]; then
	$HOME/bin/kraken-report --db /PHShome/tw786/neurogen/Tao/kraken_standard_db \
	$output_sample_dir/kraken_output > $output_sample_dir/kraken_output.report

fi


if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in kraken-report"
fi






