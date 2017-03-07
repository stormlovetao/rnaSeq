#!/bin/bash
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


if [[ $# -ne 1 ]]; then
	usage
	exit 1
fi
input_sample_name=$1

input_sample_dir=/data/neurogen/rnaseq_PD/run_output/$input_sample_name
if [[ ! -d $input_sample_dir ]]; then
	echo $input_sample_dir "doesn't exist, exit!"
	exit 1
fi
output_sample_dir=/data/neurogen/Tao/mappedKraken/
[[ -f $output_sample_dir ]] || mkdir $output_sample_dir

bam2fastx -q -Q -M -N -o $output_sample_dir/accepted_hits.fastq $input_sample_dir/accepted_hits.bam || error_exit "ERROR exists in bam2fastx commond!"


# #######################
# #Step: 
if [[ ! -f $output_sample_dir/kraken_output ]]; then
	$HOME/bin/kraken   --threads 8                                                                       \
	--db  /PHShome/tw786/neurogen/Tao/kraken_standard_db                                                  \
	--fastq-input $output_sample_dir/accepted_hits.fastq          \
	--unclassified-out $output_sample_dir/kraken_output_unclassified                                        \
	> $output_sample_dir/kraken_output

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