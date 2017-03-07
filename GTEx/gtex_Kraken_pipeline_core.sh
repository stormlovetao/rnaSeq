#!/bin/bash



error_exit()
{
	echo $1
	exit 1
}


##############################################
# Folder path defination
##############################################
# input_sample_dir=/PHShome/tw786/neurogen/external_download/gtex/dbGap/2016march/Amygdala_unmapped
# output_dir=/PHShome/tw786/neurogen/Tao/GTEx_output/Amygdala

sra_prefix=$1
input_sample_dir=$2
output_dir=$3

sra_input_bam=${sra_prefix}.unmapped.bam
if [[ ! -f $input_sample_dir/$sra_input_bam ]]; then
	echo "$input_sample_dir/$sra_input_bam not exists! Exit"
	exit 1
fi

output_sample_dir=$output_dir/$sra_prefix 
[[ -d $output_sample_dir ]] || mkdir $output_sample_dir
##############################################
# Bam to Fastq
##############################################

bam2fastx -q -Q -A -N -o $output_sample_dir/unmapped.fastq $input_sample_dir/$sra_input_bam || error_exit "ERROR exists in bam2fastx commond!" 
##############################################
# Quality control
# #Remove this control criteria Apr 14,2016-Filter duplicates(1,exact copy, threshold = 10)
# Filter Low complexity reads (dust, threshold=7)
# Trim 3'/5' PolyA/T
# Trim 5' low quality(threshold = 10)
# Filter sequence with length less than 20
##############################################
perl $HOME/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $output_sample_dir/unmapped.fastq     \
     -out_format 3 -out_good $output_sample_dir/unmapped_qc -out_bad null -log               \
     -lc_method dust -lc_threshold 7 -trim_qual_right 10 -trim_tail_left 5 -trim_tail_right 5 \
     -min_len 20


# Kraken
##############################################
# #######################
# #Step: 
if [[ ! -f $output_sample_dir/kraken_output ]]; then
	$HOME/bin/kraken   --threads 8                                                            \
	--db  /PHShome/tw786/neurogen/Tao/kraken_standard_db                                       \
	--fastq-input $output_sample_dir/unmapped_qc.fastq                                          \
	> $output_sample_dir/kraken_output
fi

if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in kraken"
else
	echo "kraken finished!"
fi


if [[ ! -f $output_sample_dir/kraken_output.report ]]; then
	$HOME/bin/kraken-report --db /PHShome/tw786/neurogen/Tao/kraken_standard_db \
	$output_sample_dir/kraken_output > $output_sample_dir/kraken_output.report
fi


if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in kraken-report"
else
	echo "kraken-report finished!"
fi

##############################################
# gunzip
gzip $output_sample_dir/unmapped.fastq
gzip $output_sample_dir/unmapped_qc.fastq
##############################################


