#!/bin/bash



error_exit()
{
	echo $1
	exit 1
}


##############################################
# Folder path defination
##############################################
sra_prefix=$1

input_sample_dir=/home/tw83/twang/GTEx/dbgap/Amygdala/dbGaP-6101/bam/$sra_prefix
output_dir=/home/tw83/twang/GTEx_output

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
# Filter duplicates(1,exact copy, threshold = 10)
# Filter Low complexity reads (dust, threshold=7)
# Trim 3'/5' PolyA/T
# Trim 5' low quality(threshold = 10)
# Filter sequence with length less than 20
##############################################
perl $HOME/bin/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $output_sample_dir/unmapped.fastq \
     -out_format 3 -out_good $output_sample_dir/unmapped_qc -out_bad null -log               \
     -lc_method dust -lc_threshold 7 -trim_qual_right 10 -trim_tail_left 5 -trim_tail_right 5 \
     -derep 1 -derep_min 11 -min_len 20
##############################################
# Kraken
##############################################
# #######################
# #Step: 
kraken_standard_db=/home/tw83/twang/kraken_standard_db
if [[ ! -f $output_sample_dir/kraken_output ]]; then
	$HOME/bin/kraken/kraken   --threads 8                                     \
	--db  $kraken_standard_db                                                  \
	--fastq-input $output_sample_dir/unmapped_qc.fastq                             \
	> $output_sample_dir/kraken_output

fi

if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in kraken"
fi

# if [[ ! -f $output_sample_dir/kraken_output.translate ]]; then
# 	$HOME/bin/kraken-translate                                                                 \
# 	--db $kraken_standard_db                                         \
# 	$output_sample_dir/kraken_output > $output_sample_dir/kraken_output.translate

# fi


# if [[ $? -ne 0 ]]; then
# 	error_exit "ERROR exists in kraken-translate"
# fi

if [[ ! -f $output_sample_dir/kraken_output.report ]]; then
	$HOME/bin/kraken/kraken-report --db $kraken_standard_db   \
	$output_sample_dir/kraken_output > $output_sample_dir/kraken_output.report

fi


if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in kraken-report"
fi