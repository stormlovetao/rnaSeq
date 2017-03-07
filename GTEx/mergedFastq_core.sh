#!/bin/bash
"""

"""
output_sample_dir=$1
unmapped_bam_file_dir=$2
srr_name=$3
#bam2fastq
bam2fastx -q -Q -A -N -o $output_sample_dir/unmapped.fastq $unmapped_bam_file_dir 
if [[ $? -ne 0 ]]; then
	echo "ERROR exists in bam2fastx commond!" 
	exit 1	
fi
##############################################
# Quality control
# #Remove this control criteria Apr 14,2016-Filter duplicates(1,exact copy, threshold = 10)
# Filter Low complexity reads (dust, threshold=7)
# Trim 3'/5' PolyA/T
# Trim 5' low quality(threshold = 10)
# Filter sequence with length less than 20
##############################################
prinseq_pl=$HOME/prinseq-lite-0.20.4/prinseq-lite.pl
if [[ ! -f $prinseq_pl ]]; then
	prinseq_pl=$HOME/bin/prinseq-lite-0.20.4/prinseq-lite.pl
	[[ ! -f $prinseq_pl ]] && echo "cannot find " $prinseq_pl && exit 1 
fi
perl $prinseq_pl -fastq $output_sample_dir/unmapped.fastq     \
     -out_format 3 -out_good $output_sample_dir/unmapped_qc -out_bad null -log               \
     -lc_method dust -lc_threshold 7 -trim_qual_right 10 -trim_tail_left 5 -trim_tail_right 5 \
     -min_len 20
if [[ $? -ne 0 ]]; then
	echo "ERROR in $prinseq_pl exit"
	exit 1
fi
##############################################
# fastq 2 fasta
# add tag(srr_name) to seq identity
cd $output_sample_dir
awk -v srr="$srr_name" 'NR%4==1 {printf(">%s_%s\n", srr, substr($0,2))} NR%4==2 {print $0}' unmapped_qc.fastq > unmapped_qc.fa
if [[ $? -ne 0 ]]; then
	echo "ERROR in fastq2fasta(awk) exit"
	exit 1
fi
##############################################
#zip fastq, fasta
gzip unmapped.fastq
gzip unmapped_qc.fastq
gzip unmapped_qc.fa