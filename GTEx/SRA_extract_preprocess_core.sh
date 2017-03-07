#!/bin/bash


sra_file=$1
out_dir=$2
sra_prefix=$3
work_dir=$4

error_exit()
{
	echo $1
	exit 1
}

cd $work_dir
sam-dump -r -u $sra_file | samtools view -b -S - > $out_dir/${sra_prefix}.bam
if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in sam-dump commond"
fi
samtools view -f 4 $out_dir/${sra_prefix}.bam | samtools view -b -S - > $out_dir/${sra_prefix}.unmapped.bam
[[ -f $out_dir/${sra_prefix}.bam ]] && rm $out_dir/${sra_prefix}.bam
##############################################
# Bam to Fastq
##############################################

# bam2fastx -q -Q -A -N -o $out_dir/${sra_prefix}.unmapped.fastq $out_dir/${sra_prefix}.unmapped.bam || error_exit "ERROR exists in bam2fastx commond!" 

##############################################
# Quality control
# Filter duplicates(1,exact copy, threshold = 10)
# Filter Low complexity reads (dust, threshold=7)
# Trim 3'/5' PolyA/T
# Trim 5' low quality(threshold = 10)
# Filter sequence with length less than 20
##############################################
# perl $HOME/bin/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $out_dir/${sra_prefix}.unmapped.fastq \
#      -out_format 3 -out_good $out_dir/${sra_prefix}.unmapped_qc -out_bad null -log               \
#      -lc_method dust -lc_threshold 7 -trim_qual_right 10 -trim_tail_left 5 -trim_tail_right 5     \
#      -derep 1 -derep_min 11 -min_len 20	