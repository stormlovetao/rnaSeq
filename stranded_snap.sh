#!/bin/bash

hg19_index=/PHShome/tw786/neurogen/Tao/hg19/snap_index
sample_name=$1
work_dir=/PHShome/tw786/neurogen/Tao/stranded_SNAP/$sample_name
sample_dir=/PHShome/tw786/neurogen/rnaseq_PD/filtered
R1=$sample_dir/${sample_name}.R1.fastq.gz
R2=$sample_dir/${sample_name}.R2.fastq.gz
if [[ ! -f $R1 ]]; then
	echo "$R1 doesn't exist!"
	exit
fi
if [[ ! -f $R2 ]]; then
	echo "$R2 doesn't exist!"
	exit
fi
[[ -d $work_dir ]] || mkdir $work_dir
bsub -q big-multi -M 64000 -R 'rusage[mem=64000]' -oo $work_dir/job_out \
snap-aligner paired $hg19_index $R1 $R2 -o $work_dir/snap_output.bam  