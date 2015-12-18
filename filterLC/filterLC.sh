#!/bin/bash

#################################################
#
# After mapping unmapped.bam to the viruses genomes, we found huge mount of low-complexity reads occur in the snap_output.sam file.
# According to Prinseq: Low-complexity sequences are defined as having commonly found stretches of nucleotides with limited information content 
# E.g. the dinucleotide repeat CACACACACA). Such sequences can produce a large number of high-scoring but biologically insignificant results 
# in database searches. 
#
# This script aims to filter out the low-complexity reads from snap_output.sam files, and see how many reads left there.
# The results can be non-promising, but who knows without a try.
#
# Here we will use a tool named Prinseq, which require fastq file as input, so steps will be as follows:
#
# Step1: convert sam file to fastq files.
#################################################

module load tophat/default
module load samtools/default


# cmd_fun(){
# 	samtools view -b -S $1 > $2
# 	bam2fastx -q -Q -A -N -o $3 $2
# 	perl $HOME/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $3 -out_format 3 -out_good $4 -out_bad $5 -log -lc_method dust -lc_threshold 7 
# }

samples_root_dir=/data/neurogen/Tao/run_output

viruses_list_file=filterLC_viruses_list.txt



for sample_dir in $samples_root_dir/*
do
	if [[ -d $sample_dir ]]; then
		sample_name=${sample_dir##*/}
		
		if [[ $sample_name == "PD"* ]]; then
			
			while read viruse_name; do
				if [[ -d $sample_dir/$viruse_name ]]; then
					if [[ -f $sample_dir/$viruse_name/snap_output.sam ]]; then
						#samtools view -b -S $sample_dir/$viruse_name/snap_output.sam > $sample_dir/$viruse_name/snap_output.bam
						#bsub -q "short" -o $sample_dir/$viruse_name/sam2bam.jobout
						#bam2fastx -q -Q -A -N -o  $sample_dir/$viruse_name/snap_output.fastq 
						#bsub -q "medium" -oo $sample_dir/$viruse_name/lc_jobout cmd_fun $sample_dir/$viruse_name/snap_output.sam $sample_dir/$viruse_name/snap_output.bam $sample_dir/$viruse_name/snap_output.fastq \
						#$sample_dir/$viruse_name/snap_output.good $sample_dir/$viruse_name/snap_output.bad > $sample_dir/$viruse_name/lc_stdout
						bsub -q "medium" -oo $sample_dir/$viruse_name/lc_jobout ./filterLC_core.sh $sample_dir/$viruse_name/snap_output.sam $sample_dir/$viruse_name/snap_output.bam $sample_dir/$viruse_name/snap_output.fastq \
						$sample_dir/$viruse_name/snap_output.good $sample_dir/$viruse_name/snap_output.bad > $sample_dir/$viruse_name/lc_stdout
						
					fi
				fi
			done < "$viruses_list_file"
			
			
		fi
		
	fi
	
done

