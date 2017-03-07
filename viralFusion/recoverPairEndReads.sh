#!/bin/bash

# """
# Well, let's start working on the viral intergration detection
# First, I will try to use ViralFusionSeq
# The input of VFS can be single end reads or pair end reads, and the latter is better.
# So our purpose here is to find the missing mate pair of reads Kraken assigned to some virus.
# """

# step1 unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam --> *.fastq
# step2 find the missing pair
usage()
{
	echo "usage: `basename $0` [sample_name] "
}
#######################
#Input
if [[ $# -ne 1 ]]; then
	usage
	exit 1
fi
input_sample_name=$1

input_sample_dir=/PHShome/tw786/neurogen/Tao/kraken_output/$input_sample_name
if [[ ! -f $input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam ]]; then
	echo $input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam "don't exist, exit!"
	exit 1
fi
samtools view $input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam | cut -f 1 | sed 's/\/[1-2]*//g' | sort | uniq \
			> $input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam.uniqReadsID

if [[ ! -f /data/neurogen/rnaseq_PD/rawfiles/$input_sample_name.R1.fastq.gz ]]; then
	echo /data/neurogen/rnaseq_PD/rawfiles/$input_sample_name.R1.fastq.gz "don't exist, exit!"
	exit 1
fi

seqtk subseq /data/neurogen/rnaseq_PD/rawfiles/$input_sample_name.R1.fastq.gz $input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam.uniqReadsID \
			> $input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.R1.fq

seqtk subseq /data/neurogen/rnaseq_PD/rawfiles/$input_sample_name.R2.fastq.gz $input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam.uniqReadsID \
			> $input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.R2.fq