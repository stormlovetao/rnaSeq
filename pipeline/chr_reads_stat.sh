#!/bin/bash
######################
# After mapping reads to hg38, we
# Count how many reads are mapped to every chromosome
#
######################
usage()
{
	echo "usage: sh  `basename $0` [bam_file_path]"
}


if [[ $# -ne 1 ]]; then
	usage 
	exit 1
fi
if [[ ! -f $1 ]]; then
	echo "$1 don't exist!"
	exit 1
fi

chr_size_file=`dirname $0`/hg38.chrom.sizes

if [[ ! -f $chr_size_file ]]; then
	echo "$chr_size_file doesn't exist!"
	exit 1
fi


output_file_path=$(dirname $1)/chr_reads_stat.txt

samtools view $1 | cut -f 3 > $(dirname $1)/chr_reads.tmp

while  read line ; do

	chr_name=$(echo $line | cut -d ' ' -f1)
	reads_count=fgrep $chr_name $(dirname $1)/chr_reads.tmp | wc -l
	echo $chr_name $reads_count >> $output_file_path

done < "$chr_size_file"

rm $(dirname $1)/chr_reads.tmp