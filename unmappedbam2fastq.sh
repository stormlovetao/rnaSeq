#!/bin/bash

#########################
#
#Transform unmapped Bam file to fastq file.
#
#########################
module load tophat/2.1.0 #newest on server

#bsub bam2fastx -q -Q -A -N -o unmapped.fastq unmapped.bam

input_dir=/data/neurogen/rnaseq_PD/run_output

output_dir=/data/neurogen/Tao/SamplesFastq
[ -d $output_dir ] || mkdir $output_dir

cd $input_dir
for folder in $input_dir/*;
do
	if test -d $folder
		then
		folder_name=${folder##*/}
		#echo $folder
		#echo $folder_name
		new_sample_dir=$output_dir/$folder_name
		#echo $new_sample_dir
		
		#echo ${folder_name}"_unmapped.fastq"
		[ -d $new_sample_dir ] || mkdir $new_sample_dir
		if test -f $new_sample_dir/${folder_name}"_unmapped.fastq"
			then 
			continue
		else
			bsub -q "short" -o $new_sample_dir/bam2fq_jobout bam2fastx -q -Q -A -N -o $new_sample_dir/${folder_name}"_unmapped.fastq" $folder/unmapped.bam
		fi
	fi
	
	
done
exit