#!/bin/bash
# write the unmapped reads number into unmapped_reads.txt
sample_dir=/PHShome/tw786/neurogen/external_download/gtex/dbGap/2016march/unmapped
out_file=/PHShome/tw786/neurogen/external_download/gtex/dbGap/2016march/unmapped_reads.txt
cd $sample_dir
echo "Sample unmapped_reads_number" > $out_file
for sample_path in $sample_dir/*;do
	sample_name=${sample_path##*/}
	reads_num=`samtools view -c $sample_name`
	echo $sample_name $reads_num >> $out_file
done
	