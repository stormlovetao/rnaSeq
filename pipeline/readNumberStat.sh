#!/bin/bash

# outputfile: readsCake.numstat

assert()
{
	if [[ ! -f $1 ]]; then
		echo $1 "does not exist"
		exit 1
	fi
}
#Input
if [[ $# -ne 1 ]]; then
	usage
	exit 1
fi
input_sample_name=$1
input_sample_dir=/data/neurogen/Tao/kraken_output/$input_sample_name
if [[ ! -d $input_sample_dir ]]; then
	echo $input_sample_dir "doesn't exist, exit!"
	exit 1
fi

unmapped_fastq_log=$input_sample_dir/unmapped.fastq.log
assert $unmapped_fastq_log

unmapped_filterLC_filterPhiX_bam=$input_sample_dir/unmapped.filterLC.filterPhiX.bam
assert $unmapped_filterLC_filterPhiX_bam

unmapped_filterLC_filterPhiX_filterHg38_bam=$input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.bam
assert $unmapped_filterLC_filterPhiX_filterHg38_bam

unmapped_filterLC_filterPhiX_filterHg38_filterLSU_SSU_bam=$input_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam
assert $unmapped_filterLC_filterPhiX_filterHg38_filterLSU_SSU_bam

kraken_output_report=$input_sample_dir/kraken_output.report
assert $kraken_output_report

total_reads=`grep 'Input sequences:.*' $unmapped_fastq_log -o | awk '{ print $3 }'|sed 's/,//g' `

LCS=`grep 'Bad sequences:.*' $unmapped_fastq_log -o | awk '{ print $3 }'|sed 's/,//g' `

filter_phiX_left=`samtools view -c $unmapped_filterLC_filterPhiX_bam`

filter_hg38_left=`samtools view -c $unmapped_filterLC_filterPhiX_filterHg38_bam`

filter_rrna_left=`samtools view -c $unmapped_filterLC_filterPhiX_filterHg38_filterLSU_SSU_bam`

kraken_unclassified=`grep 'unclassified$' $kraken_output_report | awk '{print $2}'` 

Homo_sapiens=`grep 'Homo sapiens' $kraken_output_report | awk '{print $2}' `

Bacteria=`grep 'D.*Bacteria$' $kraken_output_report | awk '{print $2}'`

Viruses=`grep 'D.*Viruses' $kraken_output_report | awk '{print $2}'`

output_stat_file=$input_sample_dir/readsCake.numstat

echo -e 'samtools view -c' $unmapped_filterLC_filterPhiX_bam '\t' $filter_phiX_left > $output_stat_file
echo -e 'samtools view -c' $unmapped_filterLC_filterPhiX_filterHg38_bam '\t' $filter_hg38_left >> $output_stat_file
echo -e 'samtools view -c' $unmapped_filterLC_filterPhiX_filterHg38_filterLSU_SSU_bam '\t' $filter_rrna_left >> $output_stat_file

echo -e 'total_reads\t'$total_reads >> $output_stat_file

echo -e 'LCS_reads\t'$LCS >> $output_stat_file  

((phiX = $total_reads - $LCS - $filter_phiX_left))

echo -e 'PhiX_reads\t'$phiX >> $output_stat_file

(( hg38 = $filter_phiX_left - $filter_hg38_left))

echo -e 'hg38_reads\t'$hg38 >> $output_stat_file

((rrna = $filter_hg38_left - filter_rrna_left))

echo -e 'rrna_reads\t'$rrna >> $output_stat_file

echo -e 'kraken_unclassified_reads\t'$kraken_unclassified >> $output_stat_file
echo -e 'kraken_Homo_sapiens_reads\t'$Homo_sapiens >> $output_stat_file
echo -e 'kraken_Bacteria_reads\t'$Bacteria >> $output_stat_file
echo -e 'kraken_Viruses_reads\t'$Viruses >> $output_stat_file

kraken_unclassified_cirRna=$input_sample_dir/kraken.unclassified.cirRna.txt
if [[ -f $kraken_unclassified_cirRna ]]; then
	kraken_unclassified_cirRna_reads=`wc -l $kraken_unclassified_cirRna | awk '{print $1}' `
	echo -e 'kraken_unclassified_cirRna_reads\t'$kraken_unclassified_cirRna_reads >> $output_stat_file
else
	echo -e 'kraken_unclassified_cirRna_reads\t-1' >> $output_stat_file
fi



kraken_unclassified_fusion=$input_sample_dir/kraken.unclassified.fusion.txt

if [[ -f $kraken_unclassified_fusion ]]; then
	kraken_unclassified_fusion_reads=`wc -l $kraken_unclassified_fusion | awk '{print $1}'`
	echo -e 'kraken_unclassified_fusion_reads\t'$kraken_unclassified_fusion_reads >> $output_stat_file
else
	echo -e 'kraken_unclassified_fusion_reads\t-1' >> $output_stat_file
fi


kraken_human9606_cirRna=$input_sample_dir/kraken.human9606.cirRna.txt
if [[ -f $kraken_human9606_cirRna ]]; then
	kraken_human9606_cirRna_reads=`wc -l $kraken_human9606_cirRna | awk '{print $1}' `
	echo -e 'kraken_human9606_cirRna_reads\t'$kraken_human9606_cirRna_reads >> $output_stat_file
else
	echo -e 'kraken_human9606_cirRna_reads\t-1' >> $output_stat_file
fi


kraken_human9606_fusion=$input_sample_dir/kraken.human9606.fusion.txt
if [[ -f $kraken_human9606_fusion ]]; then
	kraken_human9606_fusion_reads=`wc -l $kraken_human9606_fusion | awk '{print $1}'`
	echo -e 'kraken_human9606_fusion_reads\t'$kraken_human9606_fusion_reads >> $output_stat_file
else
	echo -e 'kraken_human9606_fusion_reads\t-1' >> $output_stat_file
fi





















