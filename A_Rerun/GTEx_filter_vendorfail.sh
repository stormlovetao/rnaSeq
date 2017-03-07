#!/bin/bash
# for each sample, filter vendor failed reads from unmapped_qc.fastq.gz
prog_dir=/home/tw83/Myscript
bam_dir=/home/tw83/twang/GTEx/dbgap/unmapped
fastq_dir=/home/tw83/twang/GTEx_output
tissue=$1
SRR_name=$2

function error_exit
{
	echo "$1" 1>&2
	exit 1
}

bam_file=$bam_dir/$tissue/${SRR_name}.unmapped.bam
fastq_gz=$fastq_dir/$tissue/$SRR_name/unmapped_qc.fastq.gz

if [[ ! -f $bam_file ]]; then
	error_exit "Not found:" $bam_file
fi
if [[ ! -f $fastq_gz ]]; then
	error_exit "Not found:" $fastq_gz
fi

cd $fastq_dir/$tissue/$SRR_name

# #extract vendor-failed reads
# samtools view -f 516 $bam_file | samtools view -b -S - > vendor_failed.bam
# if [[ $? -ne 0 ]]; then
# 	error_exit "samtools view failed: cannot output vendor_failed.bam"
# fi
# bam2fastx -q -Q -A -N -o vendor_failed.fastq vendor_failed.bam
# if [[ $? -ne 0 ]]; then
# 	error_exit "bam2fastx failed: cannot convert vendor_failed.bam to vendor_failed.fastq"
# fi

# awk 'NR%4==1' vendor_failed.fastq | sed 's/@//' > vendor_failed_list

# [[ -f vendor_failed.bam ]] && rm vendor_failed.bam
# [[ -f vendor_failed.fastq ]] && rm vendor_failed.fastq

#Transfer unmapped_qc.fastq.gz to fa, during which remove vendor failed reads
zcat $fastq_gz | python $prog_dir/GTEx_filter_vendorfail_fq2fa.py vendor_failed_list > unmapped_qc.fa
if [[ $? -ne 0 ]]; then
	error_exit "GTEx_filter_vendorfail_fq2fa.py failed"
fi

[[ -f unmapped_qc.fa.gz ]] && rm unmapped_qc.fa.gz
gzip unmapped_qc.fa
if [[ $? -ne 0 ]]; then
	error_exit "gzip unmapped_qc.fa failed"
else
	touch _GTEx_filter_vendorfail_finished
fi
# zcat unmapped_qc.fastq.gz | \
# awk '
# NR%4==1 { readID = substr($0,2);
# 		  if (grep readID ) {printf(">%s_%s\n", srr, substr($0,2))} } 
# NR%4==2 {print $0}
# ' > unmapped_qc.fa

