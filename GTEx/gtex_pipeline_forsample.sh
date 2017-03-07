#!/bin/bash
#bsub command                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
# while read -r line;do if [[ ! -d /home/tw83/twang/GTEx_output/Liver/$line ]];then bsub -q short -W 12:00 bash /home/tw83/Myscript/gtex_pipeline_forsample.sh $line ; fi;break; done < SRR_list.txt
SRR_name=$1

work_dir=/home/tw83/twang/GTEx/dbgap/Liver/dbGaP-6101
log_dir=/home/tw83/twang/GTEx_output/Liver/Log
out_dir=/home/tw83/twang/GTEx_output/Liver/$SRR_name
[[ -d $out_dir ]] || mkdir $out_dir
[[ -d $log_dir ]] || mkdir $log_dir
cd $work_dir 
####DOWNLOAD####
if [[ ! -f $work_dir/sra/${SRR_name}.sra ]]; then
	prefetch -X 30000000 $SRR_name
	if [[ $? -ne 0 ]]; then
		echo "Error: prefetch" > $log_dir/${SRR_name}.log
		exit 1
	else
		echo "DOWNLOAD Success!"
	fi
fi


####DECOMPRESS####
if [[ ! -f $out_dir/${SRR_name}.unmapped.bam ]]; then
	sam-dump -r -u ./sra/${SRR_name}.sra | samtools view -b -S - > $out_dir/${SRR_name}.bam
	if [[ $? -ne 0 ]]; then
		echo "Error: sam-dump" > $log_dir/${SRR_name}.log
		exit 1
	fi
	samtools view -f 4 $out_dir/${SRR_name}.bam | samtools view -b -S - > $out_dir/${SRR_name}.unmapped.bam
	if [[ $? -ne 0 ]]; then
		echo "Error: samtools view" > $log_dir/${SRR_name}.log
		exit 1
	else
		echo "DECOMPRESS Success!"
	fi
	[[ -f $out_dir/${SRR_name}.bam ]] && rm $out_dir/${SRR_name}.bam
fi


##############################################
# Bam to Fastq
##############################################
cd $out_dir
if [[ ! -f unmapped.fastq ]]; then
	bam2fastx -q -Q -A -N -o unmapped.fastq ${SRR_name}.unmapped.bam 
	if [[ $? -ne 0 ]]; then
		echo "Error: bam2fastx" > $log_dir/${SRR_name}.log
		exit 1
	else
		echo "Bam2Fastq Success!"
	fi
fi


##############################################
# Quality control
# #Remove this control criteria Apr 14,2016-Filter duplicates(1,exact copy, threshold = 10)
# Filter Low complexity reads (dust, threshold=7)
# Trim 3'/5' PolyA/T
# Trim 5' low quality(threshold = 10)
# Filter sequence with length less than 20
##############################################
if [[ ! -f unmapped_qc.fastq ]] && [[ ! -f unmapped_qc.fastq.gz ]]; then
	prinseq_pl=$HOME/prinseq-lite-0.20.4/prinseq-lite.pl
	if [[ ! -f $prinseq_pl ]]; then
		prinseq_pl=$HOME/bin/prinseq-lite-0.20.4/prinseq-lite.pl
		[[ ! -f $prinseq_pl ]] && echo "cannot find " $prinseq_pl && exit 1 
	fi
	if [[ ! -f unmapped.fastq ]]; then
		echo "Error: unmapped.fastq doesn't exist" > $log_dir/${SRR_name}.log
		exit 1
	fi
	perl $prinseq_pl -fastq unmapped.fastq \
	    -out_format 3 -out_good unmapped_qc -out_bad null -log \
	    -lc_method dust -lc_threshold 7 -trim_qual_right 10 -trim_tail_left 5 -trim_tail_right 5 \
	    -min_len 20
	if [[ $? -ne 0 ]]; then
		echo "Error: prinseq_pl" > $log_dir/${SRR_name}.log
		exit 1
	else
		echo "QC Success!"
	fi
fi


##############################################
# fastq 2 fasta
# add tag(srr_name) to seq identity
if [[ ! -f unmapped_qc.fa ]] && [[ ! -f unmapped_qc.fa.gz ]]; then
	awk -v srr="$SRR_name" 'NR%4==1 {printf(">%s_%s\n", srr, substr($0,2))} NR%4==2 {print $0}' unmapped_qc.fastq > unmapped_qc.fa
	if [[ $? -ne 0 ]]; then
		echo "Error: awk" > $log_dir/${SRR_name}.log
		exit 1
	else
		echo "fastq2fasta Success!"
	fi
fi

##############################################
#zip fastq, fasta
gzip unmapped.fastq
gzip unmapped_qc.fastq
gzip unmapped_qc.fa
echo "Success" > $log_dir/${SRR_name}.log
echo "Finished!"

