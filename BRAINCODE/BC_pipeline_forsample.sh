#!/bin/bash
# adopt gtex pipeline for braincode data
# add one step, filter PhiX174, using SNAP(mismatch= 4) 
#Reference: /pipeline/preprocess_pipeline_core.sh
#######################


input_sample_name=$1
input_sample_dir=/data/neurogen/rnaseq_PD/run_output/$input_sample_name
if [[ ! -d $input_sample_dir ]]; then
	echo $input_sample_dir "doesn't exist, exit!"
	exit 1
fi

#Output
kraken_output_dir=/data/neurogen/Tao/BRAINCODE_output
[ -d $kraken_output_dir ] || mkdir $kraken_output_dir

output_sample_dir=$kraken_output_dir/$input_sample_name
[ -d $output_sample_dir ] || mkdir $output_sample_dir


#######################
#Step1: transfer bam to fastq using bam2fastx( or samtools bam2fq )
#######################


if  [[ ! -f $output_sample_dir/unmapped.fastq ]]; then
	bam2fastx -q -Q -A -N -o $output_sample_dir/unmapped.fastq $input_sample_dir/unmapped.bam
	if [[ $? -ne 0 ]]; then
		echo "Error: bam2fastx"
		exit 1
	else
		echo "bam2fastx Success!"
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
cd $output_sample_dir
if [[ ! -f unmapped_qc.fastq ]] && [[ ! -f unmapped_qc.fastq.gz ]]; then
	prinseq_pl=$HOME/prinseq-lite-0.20.4/prinseq-lite.pl
	if [[ ! -f $prinseq_pl ]]; then
		prinseq_pl=$HOME/bin/prinseq-lite-0.20.4/prinseq-lite.pl
		[[ ! -f $prinseq_pl ]] && echo "cannot find " $prinseq_pl && exit 1 
	fi
	
	perl $prinseq_pl -fastq unmapped.fastq \
	    -out_format 3 -out_good unmapped_qc -out_bad null -log \
	    -lc_method dust -lc_threshold 7 -trim_qual_right 10 -trim_tail_left 5 -trim_tail_right 5 \
	    -min_len 20
	if [[ $? -ne 0 ]]; then
		echo "Error: prinseq_pl"
		exit 1
	else
		echo "QC Success!"
	fi
fi

###gzip###
[[ -f unmapped.fastq ]] && rm unmapped.fastq
[[ -f unmapped_qc.fastq ]] && gzip unmapped_qc.fastq


#######################
#Step3: Filter spiked-in phi-X reads. viral folder name: Enterobacteria_phage_phiX174_sensu_lato_uid14015. GI number: gi|9626372|ref|NC_001422.1
#######################
#Use SNAP aligner
#SNAP index path: /PHShome/tw786/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_phiX174_sensu_lato_uid14015/snap_index
if [[ ! -f unmapped_qc.filterPhiX.bam ]]; then
	# filter reads from PhiX using SNAP -mrl 20 -F u -d 4(maximum edit distance = 4)
	phiX_snap_index_path=/PHShome/tw786/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_phiX174_sensu_lato_uid14015/snap_index
	snap-aligner single $phiX_snap_index_path \
					-compressedFastq unmapped_qc.fastq.gz \
					-mrl 20 -d 4 -F u \
					-o unmapped_qc.filterPhiX.bam \
					> unmapped_qc.filterPhiX.log
	if [[ $? -ne 0 ]]; then
		echo "Error: snap-aligner failed to filter phi-X reads"
		exit 1
	else
		echo "snap-aligner filter phi-X Successfully!"
	fi
fi


#######################
#transform bam to fastq
#######################


if [[ ! -f unmapped_qc.filterPhiX.fastq ]]; then
	bam2fastx -q -Q -A -N -o  unmapped_qc.filterPhiX.fastq unmapped_qc.filterPhiX.bam
	if [[ $? -ne 0 ]]; then
		echo "Error: bam2fastx failed to transform bam to fastq"
		exit 1
	else
		echo "bam2fastx transform bam to fastq Successfully!" 
	fi
fi

########gzip#########
[[ -f unmapped_qc.filterPhiX.fastq ]] && gzip unmapped_qc.filterPhiX.fastq

######################################################
##############################################
# Kraken
##############################################
# #######################

kraken_standard_db=/PHShome/tw786/neurogen/Tao/kraken_standard_db
kraken_bin=kraken
kraken_report_bin=kraken-report
if [[ ! -f kraken_output ]]; then
	$kraken_bin  --threads 4                              \
	--db  $kraken_standard_db                              \
	--fastq-input unmapped_qc.filterPhiX.fastq.gz           \
	--gzip-compressed							             \
	> kraken_output
fi

if [[ $? -ne 0 ]]; then
	echo "ERROR exists in kraken"
fi
if [[ ! -f kraken_output.report ]]; then
	$kraken_report_bin --db $kraken_standard_db kraken_output > kraken_output.report
fi
