#!/bin/bash

# This scrip should be run on a server with more than 64 GB mem
# Preprocess a sample(unmapped.bam)
# Step1: Filtering Low-Complexity reads using prinseq-lite-0.20.4
# Step2: Filtering reads with length < 20bp, using prinseq-lite-0.20.4
# Step3: Filtering reads from PhiX174, GI number: gi|9626372|ref|NC_001422.1, using bowtie2 --very-sensitive-local
# Step4: Filtering reads from hg38 + all contigs, using SNAP -mrl 20
# Step5: Filtering reads from ribosome RNA(SSU + LSU, downloaded from SILVA DB), using SNAP -mrl 20



#######################
#module load
module load snap/1.0.20
module load tophat/2.1.0
module load bowtie2
usage()
{
	echo "usage: `basename $0` [sample_name]"
}
error_exit()
{
	echo $1
	exit 1
}
#######################
#Input
if [[ $# -ne 1 ]]; then
	usage
	exit 1
fi
input_sample_name=$1
input_sample_dir=/data/neurogen/rnaseq_PD/run_output/$input_sample_name
if [[ ! -d $input_sample_dir ]]; then
	echo $input_sample_dir "doesn't exist, exit!"
	exit 1
fi

#Output
kraken_output_dir=/data/neurogen/Tao/kraken_output
[ -d $kraken_output_dir ] || mkdir $kraken_output_dir

output_sample_dir=$kraken_output_dir/$input_sample_name
[ -d $output_sample_dir ] || mkdir $output_sample_dir

#######################



#######################
#Step1: transfer bam to fastq using bam2fastx

if [[ ! -f $output_sample_dir/unmapped.filterLC.fastq ]] && [[ ! -f $output_sample_dir/unmapped.fastq ]]; then
	bam2fastx -q -Q -A -N -o $output_sample_dir/unmapped.fastq $input_sample_dir/unmapped.bam || error_exit "ERROR exists in bam2fastx commond!"
fi


#######################

#######################
#Step2: trim low-complexity and short(<20bp) reads using Prinseq

if [[ ! -f $output_sample_dir/unmapped.filterLC.fastq ]] && [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.fastq ]] ; then
	
	perl $HOME/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $output_sample_dir/unmapped.fastq \
     -out_format 3 -out_good $output_sample_dir/unmapped.filterLC -out_bad null -log -lc_method dust -lc_threshold 7 -min_len 20

    if [[ $? -ne 0 ]]; then
    	error_exit "ERROR exists in perl prinseq-lite.pl commond"
    else
		[ -f $output_sample_dir/unmapped.fastq ] && rm $output_sample_dir/unmapped.fastq
	fi

fi

#######################
#Step3: Trim spiked-in phi-X reads. viral folder name: Enterobacteria_phage_phiX174_sensu_lato_uid14015. GI number: gi|9626372|ref|NC_001422.1
#Use SNAP aligner
#SNAP index path: /PHShome/tw786/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_phiX174_sensu_lato_uid14015/snap_index

if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.fastq ]]; then
	# firstly, filtering reads from PhiX rapidly using SNAP -mrl 20 -F u
	phiX_snap_index_path=/PHShome/tw786/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_phiX174_sensu_lato_uid14015/snap_index
	snap-aligner single $phiX_snap_index_path -fastq $output_sample_dir/unmapped.filterLC.fastq -mrl 20 -F u \
						-o $output_sample_dir/unmapped.filterLC.filterPhiX.bam
	if [[ $? -ne 0 ]]; then
		error_exit "ERROR exists in snap-aligner single ${phiX_index_path} ..."
	fi
	# Transform bam file to fastq as input of bowtie2
	# bam2fastx -q -Q -A -N -o $output_sample_dir/unmapped.filterLC.filterPhiX_snap.fastq  $output_sample_dir/unmapped.filterLC.filterPhiX_snap.bam

	# secondly, sensitivly filtering reads again from Phix using bowtie2 --very-sensitive-local
	# remove phi-X reads again using bowtie2. That's because snap cannot well handle short reads such like reads length <= 20 or a few longer than 20
	# phiX_bowtie2_index_path=/PHShome/tw786/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_phiX174_sensu_lato_uid14015/Enterobacteria_phage_phiX174_sensu_lato_uid14015
	# bowtie2 -x $phiX_bowtie2_index_path -U $output_sample_dir/unmapped.filterLC.filterPhiX_snap.fastq \
	# 		 --no-unal -S $output_sample_dir/unmapped.filterLC.bowtie2PhiX.sam --un $output_sample_dir/unmapped.filterLC.filterPhiX.fastq --very-sensitive-local
	if [[ $? -eq 0 ]]; then
		echo "Filter PhiX contamination, finished!"
		#remove some internal files
		#rm $output_sample_dir/unmapped.filterLC.filterPhiX_snap.bam
		# [ -f $output_sample_dir/unmapped.filterLC.filterPhiX_snap.fastq ] && rm $output_sample_dir/unmapped.filterLC.filterPhiX_snap.fastq
	fi


fi

# remove reads with length < 20bp
# if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.lg20.fastq ]]; then
# 	perl $HOME/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $output_sample_dir/unmapped.filterLC.filterPhiX.fastq -min_len 20 \
# 	-out_format 3 -out_good $output_sample_dir/unmapped.filterLC.filterPhiX.lg20  \
# 	-out_bad $output_sample_dir/unmapped.filterLC.filterPhiX.sht20

# fi





#######################
#Step4: Trim reads from human hg38 genome and all unaligned contigs
if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap.bam ]]; then
	hg38_index_path=/PHShome/tw786/neurogen/Tao/hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna.snap_index/
	snap-aligner single $hg38_index_path -mrl 20 \
				-bam $output_sample_dir/unmapped.filterLC.filterPhiX.bam \
	            -o $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap.bam 
	if [[ $? -ne 0 ]]; then
		error_exit "ERROR exists in snap-aligner single ${hg38_index_path} ..."
	fi
fi

# if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap.MAPQeb10.bam ]]; then
# 	samtools view $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap.bam -q 10 -b \
# 		> $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap.MAPQeb10.bam
# fi

if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.bam ]]; then
	samtools view -b -f 4 $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap.bam \
		> $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.bam
	if [[ $? -ne 0 ]]; then
		error_exit "ERROR exists in samtools view ..."
	fi
fi



#######################
#Step5: Trim reads from rRna. rRna database is downloaded from SILVA
# First, trim LSU rRna
if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_snap.bam ]]; then
	LSUrRna_index_path=/PHShome/tw786/neurogen/Tao/silva/SILVA_123_LSURef_tax_silva_trunc.2dna.fasta.snap_index
	snap-aligner single $LSUrRna_index_path  -mrl 20 \
				-bam $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.bam \
				-o $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_snap.bam

	samtools view -b -f 4 $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_snap.bam \
		> $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_snap_unmapped.bam
fi
if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in trimming LSU rRna "
fi

# Second, trim SSU rRna

if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap.bam ]]; then
	SSUrRna_index_path=/PHShome/tw786/neurogen/Tao/silva/SILVA_123_SSURef_Nr99_tax_silva_trunc.2dna.fasta.snap_index
	snap-aligner single $SSUrRna_index_path -mrl 20 \
				-bam $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_snap_unmapped.bam \
				-o $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap.bam

	samtools view -b -f 4 $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap.bam \
		> $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap_unmapped.bam
fi
if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in trimming SSU rRna "
fi

# Note that I saparate the SSU rRna into 16 parts because the SSU database is too large for SNAP to handle
# input_bam_file=unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_snap_unmapped.bam

# if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_ssu_snap_unmapped_16.bam ]]; then
# 	for (( i = 1; i <= 16; i++ )); do
# 		SSUrRna_index_path=/PHShome/tw786/neurogen/Tao/SILVA/SILVA_123_SSUParc_tax_silva_trunc.2dna.p$i.fasta.snap_index
# 		snap-aligner single $SSUrRna_index_path -F u \
# 			-bam $output_sample_dir/$input_bam_file   \
# 			-o $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_ssu_snap_unmapped_$i.bam
# 		input_bam_file=unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_ssu_snap_unmapped_$i.bam
# 	done
# fi


if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap_unmapped.fastq ]]; then
	bam2fastx -q -Q -A -N -o $output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap_unmapped.fastq  \
			$output_sample_dir/unmapped.filterLC.filterPhiX.hg38_snap_unmapped.rrna_lsu_ssu_snap_unmapped.bam
fi

if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in bam2fastx ..."
fi








