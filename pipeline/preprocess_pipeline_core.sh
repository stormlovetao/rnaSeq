#!/bin/bash

# This script should be run on a server with more than 64 GB mem
# Preprocess a sample(unmapped.bam)
# Step1: Filtering Low-Complexity reads and reads with length < 20bp using prinseq-lite-0.20.4
# Step2: Filtering reads from PhiX174, GI number: gi|9626372|ref|NC_001422.1, using SNAP (bowtie2 --very-sensitive-local)
# Step3: Filtering reads from hg38 + all contigs, using SNAP -mrl 20
# Step4: Filtering reads from ribosome RNA(SSU + LSU, downloaded from SILVA DB), using SNAP -mrl 20

# Step5: CirRna
# Step6: gene fusion

# bsub -q big -M 64000 -R 'rusage[mem=64000]' ./preprocess_pipeline_core.sh samplename

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
#Step1: transfer bam to fastq using bam2fastx( or samtools bam2fq )
#######################


if [[ ! -f $output_sample_dir/unmapped.filterLC.fastq ]] && [[ ! -f $output_sample_dir/unmapped.fastq ]]; then
	bam2fastx -q -Q -A -N -o $output_sample_dir/unmapped.fastq $input_sample_dir/unmapped.bam || error_exit "ERROR exists in bam2fastx commond!"
fi




#######################
#Step2: trim low-complexity and short(<20bp) reads using Prinseq
#######################


if [[ ! -f $output_sample_dir/unmapped.filterLC.fastq ]]; then
	
	perl $HOME/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $output_sample_dir/unmapped.fastq \
     -out_format 3 -out_good $output_sample_dir/unmapped.filterLC -out_bad null -log -lc_method dust -lc_threshold 7 -min_len 20

    if [[ $? -ne 0 ]]; then
    	error_exit "ERROR exists in perl prinseq-lite.pl commond"
    else
    	echo "Low-Complexity reads and short reads(<20bp) are filtered!"
	fi

fi


#######################
#Step3: Filter spiked-in phi-X reads. viral folder name: Enterobacteria_phage_phiX174_sensu_lato_uid14015. GI number: gi|9626372|ref|NC_001422.1
#######################
#Use SNAP aligner
#SNAP index path: /PHShome/tw786/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_phiX174_sensu_lato_uid14015/snap_index

if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.bam ]]; then
	# firstly, filtering reads from PhiX rapidly using SNAP -mrl 20 -F u
	phiX_snap_index_path=/PHShome/tw786/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_phiX174_sensu_lato_uid14015/snap_index
	snap-aligner single $phiX_snap_index_path -fastq $output_sample_dir/unmapped.filterLC.fastq -mrl 20 -F u \
						-o $output_sample_dir/unmapped.filterLC.filterPhiX.bam
	if [[ $? -ne 0 ]]; then
		error_exit "ERROR exists in snap-aligner single ${phiX_index_path} ..."
	fi
	if [[ $? -eq 0 ]]; then
		echo "Filter PhiX contamination, finished!"
	fi


fi


#######################
#Step4: Filter reads from human hg38 genome and all unaligned contigs
#######################

if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.bam ]]; then
	hg38_index_path=/PHShome/tw786/neurogen/Tao/hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna.snap_index/
	snap-aligner single $hg38_index_path -mrl 20 -F u \
				-bam $output_sample_dir/unmapped.filterLC.filterPhiX.bam \
	            -o $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.bam 
	if [[ $? -ne 0 ]]; then
		error_exit "ERROR exists in snap-aligner single ${hg38_index_path} ..."
	else
		echo "Reads from hg38(reference + unaligned contigs) are filtered!"
	fi
fi






#######################
#Step5: Filter reads from rRna. rRna database is downloaded from SILVA
#######################

# First, filter LSU rRna
if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU.bam ]]; then
	LSUrRna_index_path=/PHShome/tw786/neurogen/Tao/silva/SILVA_123_LSURef_tax_silva_trunc.2dna.fasta.snap_index
	snap-aligner single $LSUrRna_index_path  -mrl 20 -F u \
				-bam $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.bam \
				-o $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU.bam
	if [[ $? -ne 0 ]]; then
		error_exit "ERROR exists in trimming LSU rRna "
	else
		echo "Reads from rRna(LSU) are filtered"
	fi
	
fi


# Second, filter SSU rRna

if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam ]]; then
	SSUrRna_index_path=/PHShome/tw786/neurogen/Tao/silva/SILVA_123_SSURef_Nr99_tax_silva_trunc.2dna.fasta.snap_index
	snap-aligner single $SSUrRna_index_path -mrl 20 -F u \
				-bam $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU.bam \
				-o $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam
	if [[ $? -ne 0 ]]; then
		error_exit "ERROR exists in trimming SSU rRna "
	else
		echo "Reads from rRna(SSU) are filtered"
	fi
	
fi


#######################
#Step 6 transform bam to fastq
#######################


if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.fastq ]]; then
	bam2fastx -q -Q -A -N -o $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.fastq  \
			$output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam
fi

if [[ $? -ne 0 ]]; then
	error_exit "ERROR exists in bam2fastx ..."
else
	echo "Transform Bam to fastq for Kraken process"
fi

#######################
#Step 7 prepare read lists supporting cirRna and gene fusion
# this two read lists are from xianjun's work and allessy's work
#######################


if [[ ! -f $output_sample_dir/AllCircRnaList.txt ]]; then
	python /PHShome/tw786/MyOwnScript/getCircRna.py $input_sample_name
	if [[ $? -ne 0 ]]; then
	 	echo "ERROR exists in getCircRna.py script"
	fi 
		
fi

if [[ ! -f $output_sample_dir/AllGeneFusionReads.txt ]]; then
	gene_fusion_bamfile=/PHShome/tw786/neurogen/rnaseq_PD/run_output/$input_sample_name/tophat_fusion/accepted_hits.bam
	samtools view $gene_fusion_bamfile | cut -f 1 > $output_sample_dir/AllGeneFusionReads.txt
	if [[ $? -ne 0 ]]; then
		echo "ERROR exists in get reads from tophat_fusion"
	fi
fi


##################
# Run Kraken 

bsub -q big-multi -n 8 -M 128000 -R 'rusage[mem=128000]' sh /PHShome/tw786/MyOwnScript/kraken_pipeline_core.sh $input_sample_name



































