#!/bin/bash

# - Detect mixed reads from kraken's output.
#     a. do we need to check specificity of those mixed reads?
# -  Find the virus-human pairs. Virus singleton read can be extracted from kraken output, and then check if their pairs can be mapped to human reference using tophat.
# - De novo assemble read pairs from step 1 and 2, output contigs.
# - blast those contigs from step 3 against nt database to check if there is mixed contig exist. And check its specificity.
#######################
#Input
usage()
{
	echo "usage: `basename $0` [sample_name] [virus_taxid]"
}
if [[ $# -ne 2 ]]; then
	usage
	exit 1
fi
module_name=`basename $0`
script_dir=`pwd`
sample_name=$1
homo_sapiens_taxid=9606
virus_taxid=$2
#sample_dir=/PHShome/tw786/neurogen/Tao/BRAINCODE_output/$sample_name
sample_dir=/PHShome/tw786/neurogen/Tao/dbGaP/public/sra
work_folder=virus_$virus_taxid
work_dir=${sample_dir}/$work_folder
[[ -d $work_dir ]] || mkdir $work_dir
logfile=$work_dir/log.txt
echo "#Start--->>>" > $logfile
###############################Step 1 ###############################
# Extract mixed reads from kraken's output, finished in step2
# cd $work_dir
# awk -v v1=$homo_sapiens_taxid -v v2=$virus_taxid -F '\t' '($5 ~ v1) && ($5 ~ v2)' kraken_output | \
# 		cut -f 2 > mixedReadsList.txt

###############################Step 2 ###############################
# read kraken_output, and output mixed_reads, virus_virus reads, virus_human reads, virus_unclassify reads, virus_unsure reads.

python findVirusPair.py $homo_sapiens_taxid $virus_taxid $sample_dir/kraken_output $work_dir
findVirusPair_output_dir=$work_dir/virus_pair_catogories

cd $findVirusPair_output_dir
mixed_reads_count=$(wc -l mixed_reads.txt | awk '{print $1}')
echo "mixed_reads: $mixed_reads_count" >> $logfile
virus_human_reads_count=$(wc -l virus_human_reads.txt | awk '{print $1}')
echo "virus_human_reads: $virus_human_reads_count" >> $logfile
virus_unsure_reads_left_count=$(wc -l virus_unsure_reads_left.txt | awk '{print $1}')
virus_unsure_reads_right_count=$(wc -l virus_unsure_reads_right.txt | awk '{print $1}')
virus_unsure_reads_count=$(($virus_unsure_reads_left_count + $virus_unsure_reads_right_count))
echo "virus_unsure_reads: $virus_unsure_reads_count" >> $logfile
virus_virus_reads_count=$(wc -l virus_virus_reads.txt | awk '{print $1}')
echo "virus_virus_reads: $virus_virus_reads_count" >> $logfile
cat mixed_reads.txt virus_human_reads.txt virus_unsure_reads_left.txt virus_unsure_reads_right.txt > mixed_virusHuman_virusUnsure.txt
mixed_virusHuman_virusUnsure_count=$(wc -l mixed_virusHuman_virusUnsure.txt | awk '{print $1}')
echo "mixed_virusHuman_virusUnsure_count: $mixed_virusHuman_virusUnsure_count" >> $logfile
if (($mixed_virusHuman_virusUnsure_count < 10)); then
	echo "#Exit: There are only $mixed_virusHuman_virusUnsure_count (<10) reads in mixed_virusHuman_virusUnsure.txt, Exit!" >> $logfile
	exit
fi
if [[ ! -f $sample_dir/${sample_name}_1.fastq.gz ]]; then
	echo "#Exit: there is no $sample_dir/${sample_name}_1.fastq.gz, exit" >> $logfile
	exit
fi
if [[ ! -f mixed_virusHuman_virusUnsure.txt.R1.fq ]]; then
	seqtk subseq $sample_dir/${sample_name}_1.fastq.gz mixed_virusHuman_virusUnsure.txt > mixed_virusHuman_virusUnsure.txt.R1.fq
	seqtk subseq $sample_dir/${sample_name}_2.fastq.gz mixed_virusHuman_virusUnsure.txt > mixed_virusHuman_virusUnsure.txt.R2.fq
	# gzip mixed_virusHuman_virusUnsure.txt.R1.fq
	# gzip mixed_virusHuman_virusUnsure.txt.R2.fq
fi

# get reads from source data
cd $work_dir


# Trinity
module load bowtie
Trinity --seqType fq --JM 4G --left mixed_virusHuman_virusUnsure.txt.R1.fq --right mixed_virusHuman_virusUnsure.txt.R1.fq --CPU 4
if [[ ! -f ./trinity_out_dir/Trinity.fasta ]]; then
	echo "#Exit: Trinity failed!"
	echo "#Exit: Trinity failed!" >> $logfile
	exit
fi
bowtie_PE_separate_then_join.pl --seqType fq --left TrinityInputRead.filtered_1.fastq --right TrinityInputRead.filtered_2.fastq \
								--target ./trinity_out_dir/Trinity.fasta --aligner bowtie -- -p 4 --all --best --strata -m 300
Trinity_contigs_count=$(grep '>' ./trinity_out_dir/Trinity.fasta | wc -l | awk '{print $1}')
if (($Trinity_contigs_count < 1 )); then
	echo "#Exit: There is no contigs in trinity_out_dir/Trinity.fasta, Exit!" >> $logfile
	exit
fi
# # Blast Trinity.fasta against nt database
# [[ -d blast_out ]] || mkdir blast_out
# blastn_out_format="7 qseqid sseqid sgi staxids slen pident length mismatch gapopen qstart qend sstart send evalue bitsore stitle"
# blastn -db /PHShome/tw786/neurogen/Tao/vfs/annotation/nt/nt -query ./trinity_out_dir/Trinity.fasta -outfmt "$blastn_out_format"  \
# 			-out ./blast_out/Trinity.fasta.blastnout

# # Python Script to parse blast output
# export PATH=~/anaconda_ete/bin:$PATH

# python $script_dir/parseBlast.py $homo_sapiens_taxid $virus_taxid $work_dir/blast_out/Trinity.fasta.blastnout
# if [[ -f $work_dir/blast_out/parseBlast_ouput.txt ]]; then
# 	cat $work_dir/blast_out/parseBlast_ouput.txt
# 	echo "#Successed!" >> $logfile
# 	echo "Successed!"
# else
# 	echo "#failed: parseBlast_ouput.txt is not generated!"
# fi


























