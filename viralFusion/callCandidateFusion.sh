#!/bin/bash

# de nove assemble rna reads labeled to virus and human-virus
# Prepare pair end reads

#extract reads'ID which are mixed by human part and virus part

module_name=`basename $0`
sample_name=HC_NZ-H152_MCPY_4_rep1
homo_sapiens_taxid=9606
virus_taxid=113553
work_dir=/PHShome/tw786/neurogen/Tao/kraken_output/HC_NZ-H152_MCPY_4_rep1
cd $work_dir
awk -v v1=$homo_sapiens_taxid -v v2=$virus_taxid -F '\t' '($5 ~ v1) && ($5 ~ v2)' kraken_output | cut -f 2 > ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID 
sed 's/\/[1-2]*//g' ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID > ${homo_sapiens_taxid}_${virus_taxid}.readsID

#extract reads' ID which are label to virus

awk -v v=$virus_taxid '$3 == v' kraken_output | cut -f 2 | sed 's/\/[1-2]*//g' >> ${homo_sapiens_taxid}_${virus_taxid}.readsID

sort ${homo_sapiens_taxid}_${virus_taxid}.readsID | uniq > ${homo_sapiens_taxid}_${virus_taxid}.readsID1
mv ${homo_sapiens_taxid}_${virus_taxid}.readsID1 ${homo_sapiens_taxid}_${virus_taxid}.readsID

seqtk subseq /data/neurogen/rnaseq_PD/rawfiles/$sample_name.R1.fastq.gz ${homo_sapiens_taxid}_${virus_taxid}.readsID \
				> ${homo_sapiens_taxid}_${virus_taxid}.R1.fq
seqtk subseq /data/neurogen/rnaseq_PD/rawfiles/$sample_name.R2.fastq.gz ${homo_sapiens_taxid}_${virus_taxid}.readsID \
				> ${homo_sapiens_taxid}_${virus_taxid}.R2.fq

#blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 hg38+virus113553.fa 113553-9606.txt.fa output.psl
bsub -q normal -M 6000 Trinity --seqType fq --JM 4G --left 9606_113553.R1.fq --right 9606_113553.R2.fq --CPU 6
bowtie_PE_separate_then_join.pl --seqType fq --left 9606_113553.R1.fq --right 9606_113553.R2.fq --target ./trinity_out_dir/Trinity.fasta --aligner bowtie -- -p 4 --all --best --strata -m 300