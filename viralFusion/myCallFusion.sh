#!/bin/bash

# Call virus and human fusion 
# step1, find candidate fusion reads using kraken's results
# step2, blast those candidate reads against nt database to check if the reads is specificly mapped to human and virus
# step3, blat reads outputed from step2 against a new reference which is combined by hg38 and virus genome, then find the exact fusion position and fusion pattern.

module_name=`basename $0`

sample_name=HC_NZ-H152_MCPY_4_rep1
homo_sapiens_taxid=9606
virus_taxid=113553
work_dir=/PHShome/tw786/neurogen/Tao/kraken_output/HC_NZ-H152_MCPY_4_rep1
cd $work_dir

# extract reads' ID which are labeled with $homo_sapiens_taxid and $virus_taxid
awk -v v1=$homo_sapiens_taxid -v v2=$virus_taxid -F '\t' '($5 ~ v1) && ($5 ~ v2)' kraken_output | cut -f 2 > ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID 
# remove R1/R2 tag from reads' ID
# sed 's/\/[1-2]*//g' ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID > ${homo_sapiens_taxid}_${virus_taxid}.readsID

# creat a fasta file using mixed reads' ID from source fq file (here, we extract from bam file)
samtools view unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam | grep -w -f ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID \
			| awk -F '\t' '{print ">"$1"\n"$10}' > ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID.fa
blastn_out_format="7 qseqid sseqid sgi staxids slen pident length mismatch gapopen qstart qend sstart send evalue bitsore stitle"
blastn -db /PHShome/tw786/neurogen/Tao/vfs/annotation/nt/nt -query ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID.fa -outfmt $blastn_out_format      \
			-out ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID.fa.blastnout
			"7 qseqid sseqid sgi staxids slen pident length mismatch gapopen qstart qend sstart send evalue bitsore stitle"
# blastn -db /PHShome/tw786/neurogen/Tao/vfs/annotation/nt/nt -query test1.fa -outfmt "7 qseqid sseqid sgi staxids slen pident length mismatch gapopen qstart qend sstart send evalue bitsore stitle" -out test1.fa.blastnout

# Python script to check if a read is specificly mixed by human part and virus part
# input output of blastn, ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID.fa.blastnout
# output  ${homo_sapiens_taxid}_${virus_taxid}.specificMixedReads.fa
# enable ete3 module
export PATH=~/anaconda_ete/bin:$PATH
python checkSpecificity.py $homo_sapiens_taxid $virus_taxid $work_dir/${homo_sapiens_taxid}_${virus_taxid}.specificMixedReads.fa
# need to modify the interface

# Creat blat database list(blatdblist.tmp) containing hg38 and virus genome
hg38=/PHShome/tw786/neurogen/Tao/hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna
virus_fa=/PHShome/tw786/neurogen/Tao/fna/KrakenSDBViral/virus_kraken_sdb.fa
virus_fa_header =`grep $virus_taxid $gi2taxid | cut -f 1`
chrom_id=`echo $virus_fa_header | sed 's/\(.*|\).*/\1/' | sed 's/>//' `
samtools faidx $virus_fa $chrom_id > virus_gi${virus_giid}taxid${virus_taxid}.fa
echo "generated virus_gi${virus_giid}taxid${virus_taxid}.fa"
echo $hg38 > blatdblist.tmp
echo $work_dir/virus_gi${virus_giid}taxid${virus_taxid}.fa >> blatdblist.tmp

# Python script gi2taxid.py to map virus gi to taxid
# input: /PHShome/tw786/neurogen/Tao/kraken_standard_db/taxonomy/gi_taxid_nucl.dmp
# 		/virus_fa=/PHShome/tw786/neurogen/Tao/fna/KrakenSDBViral/virus_kraken_sdb.fa
# output: /PHShome/tw786/neurogen/Tao/fna/KrakenSDBViral/virus_kraken_sdb.fa.gi2taxid
gi2taxid=/PHShome/tw786/neurogen/Tao/fna/KrakenSDBViral/virus_kraken_sdb.fa.virus2taxid


blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 blatdblist.tmp ${homo_sapiens_taxid}_${virus_taxid}.specificMixedReads.fa \
		${homo_sapiens_taxid}_${virus_taxid}.specificMixedReads.fa.blatout.psl

# Python script to analyze output of blat
























