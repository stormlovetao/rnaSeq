#!/bin/bash
prog=$0
sample_dir=$1
prog_dir=/PHShome/tw786/MyOwnScript
source_fq=unmapped_qc.filterPhiX.fastq.gz

NT=/PHShome/tw786/neurogen/Tao/vfs/annotation/nt/nt
cd $sample_dir
#abbr. M, Micribe  CF, ConFirmed
function error_exit
{
	echo "$1" 1>&2
	exit 1
}
if [[ ! -f kraken_output ]] ; then
	error_exit "not find $sample_dir/kraken_output"
fi

echo [`date`] "Program: $prog" >> rerun.log
awk '$1 == "C" && $3 != 9606' kraken_output  > kraken_output_M
awk '$4 == "S"' kraken_output.report | cut -f 5 > kraken_report_S.txt
python $prog_dir/kraken_output_extract.py  kraken_output_M kraken_report_S.txt > kraken_output_M_S

cut -f 2 kraken_output_M_S | seqtk subseq $source_fq - | seqtk seq -A - > kraken_output_M_S.fa

echo [`date`] "Blastn Start" >> rerun.log
blastn -task 'megablast' -db $NT -query kraken_output_M_S.fa -num_threads 16 -evalue 1.0e-10 -max_target_seqs 20 \
-outfmt '6 qseqid sseqid sgi staxids slen pident length mismatch gapopen qstart qend sstart send evalue bitsore'  \
-out kraken_output_M_S.blastout  

if [[ $? -ne 0 ]]; then
	error_exit "Blastn error! Exit!"
else
	echo [`date`] "Blastn Completed!" >> rerun.log
fi
[[ -f kraken_output_M_S.fa.gz ]]&& rm kraken_output_M_S.fa.gz
gzip kraken_output_M_S.fa