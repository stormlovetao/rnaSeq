#!/bin/bash

#################################################
# called by filterLC.sh
#################################################

module load tophat/default
module load samtools/default

samtools view -b -S $1 > $2
bam2fastx -q -Q -A -N -o $3 $2
perl $HOME/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $3 -out_format 3 -out_good $4 -out_bad $5 -log -lc_method dust -lc_threshold 7 
