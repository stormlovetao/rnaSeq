#!/bin/bash

sample_dir=$1
sample_name=${sample_dir##*/}
if [[ ! -f $sample_dir/kraken.human9606.fusion.txt ]]; then
	python /PHShome/tw786/MyOwnScript/KrakenOutputFindCirRnaFusion.py $sample_name
fi
if [[ ! -f $sample_dir/readsCake.numstat ]]; then
	sh /PHShome/tw786/MyOwnScript/readNumberStat.sh $sample_name
fi

			
if [[ $? -eq 0 ]]; then
	
	[[ -f $sample_dir/unmapped.fastq ]] && rm $sample_dir/unmapped.fastq
	[[ -f $sample_dir/unmapped.filterLC.fastq ]] && rm $sample_dir/unmapped.filterLC.fastq
	[[ -f $sample_dir/unmapped.filterLC.filterPhiX.bam ]] && rm $sample_dir/unmapped.filterLC.filterPhiX.bam
	[[ -f $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.bam ]] && rm $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.bam
	[[ -f $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU.bam ]] && rm $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU.bam
	[[ -f $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.fastq ]] && rm $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.fastq
	[[ -f $sample_dir/AllGeneFusionReads.txt ]] && rm $sample_dir/AllGeneFusionReads.txt
	[[ -f $sample_dir/kraken_output.translate ]] && rm $sample_dir/kraken_output.translate
				
fi