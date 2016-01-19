#######################
#Step6: Filter reads that support circRna and gene-fusion
#######################


if [[ ! -f $output_sample_dir/AllCircRnaList.txt ]]; then
	python /PHShome/tw786/MyOwnScript/getCircRna.py $input_sample_name
	if [[ $? -ne 0 ]]; then
	 	error_exit "ERROR exists in getCircRna.py script"
	fi 
		
fi

if [[ ! -f $output_sample_dir/AllGeneFusionReads.txt ]]; then
	gene_fusion_bamfile=/PHShome/tw786/neurogen/rnaseq_PD/run_output/$input_sample_name/tophat_fusion/accepted_hits.bam
	samtools view $gene_fusion_bamfile | cut -f 1 > $output_sample_dir/AllGeneFusionReads.txt
	if [[ $? -ne 0 ]]; then
		error_exit "ERROR exists in get reads from tophat_fusion"
	fi
fi





if [[ ! -f $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.filterCF.bam ]]; then
	samtools view -h $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam > \
					$output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.sam

	python filterCircRnaGeneFusion.py $input_sample_name


	if [[ -f $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.filterCF.sam ]]; then
		samtools view -b -S $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.filterCF.sam > \
							$output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.filterCF.bam
	fi
	[[ -f $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.sam ]] && rm $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.sam
	[[ -f $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.filterCF.sam ]] && rm $output_sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.filterCF.sam 
	echo "Reads support circRna and gene-fusion are filtered!"
fi