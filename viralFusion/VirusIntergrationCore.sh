#!/bin/bash

# Usage: sh VirusIntergrationCore.sh Sample_name Virus_name
# Use VirusFinder2.0


usage()
{
	echo "usage: `basename $0` [sample_name] [virus_name, e.g. 'Murine leukemia virus']"
}

######################
if [[ $# != 2 ]]; then
	usage
	exit
fi

sample_name=$1
virus_name=$2

# Check if sample_name and virus_name exist
kraken_output_dir=/PHShome/tw786/neurogen/Tao/kraken_output
sample_dir=$kraken_output_dir/$sample_name
if [[ ! -d $sample_dir ]]; then
	echo "check the sample_name you input!", $sample_name
	usage
	exit
fi
virus_fa=/PHShome/tw786/neurogen/Tao/fna/KrakenSDBViral/virus_kraken_sdb.fa

virus_gi_header=`grep "$virus_name" $virus_fa`

if [[ $virus_gi_header == \>gi* ]]; then
	echo "find this virus: " $virus_gi_header
else
	echo "check the virus_name you input!" $virus_name
	usage
	exit
fi

# make folder 
output_dir=/PHShome/tw786/neurogen/Tao/VirusFinder/$sample_name
[[ -d $output_dir ]] || mkdir $output_dir

new_virus_name=`echo $virus_name | sed 's/ /_/g'`
output_virus_dir=$output_dir/$new_virus_name
[[ -d $output_virus_dir ]] || mkdir $output_virus_dir

# prepare virus_genome sequence and config file for VirusFinder2.0
chrom_id=`echo $virus_gi_header | sed 's/\(.*|\).*/\1/' | sed 's/>//' `

samtools faidx $virus_fa $chrom_id > $output_virus_dir/${new_virus_name}.fa

echo "generated ${new_virus_name}.fa"

basic_config=/PHShome/tw786/VirusFinder2.0/basic_config.txt

sed "s/SAMPLE_NAME/$sample_name/g" $basic_config > $output_virus_dir/VirusFinderConfig.txt
echo "generated VirusFinderConfig.txt"

#prepare R1 R2 input fastq reads
if [[ ! -f $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.R1.fq ]]; then
	if [[ ! -f $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam ]]; then
		echo $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam "don't exist, exit!"
		exit 1
	fi
	samtools view $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam | cut -f 1 | sed 's/\/[1-2]*//g' | sort | uniq \
				> $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam.uniqReadsID

	if [[ ! -f /data/neurogen/rnaseq_PD/rawfiles/$sample_name.R1.fastq.gz ]]; then
		echo /data/neurogen/rnaseq_PD/rawfiles/$sample_name.R1.fastq.gz "don't exist, exit!"
		exit 1
	fi

	seqtk subseq /data/neurogen/rnaseq_PD/rawfiles/$sample_name.R1.fastq.gz $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam.uniqReadsID \
				> $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.R1.fq

	seqtk subseq /data/neurogen/rnaseq_PD/rawfiles/$sample_name.R2.fastq.gz $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.bam.uniqReadsID \
				> $sample_dir/unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.R2.fq

	echo "generated R1.fq and R2.fq"
fi

# use VirusFinder2.0 to detect virus fusion
VirusFinderScript=/PHShome/tw786/VirusFinder2.0/VirusFinder.pl
config_file=$output_virus_dir/VirusFinderConfig.txt
input_virus_fa=$output_virus_dir/$new_virus_name.fa
perl $VirusFinderScript -c $config_file -o $output_virus_dir -v $input_virus_fa











