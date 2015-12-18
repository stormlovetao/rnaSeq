#!/bin/bash
# bsub -q big  -R 'rusage[mem=8500]'       \
# -oo /PHShome/tw786/neurogen/Tao/run_output/PD_BN04-42_SNDA_5_rep1/Shamonda_virus_uid173358/kraken_jobout \
# $HOME/bin/kraken/kraken      \
# --db  /PHShome/tw786/minikraken_20141208           \
# --fastq-input  /PHShome/tw786/neurogen/Tao/run_output/PD_BN04-42_SNDA_5_rep1/Shamonda_virus_uid173358/snap_output.good.fastq \
# > /PHShome/tw786/neurogen/Tao/run_output/PD_BN04-42_SNDA_5_rep1/Shamonda_virus_uid173358/kraken_output


~/bin/kraken/kraken-translate --db ~/minikraken_20141208 ./kraken_jobout.grep > ./kraken_jobout.grep.translate

bsub -q big -R 'rusage[mem=80000]' \
-oo /PHShome/tw786/neurogen/Tao/run_output/PD_BN04-42_SNDA_5_rep1/Shamonda_virus_uid173358/kraken_standard_db_jobout \
$HOME/bin/kraken/kraken      \
--db  /PHShome/tw786/neurogen/Tao/kraken_standard_db           \
--fastq-input  /PHShome/tw786/neurogen/Tao/run_output/PD_BN04-42_SNDA_5_rep1/Shamonda_virus_uid173358/snap_output.good.fastq \

~/bin/kraken/kraken-translate --db /PHShome/tw786/neurogen/Tao/kraken_standard_db ./kraken_standard_db_jobout.grep > ./kraken_standard_db_jobout.grep.translate

samples_root_dir=/data/neurogen/rnaseq_PD/run_output 
for sample_dir in $samples_root_dir/*
do
	if [[ -d sample_dir ]]; then
		sample_name=${sample_dir##*/}
		if [[ $sample_name == "PD"* ]]; then
			if [[ ${sample_name:(-6):1} -ne "1" ]]; then

			fi
		fi
	fi
