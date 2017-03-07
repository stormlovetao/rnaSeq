#!/bin/bash
# bsub -q big-multi -n 4/8 -M 120000 -R 'rusage[mem=120000]' 
prog=$0
tissue=$1
GTExRootdir=/PHShome/tw786/neurogen/Tao/GTEx_output/$tissue
if [[ ! -d $GTExRootdir ]]; then
	echo "Not found:" $GTExRootdir
fi
kraken_standard_db=/PHShome/tw786/neurogen/Tao/Kraken_DB
kraken_bin=kraken
kraken_report_bin=kraken-report
for sample_dir in $GTExRootdir/*
do
	if [[ -d $sample_dir  ]]; then
		sample_name=${sample_dir##*/}

		cd $sample_dir
		if [[ -f unmapped_qc.fa.gz ]]; then
			if [[ ! -f kraken_output ]]; then
				if [[ -f _kraken_running ]] || [[ -f _kraken_finished ]]; then
					continue
				fi
				touch _kraken_running
				echo "Program: $prog" > rerun.log
				echo "kraken DB: $kraken_standard_db" >> rerun.log
				echo [`date`] "Start" >> rerun.log
				$kraken_bin --threads 16 --db $kraken_standard_db \
				--fasta-input unmapped_qc.fa.gz \
				--gzip-compressed > kraken_output
				if [[ $? -eq 0 ]]; then
					echo [`date`] "Kraken completed!" >>  rerun.log
				else
					echo $sample_dir "Kraken failed!"
					[[ -f _kraken_running ]]&& rm _kraken_running
					exit 1
				fi
				$kraken_report_bin --db $kraken_standard_db kraken_output > kraken_output.report
				if [[ $? -eq 0 ]]; then
					echo [`date`] "Kraken-report completed!" >>  rerun.log
					[[ -f _kraken_running ]]&& mv _kraken_running _kraken_finished
				else
					[[ -f _kraken_running ]]&& rm _kraken_running
					echo $sample_dir "Kraken-report failed!"
					exit 1
				fi
				#bsub -q big -M 30000 -R 'rusage[mem=30000]' sh BC_rerun_blast.sh $sample_dir
			fi
		fi
	fi
done



