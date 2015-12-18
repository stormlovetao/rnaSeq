"""
I submitted all the mapping jobs(map all samples against whole Viral genomes) to "short" queque on ERISOne

But some jobs failed because of swap usage limit
(e.g TERM_SWAP: job killed after reaching LSF swap usage limit.
Exited with exit code 130.)

So, this script aims to filter out the failed jobs

"""





##############################

import os
import sys

import string
import subprocess


##############################

run_output_dir = "/PHShome/tw786/neurogen/Tao/run_output/"
index_dir = "/data/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/"
samples_root_dir="/data/neurogen/rnaseq_PD/run_output/"

#failed_jobs = open("/PHShome/tw786/neurogen/Tao/run_output/failed_jobs.txt",'w')

if __name__ == '__main__':

	command = "module use /apps/modulefiles/test"
	subprocess.call(command, shell = True)
	command = "module load snap/1.0.20"
	subprocess.call(command, shell = True)

	
	samples_list = os.listdir(run_output_dir)
	for sample_folder_name in samples_list:
		sample_dir = run_output_dir + sample_folder_name + '/'
		if not os.path.exists(sample_dir):
			continue

		count = 0

		species_list = os.listdir(sample_dir)
		for species_folder_name in species_list:
			species_dir = sample_dir + species_folder_name + '/'
			if not os.path.exists(species_dir):
				continue
			flag = 0
			if not os.path.exists(species_dir + "snap_mapping_jobout"):
				#print sample_folder_name, species_folder_name, "don't have snap_mapping_jobout"
				flag = 1
			else:
				jobout = open(species_dir + "snap_mapping_jobout")
				
				while True:
					line = jobout.readline()
					if not line:
						break
					if line.startswith("Subject"):
						if line.strip().endswith("Done"):
							count += 1
							flag = 1
							# outline = sample_folder_name + '\t' + species_folder_name + '\n'
							# failed_jobs.write(outline)
							break
				jobout.close()
			# if flag == 1:
			# 	command = "bsub -q 'big'  -J " + sample_folder_name+"-mapping-"+species_folder_name + " -oo " + \
			# 	species_dir+"snap_mapping_jobout snap-aligner  single " + index_dir + species_folder_name + "/snap_index"\
			# 	+" -bam " + samples_root_dir + sample_folder_name + "/unmapped.bam -o " + species_dir+"snap_output.sam -F a > " \
			# 	+species_dir+"stdout.txt"

			# 	os.system(command)
				# print sample_folder_name, species_folder_name
				# exit(1)
				
		print sample_folder_name, count

