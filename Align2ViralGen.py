"""
This script aims to align all unmapped samples to Viral/Bacterial genomes.

version: V-0.0.0

1. Download Viral/Bacterial genomes from NCBI(ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/) to local server

2. Use Bowtie/SNAP to build index of reference genomes.

3. Use Bowtie/SNAP to map samples to indexes.

"""
#######################

import os
import sys
import string

sys_command = "module use /apps/modulefiles/test"
os.system(sys_command)
sys_command = "module load bowtie2/2.2.5"
os.system(sys_command)

#######################

def bowtie2_mapping_core(index_path, sample_path, outsam_path, log_path):
	command = "bsub bowtie2 -x " + index_path + ' -U ' + sample_path + ' --no-unal -S ' +  sam_path + ' 2> ' + log_path
	os.system(command)



if __name__ == '__main__':
	
	# make dirs

	fna_dir = "/PHShome/tw786/MyOwnScript/fna/"
	viral_dir = fna_dir + "Viral/"
	if not os.path.exists(viral_dir):
		os.mkdir(viral_dir)
	bacterial_dir = fna_dir + "Bacterial/"
	if not os.path.exists(bacterial_dir):
		os.mkdir(bacterial_dir)

	#wget genomes from NCBI ftp
		#Already download mannully

	all_viral_dirs = viral_dir + "ftp.ncbi.nlm.nih.gov/genomes/Viruses/"
	viral_folder_list = os.listdir(all_viral_dirs)

	for species_name in viral_folder_list:
		species_dir = all_viral_dirs + species_name + '/'
		#fna_list = os.listdir(species_dir)
		
		renamed_fna = species_name + ".fna"
		command = "cat " + species_dir + "*.fna > " + renamed_fna
		os.system(command)

		#bowtie2 indexing

		os.chdir(species_dir)

		command = "bsub bowtie2-build " + species_dir + renamed_fna + ' ' + species_name   # busb email??
		os.system(command)

	##indexing finished

	############# Maybe the next part -mapping- should be in the other script ###############

	samples_root_dir = "/PHShome/tw786/neurogen/rnaseq_PD/run_output/" ### NEED to transform bam to fastq first!!!
	
	samples_folder_name_list = os.listdir(samples_root_dir)

	for species_name in viral_folder_list:
		species_dir = all_viral_dirs + species_name + '/'

		for sample_folder_name in samples_folder_name_list:
			sample = sample_folder_name + "_unmapped.fastq"
			sample_path = samples_root_dir + sample_folder_name + '/' + sample
			os.chdir(species_dir)

			index_path = all_viral_dirs + species_name + '/' + species_name

			outsam_path = samples_root_dir + sample_folder_name + '/' + species_name + '.sam'
			log_path = outsam_path + '.log'

			bowtie2_mapping_core(index_path, sample_path, outsam_path, log_path)
			














		





