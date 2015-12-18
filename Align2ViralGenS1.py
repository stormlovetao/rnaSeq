"""
This script aims to align all unmapped samples to Viral/Bacterial genomes.

version: V-0.0.0

1. Download Viral/Bacterial genomes from NCBI(ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/) to local server

2. Use Bowtie/SNAP to build index of reference genomes.

3. Use Bowtie/SNAP to map samples to indexes.

This is step 1: indexing
"""
#######################

import os
import sys
import glob
import string



#######################




if __name__ == '__main__':
	
	# make dirs
	sys_command = "module use /apps/modulefiles/test"
	os.system(sys_command)
	sys_command = "module load bowtie2/2.2.5"
	os.system(sys_command)
	
	fna_dir = "/PHShome/tw786/MyOwnScript/fna/"
	viral_dir = fna_dir + "Viral/"
	if not os.path.exists(viral_dir): 
		os.mkdir(viral_dir)
	# bacterial_dir = fna_dir + "Bacterial/"
	# if not os.path.exists(bacterial_dir):
	# 	os.mkdir(bacterial_dir)

	#wget genomes from NCBI ftp
		#Already download mannully

	all_viral_dirs = viral_dir + "ftp.ncbi.nlm.nih.gov/genomes/Viruses/"
	viral_folder_list = os.listdir(all_viral_dirs)

	for species_name in viral_folder_list:
		species_dir = all_viral_dirs + species_name + '/'
		
		#judge if .fna file exists

		os.chdir(species_dir)
		fna_list = glob.glob('*.fna')

		if len(fna_list) == 0:
			continue
		
		renamed_fna = species_name + ".fna"
		command = "cat " + species_dir + "*.fna > " + renamed_fna
		os.system(command)

		#bowtie2 indexing
		command = "bsub -o jobout bowtie2-build " + species_dir + renamed_fna + ' ' + species_name   # busb email??
		os.system(command)
		break

	##indexing finished