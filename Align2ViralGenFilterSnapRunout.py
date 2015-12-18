"""

This scripts aims to filter out the snap runout to generate a table for each sample.
Each line will represent a viruses, and each column represents as follows:

Total Reads    Aligned, MAPQ >= 10    Aligned, MAPQ < 10     Unaligned              Too Short/Too Many Ns  %Pairs	Reads/s   Time in Aligner (s)
38,831,844     9,878 (0.03%)          55 (0.00%)             38,014,993 (97.90%)    806,918 (2.08%)        0.00%%	836,857   46

"""
import os

if __name__ == '__main__':
	
	sample_root_dir = "/PHShome/tw786/neurogen/Tao/run_output/"
	species_root_dir = "/PHShome/tw786/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/"

	sample_root_dir_folder_list = os.listdir(sample_root_dir)
	header_line = "\tTotal Reads    Aligned, MAPQ >= 10    Aligned, MAPQ < 10     Unaligned              Too Short/Too Many Ns  %Pairs	Reads/s   Time in Aligner (s)\n"

	for sample_folder_name in sample_root_dir_folder_list:
		sample_dir = sample_root_dir + sample_folder_name + '/'
		
		if not os.path.exists(sample_dir):
			print sample_dir, "is not a folder"
			continue

		outfile_dir = sample_dir + 'snap_runout_analysis.txt'
		outfile_fp = open(outfile_dir , 'w')
		outfile_fp.write(header_line)

		species_folder_list = os.listdir(sample_dir)
		for species_folder_name in species_folder_list:
			if not os.path.exists(sample_dir + species_folder_name + '/'):
				continue
			runout_dir = sample_dir + species_folder_name + '/snap_mapping_jobout'
			if not os.path.exists(runout_dir):
				outline = species_folder_name + "\tsnap_mapping_jobout Not exists\n"
				outfile_fp.write(outline)
				continue
			runout_fp = open(runout_dir)
			while True:
				line = runout_fp.readline()
				if not line:
					break
				if line.startswith("Total Reads"):
					line = runout_fp.readline()
					outline = species_folder_name + '\t' + line
					outfile_fp.write(outline)
		outfile_fp.close()








