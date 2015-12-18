"""
This scprit aims to analysis the mapping results against whole Viruses genomes.

A table will be created, with rows representing samples and columns representing Viral species.


"""


##############################

import os
import sys

import string


##############################


if __name__ == '__main__':
	

	#sample_root_dir = "/PHShome/tw786/neurogen/Tao/run_output/"
	sample_root_dir = "/PHShome/tw786/neurogen/Tao/snap_run_output/"
	species_root_dir = "/PHShome/tw786/neurogen/Tao/fna/Viral/ftp.ncbi.nlm.nih.gov/genomes/Viruses/"

	sample_root_dir_folder_list = os.listdir(sample_root_dir)

	sample_folder_name_list = []
	species_folder_name_list = []

	species_root_dir_list = os.listdir(species_root_dir)
	for species_name in species_root_dir_list:
		if os.path.exists(species_root_dir + species_name + '/'):
			species_folder_name_list.append(species_name)

	sample_dic = {}

	for sample_folder_name in sample_root_dir_folder_list:

		sample_dir = sample_root_dir + sample_folder_name + '/'
		
		if not os.path.exists(sample_dir):
			print sample_dir, "is not a folder"
			continue

		sample_folder_name_list.append(sample_folder_name)

		species_folder_list = os.listdir(sample_dir)

		species_dic = {}

		for species_folder_name in species_folder_list:
			species_dir = sample_dir + species_folder_name + '/'
			if not os.path.exists(species_dir):
				continue

			joboutfile_path = species_dir + "snap_mapping_jobout"

			if not os.path.exists(joboutfile_path):
				print sample_folder_name, species_folder_name, "Warning: no snap_mapping_jobout file! Check!"
				continue
			jobout = open(joboutfile_path, 'r')
			#print sample_folder_name, species_folder_name
			while True:
				line = jobout.readline()
				if not line:
					break
				if line.startswith("Subject") and line.strip().endswith("Exited"):
					align_rate = -3
					break
				elif line.startswith("Total Reads"):
					
					next_line = jobout.readline()
					next_linesplit = next_line.split('\t')
					align_reads_num_highMAQ = next_linesplit[0].split()[1]
					align_reads_num_lowMAQ = next_linesplit[0].split()[3]
					if ',' in align_reads_num_highMAQ:
						line_split = align_reads_num_highMAQ.split(',')
						align_reads_num_highMAQ = ''.join(line_split)
					if ',' in align_reads_num_lowMAQ:
						line_split = align_reads_num_lowMAQ.split(',')
						align_reads_num_lowMAQ = ''.join(line_split)
					align_reads_num_highMAQ = string.atoi(align_reads_num_highMAQ)
					align_reads_num_lowMAQ = string.atoi(align_reads_num_lowMAQ)
					# align_rate_highMAQ = next_linesplit[0].split()[2]
					# align_rate_lowMAQ = next_linesplit[0].split()[4]
					# align_rate_highMAQ = align_rate_highMAQ[1:len(align_rate_highMAQ)-2]
					# align_rate_lowMAQ = align_rate_lowMAQ[1:len(align_rate_lowMAQ)-2]
					# align_rate_highMAQ = string.atof(align_rate_highMAQ)
					# align_rate_lowMAQ = string.atof(align_rate_lowMAQ)
					# align_rate = align_rate_highMAQ + align_rate_lowMAQ
					align_rate = align_reads_num_highMAQ + align_reads_num_lowMAQ
			species_dic[species_folder_name] = align_rate

		sample_dic[sample_folder_name] = species_dic

	analysis_out = open("./analysis_out_V0.1", 'w')
	sample_line = '\t' + '\t'.join(sample_folder_name_list) + '\n'
	analysis_out.write(sample_line)

	# for sample_name in sample_folder_name_list:
	# 	if not sample_dic.has_key(sample_name):
	# 		print sample_name, "don't exists in the sample_dic"
	# 		continue
	# 	result_line = []
	# 	for species_name in species_folder_name_list:

	# 		if not sample_dic[sample_name].has_key(species_name):
	# 			result_line.append(-1)
	# 		else:
	# 			result_line.append(sample_dic[sample_name][species_name])

	# 	result_line_join = sample_name + '\t' + '\t'.join([str(item) for item in result_line]) + '\n'
	# 	analysis_out.write(result_line_join)

	
	for species_name in species_folder_name_list:
		result_line = []
		for sample_name in sample_folder_name_list:
			if not sample_dic.has_key(sample_name):
				print sample_name, "don't exists in the sample_dic"
				result_line.append(-1)
				
			else:
				if not sample_dic[sample_name].has_key(species_name):
					result_line.append(-2)
				else:
					result_line.append(sample_dic[sample_name][species_name])
		result_line_join = species_name + '\t' + '\t'.join([str(item) for item in result_line]) + '\n'
		analysis_out.write(result_line_join)



	analysis_out.close()










