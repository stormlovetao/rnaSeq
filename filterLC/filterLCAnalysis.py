"""
This scprit aims to analysis the Low-Complexity reads filtering results, to see how many good reads are in snap_output.sam

A table will be created, with rows representing Virueses and columns representing samples.


"""

##############################

import os
import sys

import string


##############################

if __name__ == '__main__':
	

	#sample_root_dir = "/PHShome/tw786/neurogen/Tao/run_output/"
	sample_root_dir = "/PHShome/tw786/neurogen/Tao/run_output/"
	species_dir = "/PHShome/tw786/MyOwnScript/filterLC_viruses_list.txt"

	species_list = []
	fp = open(species_dir)
	while True:
		line = fp.readline()
		if not line:
			break
		species_list.append(line.strip())
	fp.close()

	sample_folder_list = os.listdir(sample_root_dir)
	sample_list = []
	sample_dic = {}
	for sample_folder in sample_folder_list:
		if not os.path.exists(sample_root_dir+sample_folder+'/'):
			continue
		if not sample_folder.startswith('PD'):
			continue
		sample_list.append(sample_folder)

		species_dic = {}
		for species in species_list:
			result_dir = sample_root_dir + sample_folder + '/' + species + '/'
			if not os.path.exists(result_dir + 'lc_jobout'):
				print result_dir
				continue

			fp = open(result_dir + 'lc_jobout')
			while True:
				line = fp.readline()
				if not line:
					break
				line = line.lstrip()
				if line.startswith('Good sequences:'):
					linesplit = line.split(':')
					target = linesplit[1].lstrip().split(' ')[0]
					species_dic[species] = target
					continue
		sample_dic[sample_folder] = species_dic


	analysis_out = open("./analysis_out_filter", 'w')
	sample_line = '\t' + '\t'.join(sample_list) + '\n'
	analysis_out.write(sample_line)
	
	for species_name in species_list:
		result_line = []
		for sample_name in sample_list:
			if not sample_dic.has_key(sample_name):
				print sample_name, "don't exists in the sample_dic"
				result_line.append(-1)
				
			else:
				if not sample_dic[sample_name].has_key(species_name):
					result_line.append(-2)
				else:
					result_line.append(sample_dic[sample_name][species_name])
		result_line_join = species_name + '\t' + '\t'.join([item for item in result_line]) + '\n'
		analysis_out.write(result_line_join)



	analysis_out.close()











