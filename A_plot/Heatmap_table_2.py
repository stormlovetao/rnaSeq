"""
Reference: Heatmap_Table_1.py
BrainCode data
output viruses reads number >= 1
Output a table of viruses found in BrainCode
Rows are viruses
Columns are samples, labeled by sample_id.
"""

import os
import sys
import string


brain_regions = ['BRAINCODE_output']

level_mark = 'S'
root_dir = '/PHShome/tw786/neurogen/Tao/'

virus_read_thr = 1

############################# READ viruses reads FOR GTEX - EACH BRAIN REGION and Liver #############################
brain_regions_virus_count_dic = {}
sample_contain_virusnum = {}
for brain_region in brain_regions:
	virus_count = {}
	samples_root_dir = root_dir + brain_region + '/'
	if not os.path.isdir(samples_root_dir):
		print samples_root_dir, 'does not exist!'
		sys.exit()
	for sample_name in os.listdir(samples_root_dir):
		if not os.path.isdir(samples_root_dir + sample_name):
			continue
		
		kraken_output_path = samples_root_dir + sample_name + '/kraken_output.report'

		if not os.path.isfile(kraken_output_path):
			continue
		
		kraken_output_fp = open(kraken_output_path,'r')

		eukaryota_start = False
		viruses_start = False
		bacteria_start = False
		archaea_start = False
		sample_contain_virusnum[sample_name] = 0
		while True:
			line = kraken_output_fp.readline()
			if not line:
				break
			linesplit = line.split('\t')

			if linesplit[4] == '2': # taxonimic id of bacteria is 2 
				eukaryota_start = False
				viruses_start = False
				bacteria_start = True
				archaea_start = False
				fugi_start = False
				continue

			if linesplit[4] == '10239': #Virus(taxonimic id = 10239)
				eukaryota_start = False
				viruses_start = True
				bacteria_start = False
				archaea_start = False
				fugi_start = False
				continue
			if linesplit[4] == '2157': #archaea(taxonimic id = 2157)
				eukaryota_start = False
				viruses_start = False
				bacteria_start = False
				archaea_start = True
				fugi_start = False
				continue
			if linesplit[4] == '2759': #eukaryota(taxonimic id = 2759)
				eukaryota_start = True
				viruses_start = False
				bacteria_start = False
				archaea_start = False
				fugi_start = False
				continue
				
			if linesplit[4] == '4751': #Fugi(taxonimic id = 4751)
				eukaryota_start = False
				viruses_start = False
				bacteria_start = False
				archaea_start = False
				fugi_start = True
				continue

			if linesplit[3] == level_mark:
				name = linesplit[5].lstrip().strip()
				
				reads_num = string.atoi(linesplit[1].strip())
				species_taxid = linesplit[4]
			
				# ppm = round(float(reads_num*1000000)/goodseq_num_dic[sample_name], 3)
				#if viruses_start and  ppm > ppm_thr:
				if viruses_start and  reads_num >= virus_read_thr:
					if 'phage' in name:
						continue
					if 'Geobacillus virus E2' == name:
						continue
					if 'Choristoneura occidentalis granulovirus' == name:
						continue
					if not virus_count.has_key(name):
						virus_count[name] = {}
						virus_count[name][sample_name] = reads_num
					else:
						virus_count[name][sample_name] = reads_num
					sample_contain_virusnum[sample_name] += 1
	brain_regions_virus_count_dic[brain_region] = virus_count
# brain_regions_virus_count_dic[brain_region][virus][sample_name]=ppm

braincode_virus_count = brain_regions_virus_count_dic['BRAINCODE_output']
# braincode_virus_count[virus_name][sample_name]
############################# MAKE TABLE #############################
# Make a table
# Rows represnt viruses
# Columns represent samples
table_out_path = '/PHShome/tw786/localView/BRAINCODETables/braincode_table_samples_reads_%d.xls' % virus_read_thr
out_fp = open(table_out_path, 'w')


# Define row index: sort rows/viruses by the number of samples they are detected(above threshold)
virus_sampleNum_dic = {}
for virus in braincode_virus_count:
	sampleNum = len(virus_count[virus])
	if virus_sampleNum_dic.has_key(virus):
		virus_sampleNum_dic[virus] += sampleNum
	else:
		virus_sampleNum_dic[virus] = sampleNum
virus_index = sorted(virus_sampleNum_dic.items(), key = lambda x:x[1], reverse = True)
virus_index = [x[0] for x in virus_index]

# Define column index: sort columns/samples 
sample_names_list = []
for virus in braincode_virus_count:
	sample_names_list.extend(braincode_virus_count[virus])
sample_names_list = list(set(sample_names_list))
sample_names_list.sort()

# Write header line1: sample name
header_line = 'SampleName:\t'
header_line += '\t'.join(sample_names_list) + '\n'
out_fp.write(header_line)

# Write table
for virus in virus_index:
	outline = virus + '\t'
	out_seq = []
	for sample_name in sample_names_list:
		reads_num = 0
		if braincode_virus_count.has_key(virus):
			if braincode_virus_count[virus].has_key(sample_name):
				reads_num = braincode_virus_count[virus][sample_name]
		out_seq.append(reads_num)
	out_seq_str = [str(x) for x in out_seq]
	outline += '\t'.join(out_seq_str) + '\n'
	out_fp.write(outline)

out_fp.close()










