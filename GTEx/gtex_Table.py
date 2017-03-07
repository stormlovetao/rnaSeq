"""
output a table
each row represents a virus
each column represents a sample clustered by brain region
"""
import os
import sys
import string

brain_regions = ['Amygdala', 'AnteriorCingulateCortex', 'Caudate', 'CerebellarHemisphere', 'Cerebellum', 'Cortex', 'FrontalCortex', 'Hippocampus', 'Hypothalamus', 'NucleusAccumbens', 'Putamen', 'SpinalCord', 'SubstantialNigra']
brain_regions.sort()
level_mark = 'S'
root_dir = '/PHShome/tw786/neurogen/Tao/GTEx_output/'
virus_read_thr = 10

brain_regions_virus_count_dic = {}
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
				continue

			if linesplit[4] == '10239': #Virus(taxonimic id = 10239)
				eukaryota_start = False
				viruses_start = True
				bacteria_start = False
				archaea_start = False
				continue
			if linesplit[4] == '2157': #archaea(taxonimic id = 2157)
				eukaryota_start = False
				viruses_start = False
				bacteria_start = False
				archaea_start = True
				continue
			if linesplit[4] == '2759': #eukaryota(taxonimic id = 2759)
				eukaryota_start = True
				viruses_start = False
				bacteria_start = False
				archaea_start = False
				continue

			if linesplit[3] == level_mark:
				name = linesplit[5].lstrip().strip()
				
				reads_num = string.atoi(linesplit[1].strip())
				species_taxid = linesplit[4]
			
				if viruses_start and reads_num >= virus_read_thr:
					if 'phage' in name:
						continue
					if not virus_count.has_key(name):
						virus_count[name] = {}
						virus_count[name][sample_name] = reads_num
						#virus_count[name] = [(sample_name,reads_num)]
					else:
						virus_count[name][sample_name] = reads_num
						#virus_count[name].append((sample_name,reads_num))
	brain_regions_virus_count_dic[brain_region] = virus_count

# output a sorted table

# sort rows/viruses by the number of samples they are detected(above threshold)
virus_sampleNum_dic = {}
for brain_region in brain_regions_virus_count_dic:
	virus_count = brain_regions_virus_count_dic[brain_region]
	for virus in virus_count:
		sampleNum = len(virus_count[virus])
		if virus_sampleNum_dic.has_key(virus):
			virus_sampleNum_dic[virus] += sampleNum
		else:
			virus_sampleNum_dic[virus] = sampleNum

virus_index = sorted(virus_sampleNum_dic.items(), key = lambda x:x[1], reverse = True)
virus_index = [x[0] for x in virus_index]

# sort columns/samples by number of virus detected in each sample/column
sample_index_dic = {}

for brain_region in brain_regions_virus_count_dic:
	sample_virus_dic = {}
	virus_count = brain_regions_virus_count_dic[brain_region]
	for virus in virus_count:
		sample_number_dic = virus_count[virus]
		for sample_name in sample_number_dic:
			if sample_virus_dic.has_key(sample_name):
				sample_virus_dic[sample_name].append(virus)
			else:
				sample_virus_dic[sample_name] = [virus]
	sample_index = sorted(sample_virus_dic.items(), key = lambda x:len(x[1]), reverse = True)
	sample_index = [x[0] for x in sample_index]
	sample_index_dic[brain_region] = sample_index

# output the table

output_path = "/PHShome/tw786/localView/GTEx_allVirus_allSamples_Table/GTEx_allVirus_allSamples_Table.xls"
out_fp = open(output_path, 'w')

header_line = 'BrainRegion\t' 
for brain_region in sample_index_dic:
	sample_number_in_this_region = len(sample_index_dic[brain_region])
	brain_region_seq = [brain_region]*sample_number_in_this_region
	header_line = header_line +'\t'.join(brain_region_seq) + '\t'
header_line = header_line.strip()
header_line += '\n'
out_fp.write(header_line)

sample_line = 'SampleName\t'
for brain_region in sample_index_dic:
	sample_line = sample_line +'\t'.join(sample_index_dic[brain_region]) + '\t'
sample_line = sample_line.strip()
sample_line += '\n'
out_fp.write(sample_line)

for virus in virus_index:
	outline = virus + '\t'
	for brain_region in sample_index_dic:
		sample_index = sample_index_dic[brain_region]
		for sample in sample_index:
			if not brain_regions_virus_count_dic[brain_region].has_key(virus):
				reads_num = 0
			elif brain_regions_virus_count_dic[brain_region][virus].has_key(sample):
				reads_num = brain_regions_virus_count_dic[brain_region][virus][sample]
			else:
				reads_num = 0
			outline += str(reads_num) + '\t'
	outline = outline.strip()
	outline += '\n'
	out_fp.write(outline)













					