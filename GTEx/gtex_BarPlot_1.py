"""
Output the raw table of GTEx
Make a landscape barplot of viruses detected from all GTEx brain samples
(Report April 26 2016)
reads_num >= 10 or p.p.m > 2
"""

import os
import sys
import string
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

brain_regions = ['Amygdala', 'AnteriorCingulateCortex', 'Caudate', 'CerebellarHemisphere', 'Cerebellum', 'Cortex', 'FrontalCortex', 'Hippocampus', 'Hypothalamus', 'NucleusAccumbens', 'Putamen', 'SpinalCord', 'SubstantialNigra', 'Liver']
brain_regions.sort()
level_mark = 'S'
root_dir = '/PHShome/tw786/neurogen/Tao/GTEx_output/'
virus_read_thr = 10
ppmI_thr = 2
ppmII_thr = 1


goodseq_num_dic = {}
total_library_dic = {}
all_sample_summary_path = '/PHShome/tw786/neurogen/Tao/GTEx_output/GTEx_allSamples_reads_statistics.txt'
total_fp = open(all_sample_summary_path)
header = total_fp.readline()
while True:
	line = total_fp.readline()
	if not line:
		break
	linesplit = line.strip().split('\t')
	srr_id = linesplit[1]
	MBase = string.atoi(linesplit[4])
	library_reads_num = MBase * 100000/76
	total_library_dic[srr_id] = library_reads_num
	goodseq_num_dic[srr_id] = string.atoi(linesplit[2])

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
			
				#if viruses_start and reads_num >= virus_read_thr:
				ppmI = float(reads_num*1000000)/goodseq_num_dic[sample_name]
				#ppmII = float(reads_num*1000000)/total_library_dic[sample_name]
				if viruses_start and  ppmI > ppmI_thr:
				#if viruses_start and ppmII > ppmII_thr:
					if 'phage' in name:
						continue
					if 'Geobacillus virus E2' == name:
						continue
					if 'Choristoneura occidentalis granulovirus' == name:
						continue
					if not virus_count.has_key(name):
						virus_count[name] = {}
						virus_count[name][sample_name] = reads_num
						#virus_count[name] = [(sample_name,reads_num)]
					else:
						virus_count[name][sample_name] = reads_num
						#virus_count[name].append((sample_name,reads_num))
		kraken_output_fp.close()
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


virus_sampleNum_tuples = sorted(virus_sampleNum_dic.items(), key = lambda x:x[1], reverse = True)
virus_list = [x[0] for x in virus_sampleNum_tuples]
sample_num_list = [x[1] for x in virus_sampleNum_tuples]
########################################Plot####################################
fig = plt.figure(figsize =(10,12),tight_layout=True)
ind = np.arange(len(virus_list))
width = 0.35
ax1 = fig.add_subplot(2,1,1)
plt.bar(ind, sample_num_list, width=width, color='r')
plt.xticks(ind+width/2., virus_list,rotation=270, fontsize=5)
plt.xlim(0,len(virus_list)+1)
	# plt.axhline(y = 1, color = 'b')
plt.ylabel('Number of samples')
plt.title("#Samples detected virus (GTEx-BrainLiver, ppmI > %.1f)" % (ppmI_thr))

#########################################Repeat####################################
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
			
				#if viruses_start and reads_num >= virus_read_thr:
				#ppmI = float(reads_num*1000000)/goodseq_num_dic[sample_name]
				ppmII = float(reads_num*1000000)/total_library_dic[sample_name]
				#if viruses_start and  ppmI > ppmI_thr:
				if viruses_start and ppmII > ppmII_thr:
					if 'phage' in name:
						continue
					if 'Geobacillus virus E2' == name:
						continue
					if 'Choristoneura occidentalis granulovirus' == name:
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
virus_sampleNum_tuples = sorted(virus_sampleNum_dic.items(), key = lambda x:x[1], reverse = True)
virus_list = [x[0] for x in virus_sampleNum_tuples]
sample_num_list = [x[1] for x in virus_sampleNum_tuples]

ind = np.arange(len(virus_list))
width = 0.35
ax2 = fig.add_subplot(2,1,2)
plt.bar(ind, sample_num_list, width=width, color='r')
plt.xticks(ind+width/2., virus_list,rotation=270, fontsize=5)
plt.xlim(0,len(virus_list)+1)
	# plt.axhline(y = 1, color = 'b')
plt.ylabel('Number of samples')
plt.title("#Samples detected virus (GTEx-BrainLiver,ppmII>%.1f)" %(ppmII_thr) )









save_folder = "/PHShome/tw786/localView/overview"
if not os.path.isdir(save_folder):
	os.mkdir(save_folder)
savepath = save_folder + "/GTExBarplot1_AllViruses_across_AllSamples_ppmI_II.pdf"
plt.savefig(savepath)




