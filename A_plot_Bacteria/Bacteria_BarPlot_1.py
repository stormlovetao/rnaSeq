"""
8 Braincode barplot side by side
SNDA(107 samples),TCPY(10 samples),PBMC(4samples),MCPY(3 samples),FB(3)
SNDA-PD, SNDA-HC, SNDA-ILB
"""
import os
import sys
import string
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
level_mark = 'S'
braincode_root_dir = '/PHShome/tw786/neurogen/Tao/BRAINCODE_output/'
ppm_thr = 1
reads_thr = 10
################Read good qc reads################
#Reference: write_goodReadsVStotalReads.py
braincode_qc_reads_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/BRAINCODE_allSamples_reads_statistics.txt'
braincode_qc_fp = open(braincode_qc_reads_path)
braincode_qc_good_dic = {}
while True:
	line = braincode_qc_fp.readline()
	if not line:
		break
	linesplit = line.strip().split('\t')
	sample_name = linesplit[1]
	good_qc_num = string.atoi(linesplit[2])
	braincode_qc_good_dic[sample_name] = good_qc_num
braincode_qc_fp.close()

################Extract viruses###################
braincode_brain_regions = ['SNDA', 'TCPY', 'PBMC', 'MCPY', 'FB'] 
brain_regions_virus_count_dic = {}
SNDA_HC_virus_count = {}
SNDA_PD_virus_count = {}
SNDA_ILB_virus_count = {}
for brain_region in braincode_brain_regions:
	virus_count = {}

	samples_root_dir = braincode_root_dir
	if not os.path.isdir(samples_root_dir):
		print samples_root_dir, 'does not exist!'
		sys.exit()
	for sample_name in os.listdir(samples_root_dir):
		if not os.path.isdir(samples_root_dir + sample_name):
			continue
		if  brain_region not in sample_name:
			continue
		kraken_output_path = samples_root_dir + sample_name + '/kraken_output.report'

		if not os.path.isfile(kraken_output_path):
			print kraken_output_path
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
				if not braincode_qc_good_dic.has_key(sample_name):
					continue
				ppm = float(reads_num)*1000000/braincode_qc_good_dic[sample_name]
				if viruses_start and ppm >= ppm_thr:
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
					if brain_region == 'SNDA' and 'HC' in sample_name:
						if not SNDA_HC_virus_count.has_key(name):
							SNDA_HC_virus_count[name] = {}
							SNDA_HC_virus_count[name][sample_name] = reads_num
						else:
							SNDA_HC_virus_count[name][sample_name] = reads_num
					if brain_region == 'SNDA' and 'ILB' in sample_name:
						if not SNDA_ILB_virus_count.has_key(name):
							SNDA_ILB_virus_count[name] = {}
							SNDA_ILB_virus_count[name][sample_name] = reads_num
						else:
							SNDA_ILB_virus_count[name][sample_name] = reads_num
					if brain_region == 'SNDA' and 'PD' in sample_name:
						if not SNDA_PD_virus_count.has_key(name):
							SNDA_PD_virus_count[name] = {}
							SNDA_PD_virus_count[name][sample_name] = reads_num
						else:
							SNDA_PD_virus_count[name][sample_name] = reads_num
					if name == 'Murine leukemia virus':
						print sample_name
	brain_regions_virus_count_dic[brain_region] = virus_count

brain_regions_virus_count_dic['SNDA-ILB'] = SNDA_ILB_virus_count
brain_regions_virus_count_dic['SNDA-HC'] = SNDA_HC_virus_count
brain_regions_virus_count_dic['SNDA-PD'] = SNDA_PD_virus_count


################ Braincode virus/rows sorting##################
# sort rows/viruses by the number of samples they are detected(above threshold)
virus_sampleNum_dic = {}
for brain_region in braincode_brain_regions:
	virus_count = brain_regions_virus_count_dic[brain_region]
	for virus in virus_count:
		sampleNum = len(virus_count[virus])
		if virus_sampleNum_dic.has_key(virus):
			virus_sampleNum_dic[virus] += sampleNum
		else:
			virus_sampleNum_dic[virus] = sampleNum

virus_index = sorted(virus_sampleNum_dic.items(), key = lambda x:x[1])
print virus_index
virus_index = [x[0] for x in virus_index]

braincode_brain_regions.extend(['SNDA-ILB','SNDA-HC','SNDA-PD'])
################Plot################
fontsize = 8
subplot_ind = 1
vind = np.arange(0, len(virus_index), 1)
fig = plt.figure(tight_layout=False, figsize=(28,14))
width = 0.35
#SNDA(107 samples),TCPY(10 samples),PBMC(4samples),MCPY(3 samples),FB(3)
title_dic = {
	'SNDA':'SNDA(107)',
	'TCPY':'TCPY(10)',
	'PBMC':'PBMC(4)',
	'MCPY':'MCPY(3)',
	'FB':'FB(3)',
	'SNDA-PD':'SNDA-PD(18)',
	'SNDA-HC':'SNDA-HC(61)',
	'SNDA-ILB':'SNDA-ILB(28)'
}
for brain_region in braincode_brain_regions:
	virus_count_dic = brain_regions_virus_count_dic[brain_region]
	virus_sampleNum = []
	for virus in virus_index:
		if virus_count_dic.has_key(virus):
			virus_sampleNum.append(len(virus_count_dic[virus]))
		else:
			virus_sampleNum.append(0)
	ax = fig.add_subplot(1,len(braincode_brain_regions),subplot_ind)
	ax.barh(vind, virus_sampleNum, align='center', color='b')
	ax.set_title(title_dic[brain_region], fontsize=fontsize)
	if subplot_ind != 1:
		ax.yaxis.set_visible(False)
	#bar_max = np.max(virus_sampleNum)
	bar_max = 10
	ax.set_xticks(np.arange(0,bar_max+1,1))
	ax_xticklabels = [0]
	ax_xticklabels.extend(['']*(bar_max-1))
	ax_xticklabels.append(bar_max)
	ax.set_xticklabels(ax_xticklabels)
	ax.set_ylim(-1, len(virus_index)+1)
	ax.set_xlim(0,10)
	if subplot_ind == 1:
		plt.yticks(vind+width/2., virus_index, fontsize=fontsize)

	subplot_ind += 1
fig.suptitle("BRAINCODE(>1ppm)", fontsize = 'medium')

save_folder = "/PHShome/tw786/localView/overview/May16"
if not os.path.isdir(save_folder):
	os.mkdir(save_folder)
savepath = save_folder + "/BC_8barplots(1ppm).pdf"
plt.savefig(savepath)

