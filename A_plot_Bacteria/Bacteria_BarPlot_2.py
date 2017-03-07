"""
Make 17 barplots sidebyside
BRAINCODE barplot all, HC, PD, 13 GTEx brainregions barplots, 1 Liver barplot
Output GTEx bacteria and BRAINCODE bacteria
reference BarPlot_2.py
some variables are still in virus** format. Actually, we use bacteria_start to contral what we read from kraken report.
"""
import os
import sys
import string
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

brain_regions = ['SubstantialNigra', 'Amygdala', 'AnteriorCingulateCortex', 'Caudate', 'CerebellarHemisphere', 'Cerebellum', 'Cortex', 'FrontalCortex', 'Hippocampus', 'Hypothalamus', 'NucleusAccumbens', 'Putamen', 'SpinalCord', 'Liver']
#brain_regions.sort()
level_mark = 'F'

braincode_root_dir = '/PHShome/tw786/neurogen/Tao/BRAINCODE_output/'
gtex_root_dir = '/PHShome/tw786/neurogen/Tao/GTEx_output/'
virus_read_thr = 10
ppm_thr = 1

###################################################
braincode_virus_count_dic = {}
if os.path.isdir(braincode_root_dir):
	for sample_name in os.listdir(braincode_root_dir):
		if not os.path.isdir(braincode_root_dir + sample_name):
			continue
		kraken_output_path = braincode_root_dir + sample_name + '/kraken_output.report'

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
			
				if bacteria_start and reads_num >= virus_read_thr:
					if not braincode_virus_count_dic.has_key(name):
						braincode_virus_count_dic[name] = {}
						braincode_virus_count_dic[name][sample_name] = reads_num
						
					else:
						braincode_virus_count_dic[name][sample_name] = reads_num

# #########################






brain_regions_virus_count_dic = {}
for brain_region in brain_regions:
	virus_count = {}
	samples_root_dir = gtex_root_dir + brain_region + '/'
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

				if bacteria_start and reads_num >= virus_read_thr:
					if not virus_count.has_key(name):
						virus_count[name] = {}
						virus_count[name][sample_name] = reads_num
					else:
						virus_count[name][sample_name] = reads_num
	brain_regions_virus_count_dic[brain_region] = virus_count

################ GTEx virus/rows sorting##################
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

virus_index = sorted(virus_sampleNum_dic.items(), key = lambda x:x[1])
virus_index = [x[0] for x in virus_index]

################ Output GTEx viruses###############
save_folder = "/PHShome/tw786/localView/overview/Tree/"
if not os.path.isdir(save_folder):
	os.mkdir(save_folder)
gtex_virus_outpath = save_folder + 'GTEx_bacteria.txt'
fp = open(gtex_virus_outpath, 'w')
for item in virus_index:
	fp.write("%s\n" % item)
fp.close()
################Plot################
fontsize = 8
braincode_virus_count_dic_items = sorted(braincode_virus_count_dic.items(), key = lambda x:len(x[1]))
braincode_virus_list = [x[0] for x in braincode_virus_count_dic_items]
braincode_virus_sampleNum = [len(x[1]) for x in braincode_virus_count_dic_items ]

################Output BRAINCODE viruses##############
braincode_virus_outpath = save_folder + 'BRAINCODE_bacteria.txt'
fp = open(braincode_virus_outpath, 'w')
for item in braincode_virus_list:
	fp.write("%s\n" % item)
fp.close()


###add viruses only found in GTEx data###
for virus in virus_index:
	if virus not in braincode_virus_list:
		braincode_virus_list.append(virus)
		braincode_virus_sampleNum.append(0)

fig = plt.figure(tight_layout=False, figsize=(30,14))
ax1 = fig.add_subplot(1,len(brain_regions)+3,1)
ax1.set_title('BRAINCODE', fontsize=fontsize)
#bar_max = np.max(braincode_virus_sampleNum)
bar_max = 10
ax1.set_xticks(np.arange(0,bar_max+1,1))
ax_xticklabels = [0]
ax_xticklabels.extend(['']*(bar_max-1))
ax_xticklabels.append(bar_max)
#ax1.set_xticks(ax_xticklabels)
ax1.set_xticklabels(ax_xticklabels)
vind = np.arange(0, len(braincode_virus_list), 1)
# print vind
# print braincode_virus_sampleNum
width = 0.35
ax1.barh(vind, braincode_virus_sampleNum,  align='center', color='r')
ax1.set_ylim(-1, len(braincode_virus_list)+1)
ax1.set_xlim(0,10)
plt.yticks(vind+width/2., braincode_virus_list, fontsize=fontsize-2)

#####For HC######
ax2 = fig.add_subplot(1,len(brain_regions)+3,2)
HC_virus_sampleNum = []
for virus in braincode_virus_list:
	if virus not in braincode_virus_count_dic:
		HC_virus_sampleNum.append(0)
	else:
		sample_dic = braincode_virus_count_dic[virus]
		sample_list = sample_dic.keys()
		HC_count = 0
		for sample_name in sample_list:
			if sample_name.startswith('HC'):
				HC_count += 1
		HC_virus_sampleNum.append(HC_count)

ax2.set_title('BRAINCODE(HC)', fontsize=fontsize)
ax2.set_xticks(np.arange(0,bar_max+1,1))
ax2.set_xticklabels(ax_xticklabels)
ax2.barh(vind, HC_virus_sampleNum, align='center', color='r')
ax2.set_ylim(-1, len(braincode_virus_list)+1)
ax2.set_xlim(0,10)
ax2.yaxis.set_visible(False)

#####For ILB######
ax3 = fig.add_subplot(1,len(brain_regions)+3,3)
ILB_virus_sampleNum = []
for virus in braincode_virus_list:
	if virus not in braincode_virus_count_dic:
		ILB_virus_sampleNum.append(0)
	else:
		sample_dic = braincode_virus_count_dic[virus]
		sample_list = sample_dic.keys()
		ILB_count = 0
		for sample_name in sample_list:
			if sample_name.startswith('ILB'):
				ILB_count += 1
		ILB_virus_sampleNum.append(ILB_count)

ax3.set_title('BRAINCODE(ILB)', fontsize=fontsize)
ax3.set_xticks(np.arange(0,bar_max+1,1))
ax3.set_xticklabels(ax_xticklabels)
ax3.barh(vind, ILB_virus_sampleNum, align='center', color='r')
ax3.set_ylim(-1, len(braincode_virus_list)+1)
ax3.set_xlim(0,10)
ax3.yaxis.set_visible(False)

########For GTEx#######
subplot_ind = 4
for brain_region in brain_regions:
	gtex_virus_count_dic = brain_regions_virus_count_dic[brain_region]
	gtex_virus_sampleNum = []
	for virus in braincode_virus_list:
		if gtex_virus_count_dic.has_key(virus):
			gtex_virus_sampleNum.append(len(gtex_virus_count_dic[virus]))
		else:
			gtex_virus_sampleNum.append(0)
	ax = fig.add_subplot(1,len(brain_regions)+3,subplot_ind)
	ax.barh(vind, gtex_virus_sampleNum, align='center', color='b')
	ax.set_title(brain_region, fontsize=fontsize)
	ax.yaxis.set_visible(False)
	bar_max = np.max(gtex_virus_sampleNum)
	ax.set_xticks(np.arange(0,bar_max+1,1))
	# ax_xticklabels = [0]
	# ax_xticklabels.extend(['']*(bar_max-1))
	# ax_xticklabels.append(bar_max)
	ax.set_xticklabels(ax_xticklabels)
	ax.set_ylim(-1, len(braincode_virus_list)+1)
	ax.set_xlim(0,10)
	subplot_ind += 1
title_str = "BRAINCODE v.s GTEx(>=%d reads)" %(virus_read_thr)
fig.suptitle(title_str, fontsize = 'medium')

save_folder = "/PHShome/tw786/localView/overview/Jun10"
if not os.path.isdir(save_folder):
	os.mkdir(save_folder)
savepath = save_folder + "/Bacteria_17barplots(%dreads)_atFamilyLevel.pdf" %(virus_read_thr)
plt.savefig(savepath)
