"""
Make 13 GTEx barplots sidebyside

reference gtex_BarPlot_1.py gtex_BarPlot_2.py
"""


import os
import sys
import string
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

brain_regions = ['Amygdala', 'AnteriorCingulateCortex', 'Caudate', 'CerebellarHemisphere', 'Cerebellum', 'Cortex', 'FrontalCortex', 'Hippocampus', 'Hypothalamus', 'NucleusAccumbens', 'Putamen', 'SpinalCord', 'SubstantialNigra']
brain_regions.sort()
level_mark = 'S'
root_dir = '/PHShome/tw786/neurogen/Tao/GTEx_output/'
#virus_read_thr = 10

qc_goodseq_path = '/PHShome/tw786/neurogen/Tao/GTEx_output/QC_goodseqs.txt'
qc_fp = open(qc_goodseq_path)
goodseq_num_dic = {}
while True:
	line = qc_fp.readline()
	if not line:
		break
	linesplit = line.strip().split('\t')
	goodseq_num_dic[linesplit[1]] = string.atoi(linesplit[2])

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
				ppm = round(float(reads_num*1000000)/goodseq_num_dic[sample_name], 3)
				if viruses_start and ppm > 2:
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


################Plot################
fontsize = 8
width = 0.35

fig = plt.figure(tight_layout=False, figsize=(20,10))
vind = np.arange(0, len(virus_index), 1)

plt.ylim(-1, len(virus_index)+1)
plt.yticks(vind, virus_index, fontsize=fontsize)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')
subplot_ind = 1
for brain_region in brain_regions:
	gtex_virus_count_dic = brain_regions_virus_count_dic[brain_region]
	gtex_virus_sampleNum = []
	for virus in virus_index:
		if gtex_virus_count_dic.has_key(virus):
			gtex_virus_sampleNum.append(len(gtex_virus_count_dic[virus]))
		else:
			gtex_virus_sampleNum.append(0)
	ax = fig.add_subplot(1,13,subplot_ind)
	ax.barh(vind, gtex_virus_sampleNum, align='center', color='g')
	ax.set_title(brain_region, fontsize=fontsize)
	ax.yaxis.set_visible(False)
	bar_max = np.max(gtex_virus_sampleNum)
	print bar_max
	ax.set_xticks(np.arange(0,bar_max+1,1))
	ax_xticklabels = [0]
	ax_xticklabels.extend(['']*(bar_max-1))
	ax_xticklabels.append(bar_max)
	ax.set_xticklabels(ax_xticklabels, fontsize='small')
	ax.set_ylim(-1, len(virus_index)+1)
	subplot_ind += 1
fig.suptitle("GTEx(p.p.m>2)", fontsize = 'medium')


save_folder = "/PHShome/tw786/localView/overview/test"
if not os.path.isdir(save_folder):
	os.mkdir(save_folder)
savepath = save_folder + "/13barplots.pdf"
plt.savefig(savepath)

