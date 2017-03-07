"""
Boxplot
Relationship between 10 reads and 1 p.p.m
for each sample: 10 reads = ? ppm , 1 ppm = ? reads
ppm = number of viral reads * 1M / number of good unmapped reads
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
level_mark = 'S'
braincode_root_dir = '/PHShome/tw786/neurogen/Tao/BRAINCODE_output/'
gtex_root_dir = '/PHShome/tw786/neurogen/Tao/GTEx_output/'
braincode_reads_summary_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/BRAINCODE_allSamples_reads_statistics.txt'
gtex_reads_summary_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/GTEx_allSamples_reads_statistics.txt'

braincode_fp = open(braincode_reads_summary_path)
braincode_good_unmapped_reads_dic = {}
while True:
	line = braincode_fp.readline()
	if not line:
		break
	linesplit = line.strip().split('\t')
	sample_name = linesplit[1]
	good_unmapped_reads_number = string.atoi(linesplit[2])
	braincode_good_unmapped_reads_dic[sample_name] = good_unmapped_reads_number

braincode_good_unmapped_reads_list = []
for sample_name in os.listdir(braincode_root_dir):
	if braincode_good_unmapped_reads_dic.has_key(sample_name):
		braincode_good_unmapped_reads_list.append(braincode_good_unmapped_reads_dic[sample_name])
	else:
		print sample_name
bc_how_many_ppm_of_10reads = []
for good_unmapped_number in braincode_good_unmapped_reads_list:
	ppm_value = 10000000/float(good_unmapped_number)
	bc_how_many_ppm_of_10reads.append(ppm_value)
bc_how_many_reads_of_1ppm = []
for good_unmapped_number in braincode_good_unmapped_reads_list:
	reads_number = good_unmapped_number/1000000
	bc_how_many_reads_of_1ppm.append(reads_number)

gtex_fp = open(gtex_reads_summary_path)
header = gtex_fp.readline()
gtex_how_many_ppm_of_10reads = []
gtex_how_many_reads_of_1ppm = []
while True:
	line = gtex_fp.readline()
	if not line:
		break
	linesplit = line.split('\t')
	good_unmapped_number = string.atoi(linesplit[2])
	ppm_value = 10000000/float(good_unmapped_number)
	gtex_how_many_ppm_of_10reads.append(ppm_value)
	reads_number = good_unmapped_number/1000000
	gtex_how_many_reads_of_1ppm.append(reads_number)




fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
plt.boxplot([bc_how_many_ppm_of_10reads, gtex_how_many_ppm_of_10reads])
ax1.set_xticks([1,2])
ax1.set_xticklabels(['BRAINCODE','GTEx'])
ax1.set_title("10 reads = ? ppm")

ax2 = fig.add_subplot(1,2,2)
plt.boxplot([bc_how_many_reads_of_1ppm, gtex_how_many_reads_of_1ppm])
ax2.set_xticks([1,2])
ax2.set_xticklabels(['BRAINCODE','GTEx'])
ax2.set_title("? reads = 1 ppm")
fig_save_path = '/PHShome/tw786/localView/overview/ppmVS10reads.pdf'
plt.savefig(fig_save_path)

