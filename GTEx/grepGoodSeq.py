
import os
import sys
import string
root_dir = '/home/tw83/twang/GTEx_output/'

brain_regions = ['Amygdala', 'AnteriorCingulateCortex', 'Caudate', 'CerebellarHemisphere', 'Cerebellum', 'Cortex', 'FrontalCortex', 'Hippocampus', 'Hypothalamus', 'NucleusAccumbens', 'Putamen', 'SpinalCord', 'SubstantialNigra']
file_out_path = root_dir + 'QC_goodseqs.txt'
out_fp = open(file_out_path,'w')
for brain_region in brain_regions:
	brain_region_dir = root_dir + brain_region + '/'
	if not os.path.isdir(brain_region_dir):
		print brain_region_dir, 'does not exist'
		sys.exit()
	for sample_name in os.listdir(brain_region_dir):
		if not sample_name.startswith('SRR'):
			continue
		qc_logfile_path = brain_region_dir+ sample_name+'/unmapped.fastq.log'
		if not os.path.isfile(qc_logfile_path):
			print qc_logfile_path, 'does not exist'
			continue
		fp = open(qc_logfile_path,'r')
		while True:
			line = fp.readline()
			if not line:
				break
			if 'Good sequences' in line:
				linesplit = line.split(' ')
				good_num = linesplit[5]
				good_num = good_num.replace(',','')
				good_num = string.atoi(good_num)
				if good_num > 0:
					outline = brain_region + '\t' + sample_name + '\t' + str(good_num) + '\n'
					out_fp.write(outline)
				break
