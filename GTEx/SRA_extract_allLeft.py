#!/bin/python
import os
sra_fp = open('/home/tw83/twang/GTEx/dbgap/allLeft/SraRunTable-6.txt', 'r')
header = sra_fp.readline()
body_site_dic = {}
while True:
	line = sra_fp.readline()
	if not line:
		break
	linesplit = line.split('\t')

	srr_name = linesplit[14]
	body_site = linesplit[18]
	body_site = body_site[8:]
	body_site_split = body_site.split(' ')
	body_site = '_'.join(body_site_split)
	body_site_split = body_site.split('_(')
	body_site = body_site_split[0]
	if not body_site_dic.has_key(body_site):
		body_site_dic[body_site] = [srr_name]
	else:
		body_site_dic[body_site].append(srr_name)
progress_dir = '/home/tw83/twang/GTEx/dbgap/allLeft/dbGaP-6101/extract_progress'
for body_site in body_site_dic:
	if body_site == 'Cerebellum':
		continue
	if body_site == 'Frontal_Cortex':
		continue
	if body_site != 'Hippocampus':
		continue

	print body_site, len(body_site_dic[body_site])
	for srr_name in body_site_dic[body_site]:
		if not os.path.isfile(progress_dir+ '/' + srr_name):
			command = "bsub -q short -W 12:00 -M 8000 -R 'rusage[mem=8000]' bash SRA_extract_allLeft_core.sh %s %s" % (body_site, srr_name)
			#print command
			os.system(command)
			command = "touch %s/%s" %(progress_dir, srr_name) 
			os.system(command) 
		#break
	#break



