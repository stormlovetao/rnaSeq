"""
Reference: Heatmap_Table.py
GTEx data
output bacteria reads number >= 1
Output a table of bacteria found in gtex
Rows are bacteria
Columns are samples, labeled by sample_id, srr_id, and grouped by subject id.


"""
import os
import sys
import string


brain_regions = ['Amygdala', 'AnteriorCingulateCortex', 'Caudate', 'CerebellarHemisphere', 'Cerebellum', 'Cortex', 'FrontalCortex', 'Hippocampus', 'Hypothalamus', 'NucleusAccumbens', 'Putamen', 'SpinalCord', 'SubstantialNigra','Liver']
brain_regions.sort()
level_mark = 'F'
root_dir = '/PHShome/tw786/neurogen/Tao/GTEx_output/'

virus_read_thr = 1



############################# READ SUBJECTS #############################
brain_region_mapping_dic = {
	'Brain - Amygdala' : 'Amygdala',
	'Brain - Anterior cingulate cortex (BA24)' : 'AnteriorCingulateCortex',
	'Brain - Caudate (basal ganglia)' : 'Caudate',
	'Brain - Cerebellar Hemisphere' : 'CerebellarHemisphere',
	'Brain - Cerebellum' : 'Cerebellum',
	'Brain - Cortex' : 'Cortex',
	'Brain - Frontal Cortex (BA9)' : 'FrontalCortex',
	'Brain - Hippocampus' : 'Hippocampus',
	'Brain - Hypothalamus' : 'Hypothalamus',
	'Brain - Nucleus accumbens (basal ganglia)' : 'NucleusAccumbens',
	'Brain - Putamen (basal ganglia)' : 'Putamen',
	'Brain - Spinal cord (cervical c-1)' : 'SpinalCord',
	'Brain - Substantia nigra' : 'SubstantialNigra',
	'Liver' : 'Liver'
}
sra_table_path = '/PHShome/tw786/neurogen/Tao/GTEx_output/SraRunTable_GTExBrain_2016_4_25.txt'
sra_fp = open(sra_table_path)
subjects_dic = {}
#subjects_dic[subject_id][brain_region]=[srr_1, srr_2..]
class SRR:
	def __init__(self, srr_id, subject_id, brain_region):
		self.srr_id = srr_id
		self.subject_id = subject_id
		self.brain_region = brain_region

header_line = sra_fp.readline()
sample_count = 0
srrid2sampleid_dic = {}
srrid2brainregion_dic = {}
while True:
	line = sra_fp.readline()
	if not line:
		break
	linesplit = line.split('\t')
	subject_id = linesplit[26]
	srr_id = linesplit[14]
	sample_id = linesplit[25]
	srrid2sampleid_dic[srr_id] = sample_id

	brain_region = brain_region_mapping_dic[linesplit[18]]
	srrid2brainregion_dic[srr_id] = brain_region
	if not os.path.isdir(root_dir + brain_region + '/' + srr_id):
		continue
	else:
		sample_count += 1
	srr = SRR(srr_id, subject_id, brain_region)
	if not subjects_dic.has_key(subject_id):
		subjects_dic[subject_id] = {}
		subjects_dic[subject_id][brain_region] = [srr]
	else:
		if not subjects_dic[subject_id].has_key(brain_region):
			subjects_dic[subject_id][brain_region] = [srr]
		else:
			subjects_dic[subject_id][brain_region].append(srr)
sra_fp.close()
liver_sra_table_path = '/PHShome/tw786/neurogen/Tao/GTEx_output/SraRunTable_GTExLiver_2016_4_25.txt'
sra_fp = open(liver_sra_table_path)
header_line = sra_fp.readline()
while True:
	line = sra_fp.readline()
	if not line:
		break
	linesplit = line.split('\t')
	subject_id = linesplit[26]
	srr_id = linesplit[14]
	sample_id = linesplit[25]
	srrid2sampleid_dic[srr_id] = sample_id
	brain_region = brain_region_mapping_dic[linesplit[18]]
	srrid2brainregion_dic[srr_id] = brain_region
	if not os.path.isdir(root_dir + brain_region + '/' + srr_id):
		continue
	else:
		sample_count += 1
	srr = SRR(srr_id, subject_id, brain_region)
	if not subjects_dic.has_key(subject_id):
		subjects_dic[subject_id] = {}
		subjects_dic[subject_id][brain_region] = [srr]
	else:
		if not subjects_dic[subject_id].has_key(brain_region):
			subjects_dic[subject_id][brain_region] = [srr]
		else:
			subjects_dic[subject_id][brain_region].append(srr)
sra_fp.close()

print "we have %d samples" % sample_count 

############################# READ viruses reads FOR GTEX - EACH BRAIN REGION and Liver #############################
brain_regions_virus_count_dic = {}
srr_contain_virusnum = {}
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
		srr_contain_virusnum[sample_name] = 0
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
			
				# ppm = round(float(reads_num*1000000)/goodseq_num_dic[sample_name], 3)
				#if viruses_start and  ppm > ppm_thr:
				if bacteria_start and  reads_num >= virus_read_thr:
					
					if not virus_count.has_key(name):
						virus_count[name] = {}
						virus_count[name][sample_name] = reads_num
					else:
						virus_count[name][sample_name] = reads_num
					srr_contain_virusnum[sample_name] += 1
	brain_regions_virus_count_dic[brain_region] = virus_count
# brain_regions_virus_count_dic[brain_region][virus][sample_name]=ppm


############################# MAKE TABLE #############################
# Make a table
# Rows represnt viruses
# Columns represent subjects, subdivied by brain regions

# define row index: sort rows/viruses by the number of samples they are detected(above threshold)
table_out_path = '/PHShome/tw786/localView/GTExTables/gtex_table_subjects_BacteriaFamilyLevel_reads_%d_withSampleID.xls' % virus_read_thr
out_fp = open(table_out_path, 'w')

#######clean subjects with no virus##########
subject_has_virus = []
for subject_id in subjects_dic:
	tag = False
	for brain_region in subjects_dic[subject_id]:
		if tag:
			break
		for srr in subjects_dic[subject_id][brain_region]:
			# if not srr_contain_virusnum.has_key(srr_id):
			# 	continue
			srr_id = srr.srr_id
			if srr_contain_virusnum[srr_id] > 0:
				tag = True
				break
	if tag:
		subject_has_virus.append(subject_id)


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

###print header line1 SUBJECTs
subject_sample_count = {}
for subject_id in subjects_dic:
	count = 0
	for brain_region in subjects_dic[subject_id]:
		count += len(subjects_dic[subject_id][brain_region])
	subject_sample_count[subject_id] = count

header_line = 'Subjects:\t'
out_seq = []
for subject_id in subject_has_virus:
	print_num = subject_sample_count[subject_id]
	out_seq.extend([subject_id]*print_num)
header_line += '\t'.join(out_seq)
header_line += '\n'
out_fp.write(header_line)
###print header line2 submitted_sample_id
# Here we output all samples while in Heatmap_table.py we merged some samples sometime when one brainregion has several samples
header_line = 'Sample_id:\t'
out_seq = []
for subject_id in subject_has_virus:
	for brain_region in subjects_dic[subject_id]:
		for srr in subjects_dic[subject_id][brain_region]:
			srr_id = srr.srr_id
			out_seq.append(srrid2sampleid_dic[srr_id])
		
header_line += '\t'.join(out_seq)
header_line += '\n'
out_fp.write(header_line)
###print header line3 brain_region
header_line = 'Brain_regions:\t'
out_seq = []
for subject_id in subject_has_virus:
	for brain_region in subjects_dic[subject_id]:
		out_seq.extend([brain_region]*len(subjects_dic[subject_id][brain_region]))
header_line += '\t'.join(out_seq)
header_line += '\n'
out_fp.write(header_line)
####print Rows
for virus in virus_index:
	outline = virus + '\t'
	out_seq = []
	for subject_id in subject_has_virus:
		for brain_region in subjects_dic[subject_id]:
			for srr in subjects_dic[subject_id][brain_region]:
				reads_num = 0
				srr_id = srr.srr_id
				if brain_regions_virus_count_dic[brain_region].has_key(virus):
					if brain_regions_virus_count_dic[brain_region][virus].has_key(srr_id):
						reads_num = brain_regions_virus_count_dic[brain_region][virus][srr_id]
				out_seq.append(reads_num)

	outline = outline + '\t'.join([str(x) for x in out_seq]) + '\n'
	out_fp.write(outline)

