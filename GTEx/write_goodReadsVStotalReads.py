"""
scp to Orchestra and Erisone
For GTEx data:
Plot good unmapped reads number v.s. total unmapped reads number
Plot good unmapped reads number v.s. total library reads number
Plot total unmapped reads number v.s. total library reads number

For BRAINCODE data:
Plot good unmapped reads number v.s. total unmapped reads number
Plot good unmapped reads number v.s. total library reads number
Plot total unmapped reads number v.s. total library reads number
"""


import os
import sys
import string

if os.path.isdir('/PHShome/tw786'):
	root_dir = '/PHShome/tw786/neurogen/Tao/BRAINCODE_output'
	readsinfo_from_googledoc_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/Reads_info_from_googledoc.txt'
	file_out_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/BRAINCODE_allSamples_reads_statistics.txt'
	out_fp = open(file_out_path, 'w')
	header_line = "BRAINCODE\tSampleName\tReadsCount(unmapped pass QC)\tReadsCount(total unmapped)\tReadsCount(raw data)\tReadsCount(raw data pass fastqc)\n"
	# class BrainSample:
	# 	def __init__(self, sample_name, total_library_reads_number, qcpass_library_reads_number, unmapped_reads_number):
	# 		self.sample_name = sample_name
	# 		self.total_library_reads_number = total_library_reads_number
	# 		self.qcpass_library_reads_number = qcpass_library_reads_number
	# 		self.unmapped_reads_number = unmapped_reads_number
	google_dic = {}
	google_fp = open(readsinfo_from_googledoc_path)
	header_line = google_fp.readline()
	header_line = google_fp.readline()
	header_line = google_fp.readline()
	while True:
		line = google_fp.readline()
		if not line:
			break
		linesplit = line.strip().split('\t')
		sample_name = linesplit[0]
		raw_reads_number = linesplit[2]
		raw_reads_number = raw_reads_number.replace(',','')
		raw_reads_number = string.atoi(raw_reads_number)
		qcpass_reads_number =  linesplit[3]
		qcpass_reads_number = qcpass_reads_number.replace(',','')
		qcpass_reads_number = string.atoi(qcpass_reads_number)
		unmapped_reads_number = linesplit[4]
		unmapped_reads_number = unmapped_reads_number.replace(',','')
		unmapped_reads_number = string.atoi(unmapped_reads_number)
		google_dic[sample_name] = [raw_reads_number, qcpass_reads_number, unmapped_reads_number]
	for sample_name in os.listdir(root_dir):
		sample_dir = root_dir + '/' + sample_name
		if not os.path.isdir(sample_dir):
			continue
		# if not google_dic.has_key(sample_name):
		# 	continue
		qc_logfile_path = sample_dir + '/unmapped.fastq.log'
		if not os.path.isfile(qc_logfile_path):
			print qc_logfile_path, "doesnot exist"
			continue
		else:
			fp = open(qc_logfile_path,'r')
			while True:
				line = fp.readline()
				if not line:
					break
				good_num = 0
				if 'Good sequences' in line:
					linesplit = line.split(' ')
					good_num = linesplit[5]
					good_num = good_num.replace(',','')
					good_num = string.atoi(good_num)
				if 'Input sequences' in line:
					linesplit = line.split(' ')
					input_num = linesplit[5]
					input_num = input_num.replace(',','')
					input_num = string.atoi(input_num)

				if good_num > 0:
					if not google_dic.has_key(sample_name):
						google_dic_raw_reads_number = -1
						google_dic_qcpass_reads_num = -1
					else:
						google_dic_raw_reads_number = google_dic[sample_name][0]
						google_dic_qcpass_reads_num = google_dic[sample_name][1]
					outline = "BRAINCODE"+ '\t' + sample_name + '\t' + str(good_num) + '\t' + str(input_num) + '\t'+ str(google_dic_raw_reads_number) +'\t'+ str(google_dic_qcpass_reads_num) +'\n'
					out_fp.write(outline)
					break



if os.path.isdir('/home/tw83/'):
	root_dir = '/home/tw83/twang/GTEx_output/'
	sra_table_path = '/home/tw83/twang/GTEx_output/SraRunTable-allGTExBrainLiverSamples.txt'
	srr_total_library_Mbase_dic = {}
	sra_fp = open(sra_table_path)
	header_line = sra_fp.readline()
	while True:
		line = sra_fp.readline()
		if not line:
			break
		linesplit = line.split('\t')
		srr_id = linesplit[14]
		Mbase = linesplit[11]
		srr_total_library_Mbase_dic[srr_id] = Mbase
	sra_fp.close()
	brain_regions = ['Amygdala', 'AnteriorCingulateCortex', 'Caudate', 'CerebellarHemisphere', 'Cerebellum', 'Cortex', 'FrontalCortex', 'Hippocampus', 'Hypothalamus', 'NucleusAccumbens', 'Putamen', 'SpinalCord', 'SubstantialNigra','Liver']
	file_out_path = root_dir + 'GTEx_allSamples_reads_statistics.txt'
	out_fp = open(file_out_path,'w')
	header_line = 'BrainRegions\tRunName\tReadsCount pass QC\tReadsCount(total unmapped)\tMBaseCount(M,total library)\n'
	out_fp.write(header_line)
	input_seqnum_list = []
	good_seqnum_list = []
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
				good_num = 0
				if 'Good sequences' in line:
					linesplit = line.split(' ')
					good_num = linesplit[5]
					good_num = good_num.replace(',','')
					good_num = string.atoi(good_num)
				if 'Input sequences' in line:
					linesplit = line.split(' ')
					input_num = linesplit[5]
					input_num = input_num.replace(',','')
					input_num = string.atoi(input_num)

				if good_num > 0:
					if not srr_total_library_Mbase_dic.has_key(sample_name):
						out_mbase = 0
						print sample_name
					else:
						out_mbase = srr_total_library_Mbase_dic[sample_name]
					outline = brain_region + '\t' + sample_name + '\t' + str(good_num) + '\t' + str(input_num) + '\t'+ str(out_mbase) +'\n'
					out_fp.write(outline)
					
					break








