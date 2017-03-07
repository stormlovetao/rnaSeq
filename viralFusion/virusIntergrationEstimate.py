"""
In kraken_output, try to find reads containing human(9606) k-mer and Virus k-mer, e.g. Human Coronavirus virus (11137)
"""

import os
import sys
import string
import xlwt
import argparse

homos_taxid = '9606'
def worksheet_write_header(worksheet, sample_list):
	worksheet.write(0,0,'Samples')
	i = 1
	for sample_name in sample_list:
		worksheet.write(0, i, sample_name)
		i+=1

def worksheet_write(worksheet, r, tmp_list):
	c = 0
	for x in tmp_list:
		str_x = str(x)
		worksheet.write(r, c, str_x)
		c += 1
def findEvidence(sample_name, viruses_evidence_dic):
	kraken_output = "/PHShome/tw786/neurogen/Tao/kraken_output/" + sample_name + "/kraken_output"
	fp = open(kraken_output, 'r')
	while True:
		line = fp.readline()
		if not line:
			break
		linesplit = line.strip().split('\t')
		pattern = linesplit[4]
		if homos_taxid not in pattern:
			continue
		pattern_split = pattern.split()
		taxid_list = [x.split(':')[0] for x in pattern_split]
		taxid_list = list(set(taxid_list))
		for taxid in taxid_list:
			if taxid == '0' or taxid == homos_taxid:
				continue
			else:
				if not viruses_evidence_dic.has_key(taxid):
					viruses_evidence_dic[taxid] = 1
				else:
					viruses_evidence_dic[taxid] += 1 
	fp.close()


if __name__ == '__main__':



	cmd_parser = argparse.ArgumentParser()
	cmd_parser.add_argument('-l', action = 'store', dest = "level_mark",
							help = "essential, S/G/F")
	cmd_parser.add_argument('-p', action = 'store', dest = "ppm_thr", type = float,
							help = "essential, global reads threshold")
	cmd_parser.add_argument('-s', action = 'store', dest = "samples_thr", type = int,
							help = "essential, global samples threshold")
	cmd_parser.add_argument('-vp', action = 'store', dest = "virus_ppm_thr",type = float,
							help = "reads number threshold set for viruses")
	cmd_parser.add_argument('-bp', action = 'store', dest = "bacteria_ppm_thr",type = float,
							help = "reads number threshold set for bacteria")
	cmd_parser.add_argument('-vs', action = 'store', dest = "virus_samples_thr",type = int,
							help = "samples number threshold set for viruses")
	cmd_parser.add_argument('-bs', action = 'store', dest = "bacteria_samples_thr",type = int,
							help = "samples number threshold set for bacteria")

	result = cmd_parser.parse_args()

	if result.level_mark == None:
		print "No level specified, use argument '-h' for help!"
		sys.exit(1)
	else:
		level_mark = result.level_mark
	if result.ppm_thr == None :
		print "No reads threshold specified, use argument '-h' for help!"
		sys.exit(1)
	else:
		ppm_thr = result.ppm_thr
	if result.samples_thr == None :
		print "No samples threshold specified, use argument '-h' for help!"
		sys.exit(1)
	else:
		samples_thr = result.samples_thr
	#set XXX_samples_thr
	virus_samples_thr = samples_thr
	bacteria_samples_thr = samples_thr
	archaea_samples_thr = samples_thr
	eukaryota_samples_thr = samples_thr

	virus_ppm_thr = ppm_thr
	bacteria_ppm_thr = ppm_thr
	archaea_ppm_thr = ppm_thr
	eukaryota_ppm_thr = ppm_thr

	if result.virus_samples_thr != None:
		virus_samples_thr = string.atof(result.virus_samples_thr)
	if result.bacteria_samples_thr != None:
		bacteria_samples_thr = string.atof(result.bacteria_samples_thr)
	if result.virus_ppm_thr != None:
		virus_ppm_thr = string.atof(result.virus_ppm_thr)
	if result.bacteria_ppm_thr != None:
		bacteria_ppm_thr = string.atof(result.bacteria_ppm_thr)


	type_samples_thr_dic = {'Virus':virus_samples_thr, 'Bacteria':bacteria_samples_thr, 'Archaea':archaea_samples_thr, 'Eukaryota': eukaryota_samples_thr}
	
	

	if level_mark != 'S' and level_mark != 'G' and level_mark != 'F':
		print "level_mark should be S or G or F"
		print "usage: python " + sys.argv[0] + "-h for help"
		sys.exit(1)

	
	samples_root_dir = '/PHShome/tw786/neurogen/Tao/kraken_output/'
	samples_dic = {}
	samples_evidence_dic = {}

	sample_list = []
	viruses_list = []
	bacteria_list = []
	eukaryota_list = []
	archaea_list = []

	viruses_possible_intergation_dic = {}

	sample_reads_filepath = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/sampleIDtotal_sequenced_reads.txt'
	sample_reads_dic = {}
	sample_reads_fp = open(sample_reads_filepath, 'r')
	for i in range(3):
		line = sample_reads_fp.readline()
	while True:
		line = sample_reads_fp.readline()
		if not line:
			break
		linesplit = line.strip().split()
		sample_name = linesplit[0].strip()
		sample_reads_count = string.atoi(linesplit[1].replace(',',''))
		sample_reads_dic[sample_name] = sample_reads_count

	it = 1
	for x in os.listdir(samples_root_dir):
		if not os.path.isdir(samples_root_dir + x):
			continue
		if not sample_reads_dic.has_key(x):
			print x,"don't know its total reads number"
			continue
		x_reads_num = sample_reads_dic[x]
		x_reads_per_million = float(x_reads_num)/1000000
		
		kraken_output_path = samples_root_dir + x + '/kraken_output.report'
		if not os.path.isfile(kraken_output_path):
			continue

		print x
		viruses_evidence_dic = {} # taxid -> evidence#
		findEvidence(x, viruses_evidence_dic)
		# print "end of findEvidence"
		# if it < 4:
		# 	print viruses_evidence_dic
		# 	it+=1
		# else:
		# 	break
		
		sample_list.append(x)
		kraken_output_fp = open(kraken_output_path,'r')

		
		viruses_dic = {}

		eukaryota_start = False
		viruses_start = False
		bacteria_start = False
		archaea_start = False

		while True:
			line = kraken_output_fp.readline()
			if not line:
				break
			linesplit = line.split('\t')
			taxid = linesplit[4]
			if taxid == '2': # taxonimic id of bacteria is 2 
				eukaryota_start = False
				viruses_start = False
				bacteria_start = True
				archaea_start = False
				continue

			if taxid == '10239': #Virus(taxonimic id = 10239)
				eukaryota_start = False
				viruses_start = True
				bacteria_start = False
				archaea_start = False
				continue
			if taxid == '2157': #archaea(taxonimic id = 2157)
				eukaryota_start = False
				viruses_start = False
				bacteria_start = False
				archaea_start = True
				continue
			if taxid == '2759': #eukaryota(taxonimic id = 2759)
				eukaryota_start = True
				viruses_start = False
				bacteria_start = False
				archaea_start = False
				continue

			if linesplit[3] == level_mark and viruses_start:
				name = linesplit[5].lstrip().strip()
				# if name == 'Enterobacteria phage phiX174 sensu lato':
				# 	continue

				# taxid2name_dic[taxid] = name
				# name2taxid_dic[name] = taxid
				reads_num = string.atoi(linesplit[1].strip())
				ppm = round(reads_num/x_reads_per_million, 3)

		
				if not viruses_evidence_dic.has_key(taxid):
					viruses_dic[name] = (ppm, 0)
					viruses_list.append(name)
				else:
					evidence_num = viruses_evidence_dic[taxid]
					viruses_dic[name] = (ppm, evidence_num)
					viruses_list.append(name)
					# if evidence_num > 0:
					# 	print x, name, linesplit[4], evidence_num
					
		samples_dic[x] = viruses_dic
	
	

	
	#writing to tables
	viruses_list = list(set(viruses_list))
	table_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/kraken_output.VirusPpmMixedReads#.xls'
	workbook = xlwt.Workbook()
	worksheet = workbook.add_sheet('Viruses')
	worksheet_write_header(worksheet, sample_list)
	r = 1
	for virus in viruses_list:
		tmp_list = [virus]
		for sample_name in sample_list:
			if not samples_dic[sample_name].has_key(virus):
				ppm_evi = (0,0)
			else:
				ppm_evi = samples_dic[sample_name][virus]
			tmp_list.append(ppm_evi)
		worksheet_write(worksheet, r, tmp_list)
		r += 1

	workbook.save(table_path)
















		