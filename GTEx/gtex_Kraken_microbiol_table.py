"""
Reference: CallVirusBac.py
Make a table: rows indicate microbime, colomns show samples, each cell represents reads number above a threshold   
"""

import os
import sys
import string
import xlwt
import argparse

def worksheet_write_header(worksheet, sample_list, sample_list_switch):
	worksheet.write(0,0,'Samples')
	i = 1
	for sample_name in sample_list:
		if sample_list_switch[sample_name] == 0:
			continue
		worksheet.write(0, i, sample_name)
		i+=1

def worksheet_write(worksheet, r, tmp_list):
	c = 0
	for x in tmp_list:
		worksheet.write(r, c, x)
		c += 1

if __name__ == '__main__':



	cmd_parser = argparse.ArgumentParser()
	cmd_parser.add_argument('-l', action = 'store', dest = "level_mark",
							help = "essential, S/G/F")
	
	cmd_parser.add_argument('-s', action = 'store', dest = "samples_thr", type = int,
							help = "global samples threshold, default=1")
	cmd_parser.add_argument('-rt', action = 'store', dest = "reads_thr",type = int,
							help = "reads number threshold set for viruses, default=0")


	result = cmd_parser.parse_args()

	if result.level_mark == None:
		print "No level specified, use argument '-h' for help!"
		sys.exit(1)
	else:
		level_mark = result.level_mark
	
	if result.samples_thr == None :
		samples_thr = 1
	else:
		samples_thr = string.atoi(result.samples_thr)
	if result.reads_thr == None:
		reads_thr = 0
	else:
		reads_thr = string.atoi(result.reads_thr)

	#set XXX_samples_thr
	virus_samples_thr = samples_thr
	bacteria_samples_thr = samples_thr
	archaea_samples_thr = samples_thr
	eukaryota_samples_thr = samples_thr

	

	type_samples_thr_dic = {'Virus':virus_samples_thr, 'Bacteria':bacteria_samples_thr, 'Archaea':archaea_samples_thr, 'Eukaryota': eukaryota_samples_thr}
	
	

	if level_mark != 'S' and level_mark != 'G' and level_mark != 'F':
		print "level_mark should be S or G or F"
		print "usage: python " + sys.argv[0] + "-h for help"
		sys.exit(1)

	
	samples_root_dir = '/PHShome/tw786/neurogen/Tao/GTEx_output/'
	#samples_root_dir = '/PHShome/tw786/neurogen/Tao/kraken_output/'
	samples_dic = {}

	sample_list = []
	viruses_list = []
	bacteria_list = []
	eukaryota_list = []
	archaea_list = []

	
	for x in os.listdir(samples_root_dir):
		if not os.path.isdir(samples_root_dir + x):
			continue
		kraken_output_path = samples_root_dir + x + '/kraken_output.report'

		if not os.path.isfile(kraken_output_path):
			continue
		sample_list.append(x)
		kraken_output_fp = open(kraken_output_path,'r')

		VB_dic = {}
		eukaryota_dic = {}
		viruses_dic = {}
		bacteria_dic = {}
		archaea_dic = {}
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
				
				if viruses_start and reads_num >= reads_thr:
					viruses_dic[name] = reads_num
					viruses_list.append(name)
					
				if bacteria_start and reads_num >= reads_thr:
					bacteria_dic[name] = reads_num
					bacteria_list.append(name)
				
				if eukaryota_start and reads_num >= reads_thr:
					eukaryota_dic[name] = reads_num
					eukaryota_list.append(name)

				if archaea_start and reads_num >= reads_thr:
					archaea_dic[name] =  reads_num
					archaea_list.append(name)

		
		VB_dic['Virus'] = viruses_dic
		VB_dic['Bacteria'] = bacteria_dic
		VB_dic['Eukaryota'] = eukaryota_dic
		VB_dic['Archaea'] = archaea_dic

		samples_dic[x] = VB_dic
		

	

	#write into a table
	

	viruses_list = list(set(viruses_list))
	bacteria_list = list(set(bacteria_list))
	archaea_list = list(set(archaea_list))
	eukaryota_list = list(set(eukaryota_list))
	sample_list = sorted(list(set(sample_list)))

	table_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/gtexTotalReadsNumber.xls'
	#table_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/BRAINCODETotalReadsNumber.xls'
	workbook = xlwt.Workbook()

	type_dic = {'Virus':viruses_list, 'Bacteria':bacteria_list, 'Archaea':archaea_list, 'Eukaryota':eukaryota_list}
	
	
	for D_type in type_dic:
		worksheet = workbook.add_sheet(D_type)
		
		# sort colomn/sample by average/sum of each colomn
		sample_list_switch = {}
		colomn_sum_dic = {}
		for sample in sample_list:
			switch = 0
			colomn_sum = 0
			for species in samples_dic[sample][D_type]:
				if samples_dic[sample][D_type][species] > 0:
					switch += 1
					colomn_sum += samples_dic[sample][D_type][species]
			sample_list_switch[sample] = switch
			colomn_sum_dic[sample] = colomn_sum
		new_sample_list = []
		new_sample_list_speciesNum = []
		new_sample_list_speciesSum = []
		for sample in sample_list_switch:
			if sample_list_switch[sample] > 0:
				new_sample_list.append(sample)
				new_sample_list_speciesNum.append(sample_list_switch[sample])
				new_sample_list_speciesSum.append(colomn_sum_dic[sample])

		new_sample_list_speciesSum,new_sample_list=(list(x) for x in zip(*sorted(zip(new_sample_list_speciesSum, new_sample_list),reverse=True)))
		if D_type == 'Virus':
			print new_sample_list_speciesSum
			print new_sample_list
		worksheet_write_header(worksheet, new_sample_list, sample_list_switch)
		
		# sort row/species by average/sum of each row
		v_list = type_dic[D_type]
		v_list_sampleSum = []
		# how many sample contain this species
		for species in v_list:
			sample_count = 0
			for sample in new_sample_list:
				if 	samples_dic[sample][D_type].has_key(species):
					sample_count += samples_dic[sample][D_type][species]
			v_list_sampleSum.append(sample_count)
		v_list_sampleSum, v_list = (list(x) for x in zip(*sorted(zip(v_list_sampleSum, v_list),reverse=True)))



		r = 1
		for species in v_list:
			tmp_list = []
			tmp_list.append(species)
			samples_count = 0
			for sample in new_sample_list:
				if sample_list_switch[sample] == 0:
					continue
				if samples_dic[sample][D_type].has_key(species):
					samples_count += 1
					tmp_list.append(samples_dic[sample][D_type][species])
				else:
					tmp_list.append(0)
			if samples_count >= type_samples_thr_dic[D_type]:
				worksheet_write(worksheet, r, tmp_list)
				r += 1
	workbook.save(table_path)