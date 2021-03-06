"""
make a statistics what kinds of viruses and bacterias exist in all samples.

Is there any difference among samples?

we can do the analysis on different levels(species, genus, family)

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

	sample_list = []
	viruses_list = []
	bacteria_list = []
	eukaryota_list = []
	archaea_list = []

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
				if name == 'Enterobacteria phage phiX174 sensu lato':
					continue
				reads_num = string.atoi(linesplit[1].strip())
				ppm = round(reads_num/x_reads_per_million, 3)
				# if reads_num < reads_thr: # reads number threshold is appied here!
				# 	continue
				if viruses_start and ppm >= virus_ppm_thr:
					viruses_dic[name] = ppm
					viruses_list.append(name)
					
				if bacteria_start and ppm >= bacteria_ppm_thr:
					bacteria_dic[name] = ppm
					bacteria_list.append(name)
				
				if eukaryota_start and ppm >= eukaryota_ppm_thr:
					eukaryota_dic[name] = ppm
					eukaryota_list.append(name)

				if archaea_start and ppm >= archaea_ppm_thr:
					archaea_dic[name] =  ppm
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

	
	table_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/kraken_output.table_'+level_mark +'_'+ str(ppm_thr)+'.xls'
	workbook = xlwt.Workbook()

	type_dic = {'Virus':viruses_list, 'Bacteria':bacteria_list, 'Archaea':archaea_list, 'Eukaryota':eukaryota_list}
	
	
	for D_type in type_dic:
		worksheet = workbook.add_sheet(D_type)
		
		sample_list_switch = {}
		for sample in sample_list:
			switch = 0
			for species in samples_dic[sample][D_type]:
				if samples_dic[sample][D_type][species] > 0:
					switch = 1
					break
			sample_list_switch[sample] = switch

		worksheet_write_header(worksheet, sample_list, sample_list_switch)
		
		r = 1
		for species in type_dic[D_type]:
			tmp_list = []
			tmp_list.append(species)
			samples_count = 0
			for sample in sample_list:
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



	
	
	













