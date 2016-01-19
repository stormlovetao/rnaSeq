"""
make a statistics what kinds of viruses and bacterias exist in all samples.

Is there any difference among samples?

we can do the analysis on different levels(species, genus, family)

"""
import os
import sys
import string
if __name__ == '__main__':

	if len(sys.argv) != 2:
		print "usage: python " + sys.argv[0] + "S/G/F"
		sys.exit(1)


	level_mark = sys.argv[1]

	if level_mark != 'S' and level_mark != 'G' and level_mark != 'F':
		print "usage: python " + sys.argv[0] + "S/G/F"
		sys.exit(1)

	
	samples_root_dir = '/PHShome/tw786/neurogen/Tao/kraken_output/'
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
				if viruses_start:
					viruses_dic[name] = reads_num
					viruses_list.append(name)
					
					
				if bacteria_start:
					bacteria_dic[name] = reads_num
					bacteria_list.append(name)
				
				if eukaryota_start:
					eukaryota_dic[name] = reads_num
					eukaryota_list.append(name)

				if archaea_start:
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
	sample_list = list(set(sample_list))
	eukaryota_list = list(set(eukaryota_list))
	table_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/kraken_output.virus.table_'+level_mark
	table = open(table_path,'w')
	header_line = 'Samples\t' + '\t'.join(sample_list) + '\n'
	table.write(header_line)

	for virus in viruses_list:
		outline = virus
		for sample in sample_list:
			if samples_dic[sample]['Virus'].has_key(virus):
				outline = outline + '\t' + str(samples_dic[sample]['Virus'][virus])
			else:
				outline = outline + '\t' + '0'
		outline = outline + '\n'
		table.write(outline)
	table.close()

	table_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/kraken_output.bacteria.table_'+level_mark
	table = open(table_path,'w')
	header_line = 'Samples\t' + '\t'.join(sample_list) + '\n'
	table.write(header_line)
	for bac in bacteria_list:
		outline = bac
		for sample in sample_list:
			if samples_dic[sample]['Bacteria'].has_key(bac):
				outline = outline + '\t' + str(samples_dic[sample]['Bacteria'][bac])
			else:
				outline = outline + '\t' + '0'
		outline = outline + '\n'
		table.write(outline)

	table.close()

	table_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/kraken_output.archaea.table_'+level_mark
	table = open(table_path,'w')
	header_line = 'Samples\t' + '\t'.join(sample_list) + '\n'
	table.write(header_line)
	for arc in archaea_list:
		outline = arc
		for sample in sample_list:
			if samples_dic[sample]['Archaea'].has_key(arc):
				outline = outline + '\t' + str(samples_dic[sample]['Archaea'][arc])
			else:
				outline = outline + '\t' + '0'
		outline = outline + '\n'
		table.write(outline)

	table.close()




	table_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/kraken_output.eukaryota.table_'+level_mark
	table = open(table_path,'w')
	header_line = 'Samples\t' + '\t'.join(sample_list) + '\n'
	table.write(header_line)
	for homo in eukaryota_list:
		outline = homo
		for sample in sample_list:
			if samples_dic[sample]['Eukaryota'].has_key(homo):
				outline = outline + '\t' + str(samples_dic[sample]['Eukaryota'][homo])
			else:
				outline = outline + '\t' + '0'
		outline = outline + '\n'
		table.write(outline)

	table.close()














