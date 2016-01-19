#!/usr/bin/python


"""
remove reads in AllCircRnaList.txt and AllGeneFusionReads.txt from unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.sam
output new sam file unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.filterCF.sam

"""



import sys
import os

if __name__ == '__main__':
	
	if len(sys.argv) != 2:
		print "usage: python " + sys.argv[0] + " [sample_name]"
		sys.exit(1)

	sample_name = sys.argv[1]

	root_path = '/data/neurogen/Tao/kraken_output/'
	sample_path = root_path + sample_name + '/'

	if not os.path.exists(sample_path):
		print sample_path, "doesn't exist!"
		sys.exit(1)

	samfile_path = sample_path + 'unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.sam'
	if not os.path.isfile(samfile_path):
		samfile_path, "doesn't exist!"
		sys.exit(1)
	samfile = open(samfile_path, 'r')

	new_samfile_path = sample_path + 'unmapped.filterLC.filterPhiX.filterHg38.filterLSU_SSU.filterCF.sam'

	samout = open(new_samfile_path, 'w')

	circRna_list_path = sample_path + 'AllCircRnaList.txt'
	cirRna_fp = open(circRna_list_path, 'r')
	cirRna_list = []
	while True:
		line = cirRna_fp.readline()
		if not line:
			break
		line = line.strip()
		if line == '':
			continue
		cirRna_list.append(line)

	cirRna_list = set(cirRna_list)

	genefusion_list_path = sample_path + 'AllGeneFusionReads.txt'
	genefusion_fp = open(genefusion_list_path, 'r')
	genefusion_list = []
	while True:
		line = genefusion_fp.readline()
		if not line:
			break
		line = line.strip()
		if line == '':
			continue
		genefusion_list.append(line)

	genefusion_list = set(genefusion_list)

	has_cirRna = []
	has_genefusion = []

	while True:
		line = samfile.readline()
		if not line:
			break
		if not line.startswith('HISEQ:'):
			samout.write(line)
		else:
			read_id = line.split('\t')[0]
			read_id = read_id.split('/')[0]
			

			read_in_cir = False
			read_in_fusion = False

			if read_id in cirRna_list:
				read_in_cir = True
			if read_id in genefusion_list:
				read_in_fusion = True

			if read_in_cir==False and read_in_fusion==False:
				samout.write(line)
			else:
				if read_in_cir:
					has_cirRna.append(read_id)
				if read_in_fusion:
					has_genefusion.append(read_id)

	cir_fp = open(sample_path + 'reads_cirRna.txt', 'w')
	fu_fp = open(sample_path + 'reads_genefusion.txt', 'w')
	outline = '\n'.join(has_cirRna)
	cir_fp.write(outline)
	outline = '\n'.join(has_genefusion)
	fu_fp.write(outline)

	cir_fp.close()
	fu_fp.close()


	samfile.close()
	samout.close()

	#To be continued......#



