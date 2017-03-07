#!/usr/bin/python


"""
Get all read IDs that support circRna from circRna_PD/*/allchr_circlist_nostr_r2.txt 

"""



import sys
import os

if __name__ == '__main__':
	
	if len(sys.argv) != 2:
		print "usage: ./getCircRna.py [sample_name]"
		sys.exit(1)

	sample_name = sys.argv[1]

	root_path = '/data/neurogen/circRNA_PD/'
	sample_path = root_path + sample_name + '/'

	if not os.path.exists(sample_path):
		print sample_path, "doesn't exist!"
		sys.exit(1)

	circRna_path = sample_path + 'allchr_circlist_nostr_r2.txt'

	fp = open(circRna_path, 'r')
	header = fp.readline()

	out_sample_path = '/data/neurogen/Tao/kraken_output_v2/' + sample_name 
	#out_sample_path = '/data/neurogen/Tao/test'
	out_file_path = out_sample_path + '/AllCircRnaList.txt'
	out_fp = open(out_file_path, 'w')

	while True:
		line = fp.readline()
		if not line:
			break
		linesplit = line.strip().split('\t')
		if len(linesplit) < 11:
			print line
			continue
		cicr_Rna_line = linesplit[10]
		cicr_Rna_list = cicr_Rna_line.split(',')
		for rna in cicr_Rna_list:
			if rna != '':
				out_fp.write(rna.strip()+'\n')
	fp.close()
	out_fp.close()
