#!/usr/bin/python

# after reads aligned against hg38, we make a statistic that how many reads are mapped to chr, chr_alt, chr_random...

# usage: ./chr_reads_stat.py [sample_name] 

import sys
import os

if __name__ == '__main__':
	
	if len(sys.argv) != 2:
		print "usage: ./chr_reads_stat.py [sample_name]"
		sys.exit(1)
	
	sample_root_path = "/PHShome/tw786/neurogen/Tao/kraken_output/"
	bam_file_path = sample_root_path + sys.argv[1] + '/unmapped.filterLC.filterPhiX.hg38_snap.bam'
	chrom_size_path = "./hg38.chrom.sizes"
	if not os.path.isfile(chrom_size_path):
		print chrom_size_path, "not exist! exit"
		sys.exit(1)
	if not os.path.isfile(bam_file_path):
		print bam_file_path, "not exist! exit"
		sys.exit(1)

	sam_file_path = sample_root_path + sys.argv[1] + '/unmapped.filterLC.filterPhiX.hg38_snap.sam'
	if not os.path.isfile(sam_file_path):
		command = "samtools view %s > %s" %(bam_file_path, sam_file_path)
		os.system(command)
	
	chrom_dic = {}

	fp = open(sam_file_path)
	while True:
		line = fp.readline()
		if not line:
			break
		if not line.startswith('HISEQ'):
			continue
		linesplit = line.split('\t')
		chr_name = linesplit[2]
		if not chrom_dic.has_key(chr_name):
			chrom_dic[chr_name] = 1
		else:
			chrom_dic[chr_name] += 1

	out_file_path = sample_root_path + sys.argv[1] + '/unmapped.filterLC.filterPhiX.hg38_snap.sam.stat'
	fp_out = open(out_file_path, 'w')
	for chr_name in chrom_dic:
		chr_name_easy2read = chr_name.split('__')[0]
		outline = chr_name_easy2read + '\t' + `chrom_dic[chr_name]` + '\n'
		fp_out.write(outline)




	 

