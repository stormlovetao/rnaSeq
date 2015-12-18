################################
# This script abstracts reads from kraken_output which are labeled to Human and with a kind of pattern
# e.g. C       HISEQ:191:C51T5ACXX:1:1307:17805:35087/1        9606    75      9606:12 0:26 9606:7
################################

import os
import sys
import string

def usage():
	print "usage: python %s sample_name, exit!" %(sys.argv[0])
	exit(1)

if __name__ == '__main__':
	
	if len(sys.argv) != 2:
		usage()
	else:
		sample_name = sys.argv[1]

	kraken_output_root_path = "/PHShome/tw786/neurogen/Tao/kraken_output/"
	kraken_output_folder_path = kraken_output_root_path + sys.argv[1]
	kraken_output_path = kraken_output_folder_path + '/kraken_output'

	if not os.path.isfile(kraken_output_path):
		print kraken_output_path, 'not exist, exit!'
		exit(1)

	fp = open(kraken_output_path)

	overlap_pattern_seqs = []
	non_overlap_pattern_seqs = []

	while True:
		line = fp.readline()
		if not line:
			break
		linesplit = line.strip().split('\t')
		if len(linesplit) != 5:
			print linesplit, 'error!'
			exit(1)
		seq_name = linesplit[1]
		tax_id = linesplit[2]
		seq_len = string.atoi(linesplit[3])
		map_pattern = linesplit[4]

		if tax_id != '9606':
			continue

		map_pattern_split = map_pattern.split()
		if len(map_pattern_split) != 3:
			continue
		# pattern 1
		if  map_pattern_split[0].startswith('9606:') and map_pattern_split[1].startswith('0:') and map_pattern_split[2].startswith('9606:'):

			n = seq_left_count = string.atoi(map_pattern_split[0].split(':')[1])
			m = seq_right_count = string.atoi(map_pattern_split[2].split(':')[1])

			pos1 = 31 + n - 1 # to this read/seq, the segment from positon 1 to pos1 can be aligned to hg38 reference 

			pos2 = seq_len - m - 29 # seq_len - 31 + 1 - m + 1, the segment from pos2 to the right end of this read can be aligned to hg38 

			if pos2 <= pos1:
				overlap_pattern_seqs.append(line)
			else:
				non_overlap_pattern_seqs.append(line)

	#output to file
	fp_out1 = open(kraken_output_folder_path + '/kraken_output.9606_Homo_sapiens.pattern1.overlap', 'w')
	fp_out2 = open(kraken_output_folder_path + '/kraken_output.9606_Homo_sapiens.pattern1.nonoverlap', 'w')
	fp_out1.writelines(overlap_pattern_seqs)
	fp_out2.writelines(non_overlap_pattern_seqs)










