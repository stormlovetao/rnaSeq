"""
Given kraken_output, we want to find if reads that support cirRna and gene fusion exist in unclassified reads or human reads(9606)

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

	kraken_output_path = sample_path + 'kraken_output'

	kraken_output = open(kraken_output_path, 'r')


	circRna_list_path = sample_path + 'AllCircRnaList.txt'
	if not os.path.isfile(circRna_list_path):
		sys.exit(1)

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
	if not os.path.isfile(genefusion_list_path):
		sys.exit(1)
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

	unclassified_cirReads_list = []
	unclassified_fusReads_list = []
	human_cirReads_list = []
	human_fusReads_list = []
	while True:
		line = kraken_output.readline()
		if not line:
			break
		line = line.strip()
		linesplit = line.split('\t')
		read_id = linesplit[1].split('/')[0]
		if linesplit[0] == 'U':
			if read_id in cirRna_list:
				unclassified_cirReads_list.append(read_id)
			if read_id in genefusion_list:
				unclassified_fusReads_list.append(read_id)
		elif linesplit[0] == 'C' and linesplit[2] == '9606':
			if read_id in cirRna_list:
				human_cirReads_list.append(read_id)
			if read_id in genefusion_list:
				human_fusReads_list.append(read_id)

	unclassified_cirReads_list = set(unclassified_cirReads_list)
	unclassified_fusReads_list = set(unclassified_fusReads_list)
	human_cirReads_list = set(human_cirReads_list)
	human_fusReads_list = set(human_fusReads_list)

	# print len(unclassified_cirReads_list), len(unclassified_fusReads_list)
	# print len(human_cirReads_list), len(human_fusReads_list)
	unclassified_cirReads_fp = open(sample_path + 'kraken.unclassified.cirRna.txt', 'w')
	unclassified_fusReads_fp = open(sample_path + 'kraken.unclassified.fusion.txt', 'w')
	human_cirReads_fp = open(sample_path + 'kraken.human9606.cirRna.txt', 'w')
	human_fusReads_fp = open(sample_path + 'kraken.human9606.fusion.txt', 'w')

	unclassified_cirReads_fp.write('\n'.join(unclassified_cirReads_list))
	unclassified_cirReads_fp.close()

	unclassified_fusReads_fp.write('\n'.join(unclassified_fusReads_list))
	unclassified_fusReads_fp.close()

	human_cirReads_fp.write('\n'.join(human_cirReads_list))
	human_cirReads_fp.close()

	human_fusReads_fp.write('\n'.join(human_fusReads_list))
	human_fusReads_fp.close()



	kraken_output.close()
	cirRna_fp.close()
	genefusion_fp.close()






















