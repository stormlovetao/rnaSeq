"""
After filtering Low-Complexity reads from snap_output.bam(using filterLC.sh)
, a snap_output.good.fastq_path file will be generated.
I want to cut these good reads out from snap_output.bam to form a snap_output_good.bam
And then visulize it using IGV.
"""

import sys
import os

if __name__ == '__main__':
	import argparse
	cmd_parser = argparse.ArgumentParser()

	cmd_parser.add_argument('-p', action = 'store', dest = "folderpath",
							help = "Essential")

	result = cmd_parser.parse_args()

	if result.folderpath == None:
		print "No folderpath exist, use argument '-h' for help!"
		sys.exit(1)
	else:
		folderpath = result.folderpath

	#folderpath = "/PHShome/tw786/neurogen/Tao/run_output/PD_BN05-17_SNDA_5_rep1/Human_coronavirus_229E_uid14913"
	if folderpath.endswith('/'):
		folderpath = folderpath[:-1]
	fastq_path = folderpath + "/snap_output.good.fastq"
	if not os.path.isfile(fastq_path):
		print fastq_path, "doesn't exist, exit"
		sys.exit(1)

	sam_path = folderpath + "/snap_output.sam"
	if not os.path.isfile(sam_path):
		print sam_path, "doesn't exist, exit"
		sys.exit(1)

	fastq_fp = open(fastq_path)

	fastq_list = []
	while True:
		line = fastq_fp.readline()
		if not line:
			break
		if line.startswith('@'):
			line = line.strip()
			fastq_list.append(line[1:])
	fastq_fp.close()

	sam_fp = open(sam_path)

	sam_outpath = folderpath + '/snap_output_good.sam'
	samout = open(sam_outpath, 'w')

	while True:
		line = sam_fp.readline()
		if not line:
			break
		if line.startswith('@'):
			samout.write(line)

		if line.startswith('HISEQ'):
			linesplit1=line.split('\t')[0]
			if linesplit1 in fastq_list:
				samout.write(line)
	samout.close()







