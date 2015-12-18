"""
Creat a list of folders named by the samples' names.

"""


import os
import string
import sys


if __name__ == '__main__':
	
	# import argparse
	# cmd_parser = argparse.ArgumentParser()
	# cmd_parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	# cmd_parser.add_argument('-f', action = 'store', dest = "folder_names_path",
	# 						help = "Essential,The name of source gene expression file")

	# result = cmd_parser.parse_args()

	# if result.filename == None:
	# 	print "No gene expression source file exit, use argument '-h' for help!"
	# 	sys.exit(1)
	# else:
	# 	folder_names_path = result.folder_names_path
		
	folder_names_path = './folder_list.txt'
	fp = open(folder_names_path, 'r')

	
	while True:
		line = fp.readline().strip()
		if not line:
			break
		folder_name = line.strip()
		unmapped_bam_path = "../neurogene/rnaseq_PD/run_output/" + folder_name + "/unmapped.bam"
		if not os.path.isfile(unmapped_bam_path):
			print "Cannot find the bam path: " + unmapped_bam_path
			sys.exit(1)

		fastq_folder_path = './' + folder_name + '/'

		if not os.path.exists(fastq_folder_path):
			os.mkdir(fastq_folder_path)

		unmapped_fastq_path = fastq_folder_path + folder_name + "unmapped.fastq"
		shell_command = "bam2fastx -q -A -N -Q -o " + unmapped_fastq_path + " " + unmapped_bam_path 
		os.command(shell_command)

