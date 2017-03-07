# Python script to map virus gi to taxid
# input: /PHShome/tw786/neurogen/Tao/kraken_standard_db/taxonomy/gi_taxid_nucl.dmp
# 		/virus_fa=/PHShome/tw786/neurogen/Tao/fna/KrakenSDBViral/virus_kraken_sdb.fa
# output: /PHShome/tw786/neurogen/Tao/fna/KrakenSDBViral/virus_kraken_sdb.fa.gi2taxid
from __future__ import print_function
import os
import sys


def warning(*objs):
	print("WARNING: ", *objs, file=sys.stderr)

gi_taxid_nucl_path = "/PHShome/tw786/neurogen/Tao/kraken_standard_db/taxonomy/gi_taxid_nucl.dmp"
gi_taxid_nucl = open(gi_taxid_nucl_path,'r')

gi2taxid_dic = {}
while True:
	line = gi_taxid_nucl.readline()
	if not line:
		break
	linesplit = line.strip().split('\t')
	gi_id = linesplit[0]
	taxid = linesplit[1]
	if gi2taxid_dic.has_key(gi_id):
		print(line, 'repeat?')
	gi2taxid_dic[gi_id] =  taxid

virus_fa_path = "/PHShome/tw786/neurogen/Tao/fna/KrakenSDBViral/virus_kraken_sdb.fa"
output_path = "/PHShome/tw786/neurogen/Tao/fna/KrakenSDBViral/virus_kraken_sdb.fa.virus2taxid"
outputfp = open(output_path, 'w')
command = "grep '>' %s " % (virus_fa_path)
command_return = os.popen(command)
while True:
	line = command_return.readline()
	if not line:
		break
	linesplit = line.strip().split('|')
	gi = linesplit[1]
	if not gi2taxid_dic.has_key(gi):
		warning('gi2taxid_dic do not have key:', gi)
		sys.exit()
	taxid = gi2taxid_dic[gi]
	outline = line.strip() + '\t' + 'taxid:' + taxid + '\n'
	outputfp(outline)



