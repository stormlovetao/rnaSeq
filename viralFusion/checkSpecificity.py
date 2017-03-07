# This script is needed by myCallFusion.sh
# Python script to check if a read is specificly mixed by human part and virus part
# input: output of blastn, ${homo_sapiens_taxid}_${virus_taxid}.mixedReadsID.fa.blastnout
# output:  ${homo_sapiens_taxid}_${virus_taxid}.specificMixedReads.fa
import os
import sys
import string

# command = "export PATH=~/anaconda_ete/bin:$PATH"
# os.system(command)
from ete3 import NCBITaxa

homo_sapiens_taxid = "9606"
virus_taxid = "113553"

homo_sapiens_taxid = string.atoi(homo_sapiens_taxid)
virus_taxid = string.atoi(virus_taxid)

ncbi = NCBITaxa()
homo_sapiens_taxid_lineage = ncbi.get_lineage(homo_sapiens_taxid) # return a list of taxid(integer)
virus_taxid_lineage = ncbi.get_lineage(virus_taxid)

blastn_outfile = "/PHShome/tw786/neurogen/Tao/kraken_output/HC_NZ-H152_MCPY_4_rep1/test"
outputfile = "/PHShome/tw786/neurogen/Tao/kraken_output/HC_NZ-H152_MCPY_4_rep1/test.output"
output = open(outputfile, 'w')

gi_taxid_nucl_path = "/PHShome/tw786/neurogen/Tao/taxonomy/gi_taxid_nucl.dmp"
gi_taxid_nucl = open(gi_taxid_nucl_path,'r')

gi2taxid_dic = {}
while True:
	line = gi_taxid_nucl.readline()
	if not line:
		break
	linesplit = line.strip().split('\t')
	gi_id = string.atoi(linesplit[0])
	taxid = string.atoi(linesplit[1])
	if gi2taxid_dic.has_key(gi_id):
		print(line, 'repeat?')
	gi2taxid_dic[gi_id] =  taxid

blastnout = open(blastn_outfile, 'r')


taxid_list = []
readID = "NNN"
while True:
	line = blastnout.readline()
	if not line:
		break
	if line.startswith('# BLASTN') and len(taxid_list) != 0:
		taxid_list = list(set(taxid_list))
		if (homo_sapiens_taxid not in taxid_list) or (virus_taxid not in taxid_list):
			taxid_list = []
			continue
		else:
			flag = True
			for eachid in taxid_list:
				if (eachid not in homo_sapiens_taxid_lineage) and (eachid not in virus_taxid_lineage):
					flag = False
					break
			if flag:
				output.write(readID + '\n')
				taxid_list = []
		
		
	if line.startswith('# Query:'):
		linesplit = line.strip().split(': ')
		readID = linesplit[1]
		continue
	if line.startswith(readID):
		linesplit = line.split('\t')
		gi_part = linesplit[1]
		gi_part_split = gi_part.split('|')
		gi_number = string.atoi(gi_part_split[1])
		if not gi2taxid_dic.has_key(gi_number):
			print "can not find taxid of ", gi_number
			continue
		else:
			taxid = gi2taxid_dic[gi_number]
			taxid_list.append(taxid)

output.close()






















	