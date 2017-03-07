"""
This script is used in myCallFusionV1.sh
It aims at parsing the blast result of Trinity.fasta
For each contig, if its mapping result only/mostly contains homo_taxid and virus_taxid, then output it

"""
import os
import sys
import string
identity_threshold = 95.0
if len(sys.argv) != 4:
	print "python parseBlast.py [homo_taxid] [virus_taxid] [blast_outpath]"
	sys.exit()
from ete3 import NCBITaxa
ncbi = NCBITaxa()
homo_taxid = sys.argv[1]
virus_taxid = sys.argv[2]
blast_outpath = sys.argv[3]
homo_order_taxid = 9443
blastnout = open(blast_outpath, 'r')
output_dir = os.path.dirname(blast_outpath)
output_path = output_dir + '/parseBlast_ouput.txt'
output = open(output_path, 'w')
taxid_list = []
contig_header = "NNN"
contig_name = "NNN"
while True:
	line = blastnout.readline()
	if not line:
		break
	if line.startswith('# BLASTN') and len(taxid_list) != 0:
		#taxid_list = list(set(taxid_list))
		if (homo_taxid not in taxid_list) or (virus_taxid not in taxid_list):
			taxid_list = []
			continue
		else:
			total_hits = len(taxid_list)
			homo_count = taxid_list.count(homo_taxid)
			virus_count = taxid_list.count(virus_taxid)
			homo_family_count = 0
			for taxid in taxid_list:
				if taxid != homo_taxid and taxid != virus_taxid:
					if homo_order_taxid in ncbi.get_lineage(taxid):
						homo_family_count += 1
			if (homo_count + virus_count + homo_family_count)/float(total_hits) >= 0.95:
				output.write(contig_header + '\n')
				taxid_list = []
			# if contig_name == "c7_g5_i1":
			# 	print total_hits, homo_count, homo_family_count, virus_count
			# 	for taxid in taxid_list:
			# 		if taxid != homo_taxid and taxid != virus_taxid and (homo_order_taxid not in ncbi.get_lineage(taxid)):
			# 			print taxid, 
		
	if line.startswith('# Query:'):
		linesplit = line.strip().split(': ')
		contig_header = linesplit[1]
		contig_header_split = contig_header.split()
		contig_name = contig_header_split[0]
		continue
	if line.startswith(contig_name):
		linesplit = line.split('\t')
		taxid = linesplit[3]
		identity = string.atof(linesplit[5])
		if identity < identity_threshold:
			continue
		else:
			taxid_list.append(taxid)

output.close()












	