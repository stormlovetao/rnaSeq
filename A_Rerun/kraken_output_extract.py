# kraken_output_extract.py  kraken_output_M kraken_report_S.txt
"""
Extract records from kraken_output based on a list of taxid
output to stdout
"""
import os
import sys
import string

if len(sys.argv) != 3:
	print "Error input"
	print "python " + sys.argv[0] + " kraken_output  taxid_list" 

kraken_output_path = sys.argv[1]
taxid_list_path = sys.argv[2]

kraken_output_fp = open(kraken_output_path)
taxid_fp = open(taxid_list_path)

taxid_list = []
while True:
	line = taxid_fp.readline()
	if not line:
		break
	line = line.strip()
	taxid_list.append(line)

while True:
	line = kraken_output_fp.readline()
	if not line:
		break
	line = line.strip()
	linesplit = line.split('\t')
	taxid = linesplit[2]
	if taxid in taxid_list:
		print line
