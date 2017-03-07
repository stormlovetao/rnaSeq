"""
Call by BC_rerun_blast.sh
Input: kraken_output_M_S, kraken_output_M_S.blastout
Output: kraken_output_M_S_CF
Description: read assigned taxon by kraken, is confirmed by blast against nt database. If the taxon assigned by kraken is 
			 the best/among the most significant hits (according to e-value) found by blastn(megablast).
"""
import os
import sys

if len(sys.argv) != 4:
	print "Error input"
	print "python " + sys.argv[0] + " kraken_output_M_S  kraken_output_M_S.blastout output_path" 

kraken_output_M_S_path = sys.argv[1]
kraken_output_M_S_blastout_path = sys.argv[2]
output_path = sys.argv[3]

kraken_fp = open(kraken_output_M_S_path)
kraken_fp_start = kraken_fp.tell()
blast_fp = open(kraken_output_M_S_blastout_path)

# read kraken_output_M_S, make a dict
reads2taxid_dic = {}
while True:
	line = kraken_fp.readline()
	if not line:
		break
	linesplit = line.split('\t')
	read_id = linesplit[1]
	taxon_id = linesplit[2]
	reads2taxid_dic[read_id] = taxon_id

current_read_id = ""
confirm_reads_list = []
while True:
	line = blast_fp.readline()
	if not line:
		break
	line = line.strip()
	linesplit = line.split('\t')
	taxon_id = linesplit[3]
	read_id = linesplit[0]
	if read_id == current_read_id:
		continue
	else:
		current_read_id = read_id
		if taxon_id == reads2taxid_dic[read_id]:
			confirm_reads_list.append(read_id)

confirm_reads_list = set(confirm_reads_list)
out_fp = open(output_path, 'w')
kraken_fp.seek(kraken_fp_start)
while True:
	line = kraken_fp.readline()
	if not line:
		break
	linesplit = line.split('\t')
	read_id = linesplit[1]
	if read_id in confirm_reads_list:
		out_fp.write(line)
kraken_fp.close()
blast_fp.close()
out_fp.close()

