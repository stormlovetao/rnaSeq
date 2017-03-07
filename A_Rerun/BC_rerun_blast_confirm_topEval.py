"""
Call by BC_rerun_blast.sh
Input: kraken_output_M_S, kraken_output_M_S.blastout
Output: kraken_output_M_S_CF
Description: read assigned taxon by kraken, is confirmed by blast against nt database. If the taxon assigned by kraken is 
			 the best/among the most significant hits (according to e-value) found by blastn(megablast).
"""
import os
import sys
import string

if len(sys.argv) != 4:
	print "Error input"
	print "python " + sys.argv[0] + " kraken_output_M_S  kraken_output_M_S.blastout output_path" 

kraken_output_M_S_path = sys.argv[1]
kraken_output_M_S_blastout_path = sys.argv[2]
output_path = sys.argv[3]

kraken_fp = open(kraken_output_M_S_path)
kraken_fp_start = kraken_fp.tell()
blast_fp = open(kraken_output_M_S_blastout_path)
out_fp = open(output_path, 'w')

blast_reads_taxid_dic = {}
best_eval_dic = {}
while True:
	line = blast_fp.readline()
	if not line:
		break
	linesplit = line.split('\t')
	taxon_id = linesplit[3]
	read_id = linesplit[0]
	e_value = linesplit[13].strip()
	e_value = string.atof(e_value)
	if blast_reads_taxid_dic.has_key(read_id):
		best_eval = best_eval_dic[read_id]
		if e_value <= best_eval:
			blast_reads_taxid_dic[read_id].append(taxon_id)
	else:
		blast_reads_taxid_dic[read_id] = [taxon_id]
		best_eval_dic[read_id] = e_value

while True:
	line = kraken_fp.readline()
	if not line:
		break
	linesplit = line.split('\t')
	read_id = linesplit[1]
	taxon_id = linesplit[2]
	if not blast_reads_taxid_dic.has_key(read_id): # blast cannot hit any target
		continue
	else:
		if taxon_id in blast_reads_taxid_dic[read_id]:
			out_fp.write(line)
			

