
import os
import sys

body_site = sys.argv[1]

merged_kraken_out_dir = "/PHShome/tw786/neurogen/Tao/GTEx_output/MergedFaKraken"
merged_kraken_out_path =  merged_kraken_out_dir + '/' + body_site + "_merged_kraken_output"
if not os.path.isfile(merged_kraken_out_path):
	print merged_kraken_out_path, "doesn't exist"
	sys.exit()

merged_kraken_fp = open(merged_kraken_out_path, 'r')
split_kraken_out_dir = "/PHShome/tw786/neurogen/Tao/GTEx_output/" + body_site
if not os.path.isdir(split_kraken_out_dir):
	os.mkdir(split_kraken_out_dir)

srr_fp_dic = {}

while True:
	line = merged_kraken_fp.readline()
	if not line:
		break
	linesplit = line.split('\t')
	read_id = linesplit[1]
	read_id_split = read_id.split('_')
	srr_id = read_id_split[0]
	linesplit[1] = read_id_split[1]
	outline = '\t'.join(linesplit)
	srr_dir = split_kraken_out_dir + '/' + srr_id
	if not os.path.isdir(srr_dir):
		os.mkdir(srr_dir)
	if srr_fp_dic.has_key(srr_id):
		srr_fp_dic[srr_id].write(outline)
	else:
		fp = open(srr_dir+'/kraken_output', 'w')
		srr_fp_dic[srr_id] = fp

for srr_id in srr_fp_dic:
	srr_fp_dic[srr_id].close()



