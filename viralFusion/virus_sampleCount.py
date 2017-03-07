import os
import sys
import string
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':

	if len(sys.argv) != 2:
		print sys.argv[0], 'virus_readsNum_thr'
		sys.exit()

	level_mark = 'S'
	
	samples_root_dir = '/PHShome/tw786/neurogen/Tao/BRAINCODE_output/'
	virus_thr = string.atof(sys.argv[1])
	# sample_reads_filepath = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/sampleIDtotal_sequenced_reads.txt'
	# sample_reads_dic = {}
	# sample_reads_fp = open(sample_reads_filepath, 'r')
	# for i in range(3):
	# 	line = sample_reads_fp.readline()
	# while True:
	# 	line = sample_reads_fp.readline()
	# 	if not line:
	# 		break
	# 	linesplit = line.strip().split()
	# 	sample_name = linesplit[0].strip()
	# 	sample_reads_count = string.atoi(linesplit[1].replace(',',''))
	# 	sample_reads_dic[sample_name] = sample_reads_count

	virus_count = {}
	for sample_name in os.listdir(samples_root_dir):
		if not os.path.isdir(samples_root_dir + sample_name):
			continue
		
		# if not sample_reads_dic.has_key(sample_name):
		# 	print sample_name,"don't know its total reads number"
		# 	continue
		# x_reads_num = sample_reads_dic[sample_name]
		# x_reads_per_million = float(x_reads_num)/1000000
		
		kraken_output_path = samples_root_dir + sample_name + '/kraken_output.report'

		if not os.path.isfile(kraken_output_path):
			continue
		
		kraken_output_fp = open(kraken_output_path,'r')

		eukaryota_start = False
		viruses_start = False
		bacteria_start = False
		archaea_start = False

		while True:
			line = kraken_output_fp.readline()
			if not line:
				break
			linesplit = line.split('\t')

			if linesplit[4] == '2': # taxonimic id of bacteria is 2 
				eukaryota_start = False
				viruses_start = False
				bacteria_start = True
				archaea_start = False
				continue

			if linesplit[4] == '10239': #Virus(taxonimic id = 10239)
				eukaryota_start = False
				viruses_start = True
				bacteria_start = False
				archaea_start = False
				continue
			if linesplit[4] == '2157': #archaea(taxonimic id = 2157)
				eukaryota_start = False
				viruses_start = False
				bacteria_start = False
				archaea_start = True
				continue
			if linesplit[4] == '2759': #eukaryota(taxonimic id = 2759)
				eukaryota_start = True
				viruses_start = False
				bacteria_start = False
				archaea_start = False
				continue

			if linesplit[3] == level_mark:
				name = linesplit[5].lstrip().strip()
				
				reads_num = string.atoi(linesplit[1].strip())
				species_taxid = linesplit[4]
				# ppm = round(reads_num/x_reads_per_million, 3)
				# if reads_num < reads_thr: # reads number threshold is appied here!
				# 	continue
				if viruses_start and reads_num >= virus_thr:
					if 'phage' in name:
						continue
					if not virus_count.has_key(name):
						virus_count[name] = [sample_name]
					else:
						virus_count[name].append(sample_name)
					# log_line = sample_name + '\t' + species_taxid + '\t' + str(ppm) + '\t' + str(virus_ppm_thr) + '\n'
					# logfile.write(log_line)
					# commandline = "bsub -q big %s %s %s" % (core_script, sample_name, species_taxid)
					# os.system(commandline)
					
				# if bacteria_start and ppm >= bacteria_ppm_thr:
				# 	bacteria_dic[name] = ppm
				# 	bacteria_list.append(name)
				
				# if eukaryota_start and ppm >= eukaryota_ppm_thr:
				# 	eukaryota_dic[name] = ppm
				# 	eukaryota_list.append(name)

				# if archaea_start and ppm >= archaea_ppm_thr:
				# 	archaea_dic[name] =  ppm
				# 	archaea_list.append(name)

	one_count = 0
	more_than_one_count = 0
	virus_list = []
	sample_num_list = []
	for virus in virus_count:
		virus_list.append(virus)
		sample_num_list.append(len(virus_count[virus]))
	sample_num_list, virus_list = (list(x) for x in zip(*sorted(zip(sample_num_list, virus_list))))
	plt.figure(tight_layout=True)
	ind = np.arange(len(virus_list))
	width = 0.35
	plt.bar(ind, sample_num_list, width=width, color='r')
	plt.xticks(ind+width/2., virus_list,rotation=270, fontsize=5)
	plt.xlim(0,len(virus_list)+1)
	# plt.axhline(y = 1, color = 'b')
	plt.ylabel('Number of samples')
	plt.title("#Samples detected viruses")
	savepath = "//PHShome/tw786/#Samples_detect_viruses.pdf"
	plt.savefig(savepath)

	for virus in virus_count:
		if len(virus_count[virus]) == 1:
			one_count += 1
		elif len(virus_count[virus]) > 1:
			more_than_one_count += 1
	print "There are %d virus appear only in one sample " % (one_count)
	print "There are %d virus appear more than one sample" % (more_than_one_count)
	
	for virus in virus_count:
		if 'phage' not in virus:
			print virus,'\t', len(virus_count[virus])

		