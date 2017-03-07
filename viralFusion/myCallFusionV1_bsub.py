import os
import sys
import string


if __name__ == '__main__':

	if len(sys.argv) != 3:
		print sys.argv[0], 'virus_ppm_thr', '1/0(1: logoverwrite; 0: lognotoverwrite)'
		sys.exit()

	level_mark = 'S'
	log_flag = sys.argv[2]
	if log_flag != '1' and log_flag != '0':
		print "can not tell log_flag"
		print sys.argv[0], 'virus_ppm_thr', '1/0(1: logoverwrite; 0: lognotoverwrite)'
		sys.exit()
	samples_root_dir = '/PHShome/tw786/neurogen/Tao/kraken_output/'
	virus_ppm_thr = string.atof(sys.argv[1])
	sample_reads_filepath = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/sampleIDtotal_sequenced_reads.txt'
	sample_reads_dic = {}
	sample_reads_fp = open(sample_reads_filepath, 'r')
	for i in range(3):
		line = sample_reads_fp.readline()
	while True:
		line = sample_reads_fp.readline()
		if not line:
			break
		linesplit = line.strip().split()
		sample_name = linesplit[0].strip()
		sample_reads_count = string.atoi(linesplit[1].replace(',',''))
		sample_reads_dic[sample_name] = sample_reads_count


	logfile_path = '/PHShome/tw786/neurogen/Tao/kraken_output_tables/myCallFusionV1_bsub.log'
	task_done_dic = {}
	if os.path.isfile(logfile_path):
		logfile_read = open(logfile_path, 'r')
		while True:
			line = logfile_read.readline()
			if not line:
				break
			linesplit = line.strip().split('\t')
			sample_name = linesplit[0]
			virus_taxid = linesplit[1]
			if not task_done_dic.has_key(sample_name):
				task_done_dic[sample_name] = [virus_taxid]
			else:
				task_done_dic[sample_name].append(virus_taxid)

	if log_flag == '0':
		logfile = open(logfile_path, 'a')
	elif log_flag == '1':
		logfile = open(logfile_path, 'w')
	core_script = "bash /PHShome/tw786/MyOwnScript/myCallFusionV1.sh"
	for sample_name in os.listdir(samples_root_dir):
		if not os.path.isdir(samples_root_dir + sample_name):
			continue
		if not sample_reads_dic.has_key(sample_name):
			print sample_name,"don't know its total reads number"
			continue
		x_reads_num = sample_reads_dic[sample_name]
		x_reads_per_million = float(x_reads_num)/1000000
		
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
				if name == 'Enterobacteria phage phiX174 sensu lato':
					continue
				reads_num = string.atoi(linesplit[1].strip())
				species_taxid = linesplit[4]
				ppm = round(reads_num/x_reads_per_million, 3)
				# if reads_num < reads_thr: # reads number threshold is appied here!
				# 	continue
				if viruses_start and ppm >= virus_ppm_thr:
					if 'phage' in name:
						continue
					log_line = sample_name + '\t' + species_taxid + '\t' + str(ppm) + '\t' + str(virus_ppm_thr) + '\n'
					logfile.write(log_line)
					commandline = "bsub -q big %s %s %s" % (core_script, sample_name, species_taxid)
					os.system(commandline)
					
				# if bacteria_start and ppm >= bacteria_ppm_thr:
				# 	bacteria_dic[name] = ppm
				# 	bacteria_list.append(name)
				
				# if eukaryota_start and ppm >= eukaryota_ppm_thr:
				# 	eukaryota_dic[name] = ppm
				# 	eukaryota_list.append(name)

				# if archaea_start and ppm >= archaea_ppm_thr:
				# 	archaea_dic[name] =  ppm
				# 	archaea_list.append(name)

		
		