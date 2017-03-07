"""
For GTEx data:
Plot good unmapped reads number v.s. total unmapped reads number
Plot good unmapped reads number v.s. total library reads number
Plot total unmapped reads number v.s. total library reads number

For BRAINCODE data:
Plot good unmapped reads number v.s. total unmapped reads number
Plot good unmapped reads number v.s. total library reads number
Plot total unmapped reads number v.s. total library reads number
"""

import os
import sys
import string
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

########################## For BRAINCODE data#########################
if True:
	gtex_file_path = "/Users/Tao/BrainMicrobiome/Figures/BRAINCODE_allSamples_reads_statistics.txt"
	gtex_fp = open(gtex_file_path)

	header_line = gtex_fp.readline()
	unmapped_good_qc_num_list = []
	total_unmapped_list = []
	total_library_reads_list = []
	while True:
		line = gtex_fp.readline()
		if not line:
			break

		linesplit= line.strip().split('\t')
		unmapped_good_qc_num = string.atoi(linesplit[2])
		total_unmapped_number = string.atoi(linesplit[3])
		total_raw_reads = string.atoi(linesplit[4])
		unmapped_good_qc_num_list.append(unmapped_good_qc_num)
		total_unmapped_list.append(total_unmapped_number)
		total_library_reads_list.append(total_raw_reads)

	np.array(unmapped_good_qc_num_list)
	np.array(total_unmapped_list)
	fig = plt.figure(tight_layout=True, figsize=(10,15))
	ax1 = fig.add_subplot(311)
	scaled_total_unmapped_list = [np.log2(x) for x in total_unmapped_list]
	scaled_good_qc_num_list = [np.log2(x) for x in unmapped_good_qc_num_list]
	ax1.scatter(scaled_total_unmapped_list, scaled_good_qc_num_list)
	# #save_path = "/home/tw83/gtex_plot_goodqcVStotalunmapp.pdf"
	# cal trend line
	z = np.polyfit(scaled_total_unmapped_list, scaled_good_qc_num_list, 1)
	p = np.poly1d(z)
	plt.plot(scaled_total_unmapped_list, p(scaled_total_unmapped_list), "r-")
	print "y=%.6fx+(%.6f)"%(z[0],z[1])

	ax1.set_xlabel("Unmapped Reads Number(log2 scaled)")
	ax1.set_ylabel("QC-passed Unmapped Reads Number(log2 scaled)")

	scaled_total_unmapped_Reads_list = [np.log2(x) for x in total_unmapped_list]
	ax2 = fig.add_subplot(312)
	scaled_total_library_Reads_list = [np.log2(x) for x in total_library_reads_list]
	ax2.scatter(scaled_total_library_Reads_list, scaled_total_unmapped_Reads_list)
	ax2.set_xlabel("Total library Reads Number(log2 scaled)")
	ax2.set_ylabel("Unmapped Reads Number(log2 scaled)")
	#ax2.set_xlim(left = 28)

	ax3 = fig.add_subplot(313)
	ax3.scatter(scaled_total_library_Reads_list, scaled_good_qc_num_list)
	ax3.set_xlabel("Total library Reads Number(log2 scaled)")
	ax3.set_ylabel("QC-passed Unmapped Reads Number(log2 scaled)")

	plt.suptitle("BRAINCODE(%d samples)" % len(unmapped_good_qc_num_list))
	save_path = "/Users/Tao/BrainMicrobiome/Figures/BRAINCODE_plot_readsplot.pdf"
	plt.savefig(save_path)

########################For GTEx Data##############################
if False:
	gtex_file_path = "/Users/Tao/BrainMicrobiome/Figures/GTEx_allSamples_reads_statistics.txt"
	gtex_fp = open(gtex_file_path)

	header_line = gtex_fp.readline()
	good_qc_num_list = []
	total_unmapped_list = []
	total_library_Mbase_list = []
	while True:
		line = gtex_fp.readline()
		if not line:
			break

		linesplit= line.strip().split('\t')
		good_qc_num = string.atoi(linesplit[2])
		total_unmapped_number = string.atoi(linesplit[3])
		total_libray_Mbase = string.atoi(linesplit[4])
		good_qc_num_list.append(good_qc_num)
		total_unmapped_list.append(total_unmapped_number)
		total_library_Mbase_list.append(total_libray_Mbase)
	# print good_qc_num_list
	# print total_unmapped_list
	print len(good_qc_num_list)
	print len(total_unmapped_list)
	np.array(good_qc_num_list)
	np.array(total_unmapped_list)
	fig = plt.figure(tight_layout=True, figsize=(10,15))
	ax1 = fig.add_subplot(311)
	scaled_total_unmapped_list = [np.log2(x) for x in total_unmapped_list]
	scaled_good_qc_num_list = [np.log2(x) for x in good_qc_num_list]
	ax1.scatter(scaled_total_unmapped_list, scaled_good_qc_num_list)
	# #save_path = "/home/tw83/gtex_plot_goodqcVStotalunmapp.pdf"

	ax1.set_xlabel("Unmapped Reads Number(log2 scaled)")
	ax1.set_ylabel("QC-passed Unmapped Reads Number(log2 scaled)")

	scaled_total_unmapped_Reads_list = [np.log2(x) for x in total_unmapped_list]
	ax2 = fig.add_subplot(312)
	scaled_total_library_Reads_list = [np.log2(x*1000000/76) for x in total_library_Mbase_list]#reads are alway 76bp
	ax2.scatter(scaled_total_library_Reads_list, scaled_total_unmapped_Reads_list)
	ax2.set_xlabel("Total library Reads Number(log2 scaled)")
	ax2.set_ylabel("Unmapped Reads Number(log2 scaled)")
	#ax2.set_xlim(left = 28)

	ax3 = fig.add_subplot(313)
	ax3.scatter(scaled_total_library_Reads_list, scaled_good_qc_num_list)
	ax3.set_xlabel("Total library Reads Number(log2 scaled)")
	ax3.set_ylabel("QC-passed Unmapped Reads Number(log2 scaled)")

	plt.suptitle("GTEx(%d samples)" % len(good_qc_num_list))
	save_path = "/Users/Tao/BrainMicrobiome/Figures/gtex_plot_readsplot.pdf"
	plt.savefig(save_path)


