#################
#
# Kraken labels ~95% unclassified reads as human reads, which can not be linearly mapped to hg38 reference genome.
# That's because only parts of such a read can be mapped to hg38 or its sub-sequences come from different locations of hg38. 
# Here we focus on the second condition: read's subsequences(here we only focus on two subsequences) come from different locations but can cover whole read.
# We name such reads as weird reads, the subsequence near 5' end as left sub-read, the other subsequence as right sub-read
#
# So, there will be two condition for a weird read:
# condition A: two parts of this read locate in same strand(plus strand or minus strand)
# condition B: two parts of this read locate in differen strands.
# To condition A, there also will be two sub-conditions:
# 	condition A-1: there is junction between two parts(this may be resulted from normal Rna splicing)
#   	To condition A-1, there is also two sub-conditions: left sub-read locates at upstream(A-1-1); right sub-read locates at upstream(A-1-2).
# 	condition A-2: there is no junciton between two parts
# 		To condition A-2, there is also two sub-conditions: right sub-read locates at upstream(A-2-1); left sub-read locates at upstream(A-2-2).
# To condition B, there are two sub-conditons:
# 	condition B-1: there is "junction" between two parts, which means there is no overlap if we map one sub-read to its oppsite strand.
# 		B-1-1: left sub-read locates at upstream; B-1-2: right sub-read locates at upstream
# 	condition B-2: there is no such "junction"
# 		B-2-1: left sub-read locates at upstream; B-2-2: right sub-read locates at upstream
#################

import os
import sys
import string

class BlatRecord:
	def __init__(self, q_name, q_size, q_start, q_end, t_name, t_start, t_end, strand):
		self.q_name = q_name
		self.q_size = q_size
		self.q_start = q_start
		self.q_end = q_end
		self.t_name = t_name
		self.t_start = t_start
		self.t_end = t_end
		self.strand = strand

pattern = {}

for x in xrange(1,10):
	pattern[x] = []

def which_pattern(left_r, right_r):
	if left_r.q_name != right_r.q_name:
		print 'error in function which_pattern, left_r.q_name != right_r.q_name, exit'
		exit(1)
	if left_r.t_name != right_r.t_name:
		print 'error in function which_pattern, left_r.t_name != right_r.t_name, exit'
		exit(1)
	r_pair = (left_r, right_r)
	count = 0
	if left_r.strand == right_r.strand and left_r.strand == 1: # both subreads are on positive strand
		if left_r.t_end < right_r.t_start:
			pattern[1].append(r_pair)
			count += 1

		elif left_r.t_start > right_r.t_end:
			pattern[2].append(r_pair)
			count += 1
		elif left_r.t_end >= right_r.t_start and left_r.t_start < right_r.t_start and left_r.t_end < right_r.t_end:
			pattern[3].append(r_pair)
			count += 1
		elif left_r.t_start <= right_r.t_end and left_r.t_start > right_r.t_start and left_r.t_end > right_r.t_end:
			pattern[4].append(r_pair) 
			count += 1
		else:
			pattern[9].append(r_pair)
			count += 1
	elif left_r.strand == right_r.strand and left_r.strand == -1: # both subreads are on negtive strand
		if left_r.t_end > right_r.t_start:
			pattern[1].append(r_pair)
			count += 1

		elif left_r.t_start < right_r.t_end:
			pattern[2].append(r_pair)
			count += 1
		elif left_r.t_end <= right_r.t_start and left_r.t_start > right_r.t_start and left_r.t_end > right_r.t_end:
			pattern[3].append(r_pair)
			count += 1
		elif left_r.t_start >= right_r.t_end and left_r.t_start < right_r.t_start and left_r.t_end < right_r.t_end:
			pattern[4].append(r_pair)
			count += 1
		else:
			pattern[9].append(r_pair)
			count += 1
	elif left_r.strand == 1 and right_r.strand == -1:
		if left_r.t_end < right_r.t_end:
			pattern[5].append(r_pair)
			count += 1
		elif left_r.t_start > right_r.t_start:
			pattern[6].append(r_pair)
			count += 1
		elif left_r.t_end >= right_r.t_end and left_r.t_start < right_r.t_end and left_r.t_end < right_r.t_start:
			pattern[7].append(r_pair)
			count += 1
		elif left_r.t_start <= right_r.t_start and left_r.t_start > right_r.t_end and left_r.t_end > right_r.t_start:
			pattern[8].append(r_pair)
			count += 1
		else:
			pattern[9].append(r_pair)
			count += 1


	elif left_r.strand == -1 and right_r.strand == 1:
		if left_r.t_end > right_r.t_end:
			pattern[5].append(r_pair)
			count += 1
		elif left_r.t_start < right_r.t_start:
			pattern[6].append(r_pair)
			count += 1
		elif left_r.t_end <= right_r.t_end and left_r.t_start > right_r.t_end and left_r.t_end > right_r.t_start:
			pattern[7].append(r_pair)
			count += 1
		elif left_r.t_start >= right_r.t_start and left_r.t_start < right_r.t_end and left_r.t_end < right_r.t_start:
			pattern[8].append(r_pair)
			count += 1
		else:
			pattern[9].append(r_pair)
			count += 1
	return count





if __name__ == '__main__':
	
	#psl_file_path = './test-input'
	psl_file_path = '/PHShome/tw786/neurogen/Tao/kraken_output/PD_BN04-42_SNDA_5_rep1/kraken_output.9606_Homo_sapiens.pattern1.overlap.psl'

	psl = open(psl_file_path,'r')

	head_line_number = 1

	for i in range(0, head_line_number):
		psl.readline()

	ReadDic = {}
	while True:
		line = psl.readline()
		if not line:
			break

		linesplit = line.split('\t')
		flag = 0
		for x in xrange(1,8): # must be perfect match
			if linesplit[x] != '0':
				flag = 1
				break
		if flag:
			continue

		strand = linesplit[8]
		if strand == '+':
			strand = 1
		elif strand == '-':
			strand = -1

		block_count = string.atoi(linesplit[17])
		if block_count != 1:
			continue
		q_name = linesplit[9]
		q_size = string.atoi(linesplit[10])
		q_start = string.atoi(linesplit[11])
		q_end  = string.atoi(linesplit[12])
		t_name = linesplit[13]
		t_start = string.atoi(linesplit[15])
		t_end = string.atoi(linesplit[16])

		if q_start != 0 and q_end != q_size:
			continue
		if q_start == 0 and q_end == q_size:
			continue

		record = BlatRecord(q_name, q_size, q_start, q_end, t_name, t_start, t_end, strand)

		if not ReadDic.has_key(q_name):
			ChrDic = {}
			ChrDic[t_name] = ([],[])
			if record.q_start == 0:
				ChrDic[t_name][0].append(record)
			elif record.q_end == record.q_size:
				ChrDic[t_name][1].append(record)

			ReadDic[q_name] = ChrDic
		else:
			ChrDic = ReadDic[q_name]

			if ChrDic.has_key(t_name):
				if record.q_start == 0:
					ChrDic[t_name][0].append(record)
				elif record.q_end == record.q_size:
					ChrDic[t_name][1].append(record)
			else:
				ChrDic[t_name] = ([],[])
				if record.q_start == 0:
					ChrDic[t_name][0].append(record)
				elif record.q_end == record.q_size:
					ChrDic[t_name][1].append(record)

	
	for read in ReadDic:
		ChrDic = ReadDic[read]
		for chrm in ChrDic:
			left_list = ChrDic[chrm][0]
			right_list = ChrDic[chrm][1]
			if len(left_list) == 0 or len(right_list) == 0:
				continue
			for left_r in left_list:
				for right_r in right_list:
					if left_r.q_end < right_r.q_start:
						# print left_r.q_name, left_r.q_start, left_r.q_end, left_r.t_name, left_r.t_start, left_r.t_end, left_r.strand
					 	# print right_r.q_name, right_r.q_start, right_r.q_end, right_r.t_name, right_r.t_start, right_r.t_end , right_r.strand 
						continue 
					pattern_num = which_pattern(left_r, right_r)
					# if pattern_num == 1:
					# 	print left_r.q_name, left_r.q_start, left_r.q_end, left_r.t_name, left_r.t_start, left_r.t_end, left_r.strand
					# 	print right_r.q_name, right_r.q_start, right_r.q_end, right_r.t_name, right_r.t_start, right_r.t_end , right_r.strand 
	for x in xrange(1,10):
		result_path = '/PHShome/tw786/neurogen/Tao/kraken_output/PD_BN04-42_SNDA_5_rep1/kraken_output.9606_Homo_sapiens.pattern1.overlap.psl.pattern%d' %(x)
		fp = open(result_path, 'w')
		header = "Q_name\tQ_size\tT_name\tStrand_1\tQ_start_1\tQ_end_1\tT_start_1\tT_end_1\tStrand_2\tQ_start_2\tQ_end_2\tT_start_2\tT_end_2\n"
		fp.write(header)
		for pair in pattern[x]:
			left_r = pair[0]
			right_r = pair[1]
			if left_r.strand == 1:
				strand1 = '+'
			elif left_r.strand == -1:
				strand1 = '-'
			if right_r.strand == 1:
				strand2 = '+'
			elif right_r.strand == -1:
				strand2 = '-'
			outline = left_r.q_name + '\t' + str(left_r.q_size) + '\t' + left_r.t_name + '\t' + strand1 + '\t' + str(left_r.q_start) + '\t' \
						+ str(left_r.q_end) + '\t' + str(left_r.t_start) + '\t' + str(left_r.t_end) + '\t' + strand2 + '\t' \
						+ str(right_r.q_start) + '\t' + str(right_r.q_end) + '\t' + str(right_r.t_start) + '\t' + str(right_r.t_end) + '\n'
			fp.write(outline)
		fp.close()
		#print x, len(pattern[x])



	# for read in ReadDic:
	# 	for chrm in ReadDic[read]:
	# 		if chrm != 'chr19':
	# 			continue
	# 		print chrm
	# 		for x in ReadDic[read][chrm][1]:
	# 			print x.q_name, x.q_start, x.q_end, x.t_start, x.t_end, x.strand

				
				
			






		









