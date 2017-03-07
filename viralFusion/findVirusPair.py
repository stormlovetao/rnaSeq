import os
import sys
"""
Read the kraken_output, and classify virus read's pair into four catogory:
a, belong to virus
b, belong to human
c, belong to unknown class
d, unsure, which means can not find it in kraken_output, we should remap it against human reference.   
Output:
a, the mixed reads' ID
b, the virus read ID whose pair belongs to human
c, the virus read ID whose pair can not be decided(unsure)
"""
# homo_sapiens_taxid="9606"
# virus_taxid="113553"
# kraken_output_path = "/PHShome/tw786/neurogen/Tao/kraken_output/HC_NZ-H152_MCPY_4_rep1/kraken_output"

if len(sys.argv) < 5:
		print "python findVirusPair.py [homo_sapiens_taxid] [virus_taxid] [kraken_output_path] [output_path]"
		sys.exit()
module_name = sys.argv[0]
homo_sapiens_taxid = sys.argv[1]
virus_taxid = sys.argv[2]
kraken_output_path = sys.argv[3]
output_path = sys.argv[4]

def whichCategory(readID, readtag):
	"""
	Check which catogory a read belongs to
	"""
	read_name = readID + '/' + readtag
	if read_name in unclassy_reads_set:
		return "virus_unclassy"
	elif read_name in human_reads_set:
		return "virus_human"
	elif read_name in virus_reads_set:
		return "virus_virus"
	elif read_name in mixed_reads_set:
		return "virus_mixed"
	else:
		return "virus_unsure"


kraken_output = open(kraken_output_path,'r')
virus_reads_set = []
human_reads_set = []
mixed_reads_set = []
unclassy_reads_set = []
# read kraken_output
while True:
	line = kraken_output.readline()
	if not line:
		break
	linesplit = line.strip().split('\t')
	
	read_name = linesplit[1]
	read_taxid = linesplit[2]
	read_pattern = linesplit[4]
	if linesplit[0] == 'U':
		unclassy_reads_set.append(read_name)
	if (homo_sapiens_taxid in read_pattern) and (virus_taxid in read_pattern):
		mixed_reads_set.append(read_name)
	else:
		if read_taxid == virus_taxid:
			virus_reads_set.append(read_name)
		elif read_taxid == homo_sapiens_taxid:
			human_reads_set.append(read_name)
virus_reads_set = set(virus_reads_set)
human_reads_set = set(human_reads_set)
mixed_reads_set = set(mixed_reads_set)
unclassy_reads_set = set(unclassy_reads_set)


# output
folder_name = 'virus_pair_catogories'
output_dir = output_path + '/' + folder_name
if not os.path.isdir(output_dir):
	os.mkdir(output_dir)

mixed_reads_catogory = open(output_dir+'/mixed_reads.txt', 'w')
mixed_reads_catogory_Rs = open(output_dir+'/mixed_reads_Rs.txt', 'w')
virus_human_catogory = open(output_dir+'/virus_human_reads.txt' , 'w')
virus_virus_catogory = open(output_dir+'/virus_virus_reads.txt', 'w')
virus_unclassify_catogory = open(output_dir+'/virus_unclassify_reads.txt', 'w')
virus_unsure_catogory_left = open(output_dir+'/virus_unsure_reads_left.txt', 'w')
virus_unsure_catogory_right = open(output_dir+'/virus_unsure_reads_right.txt', 'w')

virus_virus_record = []
for read_name in virus_reads_set:
	linesplit = read_name.split('/')
	readID = linesplit[0]
	readtag = linesplit[1]
	flag = 0
	if readtag == '1':
		flag = whichCategory(readID, '2')	
	if readtag == '2':
		flag = whichCategory(readID, '1')
	if flag == "virus_virus":
		if readID not in virus_virus_record:
			virus_virus_record.append(readID)
			virus_virus_catogory.write(readID + '\n')
	elif flag == "virus_human":
		virus_human_catogory.write(readID + '\n')
	elif flag == "virus_unclassy":
		virus_unclassify_catogory.write(readID + '\n')
	elif flag == "virus_unsure":
		if readtag == '1':
			virus_unsure_catogory_right.write(readID + '\n')
		elif readtag == '2':
			virus_unsure_catogory_left.write(readID + '\n')

# output all the mixed reads
mixed_reads_record = []
for read_name in mixed_reads_set:
	mixed_reads_catogory_Rs.write(read_name + '\n')
	read_split = read_name.split('/')
	readID = read_split[0]
	if readID not in mixed_reads_record:
		mixed_reads_catogory.write(readID + '\n')
		mixed_reads_record.append(readID)
	

mixed_reads_catogory.close()
mixed_reads_catogory_Rs.close()
virus_human_catogory.close()
virus_virus_catogory.close()
virus_unclassify_catogory.close()
virus_unsure_catogory_left.close()
virus_unsure_catogory_right.close()





