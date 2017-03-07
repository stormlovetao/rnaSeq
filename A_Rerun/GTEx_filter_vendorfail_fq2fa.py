
import fileinput
import sys
import os
vendor_fail_list_path = sys.argv[1]

#read vendor_fail_list
fp = open(vendor_fail_list_path)
vendor_fail_list = []
while True:
	line = fp.readline()
	line = line.strip()
	if not line:
		break
	vendor_fail_list.append(line)
vendor_fail_list = set(vendor_fail_list)

switch1 = 1;
line_count = 0
for line_in in sys.stdin:

	line_count += 1
	line = line_in.strip()
	if line_count%4 == 1:
		if line.startswith('@'):
			line = line[1:]
		else:
			print "error read Id", line
			sys.exit(1)
		if line in vendor_fail_list:
			switch = 0
		else:
			switch = 1
			print '>'+line
	if line_count%4 == 2 and switch == 1:
		print line

