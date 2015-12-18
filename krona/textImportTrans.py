fp = open("./kraken_output.translate")
fp_out = open("./kraken_output.translate.txt", 'w')

while True:
	line = fp.readline()
	if not line:
		break
	linesplit = line.split('\t')
	line_second_part = linesplit[1].strip()
	line_second_part_split = line_second_part.split(';')
	outline = '1\t'+'\t'.join(line_second_part_split)+'\n'
	fp_out.write(outline)
