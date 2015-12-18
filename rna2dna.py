
import sys


if __name__ == '__main__':
	if len(sys.argv) < 3:
		print "python rna2dna.py [rna fasta path] [dna fasta path]"
		sys.exit()
	rna_path = sys.argv[1]
	dna_path = sys.argv[2]
	print rna_path, dna_path
	rna_fp = open(rna_path)
	dna_fp = open(dna_path,'w')

	while True:
		line = rna_fp.readline()
		if not line:
			break
		if line.startswith('>'):
			dna_fp.write(line)
			continue
		else:
			line = line.upper()
			if 'U' in line:
				newline = line.replace('U' , 'T')
			dna_fp.write(newline)

	dna_fp.close()
	rna_fp.close()


