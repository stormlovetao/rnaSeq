
fp = open("./Align2ViralGenFilterSnapRunout.py")

while True:
	line = fp.readline()
	if not line:
		break
	if line.startswith("Total Reads"):
		#print line.split("\t")
		print line.split("\t")[0].split("   "), line.split("\t")[1].split("  ")
		line = fp.readline()
		#print line.split("\t")
		print line.split("\t")[0].split("   "), line.split("\t")[1].split("  ")