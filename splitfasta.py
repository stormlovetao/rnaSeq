import sys
if __name__ == '__main__':
	


	fp_in = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.fasta')
	fp_out1 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p1.fasta','w')
	fp_out2 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p2.fasta','w')
	fp_out3 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p3.fasta','w')
	fp_out4 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p4.fasta','w')
	fp_out5 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p5.fasta','w')
	fp_out6 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p6.fasta','w')
	fp_out7 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p7.fasta','w')
	fp_out8 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p8.fasta','w')

	fp_out9 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p9.fasta','w')
	fp_out10 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p10.fasta','w')
	fp_out11 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p11.fasta','w')
	fp_out12 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p12.fasta','w')
	fp_out13 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p13.fasta','w')
	fp_out14 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p14.fasta','w')
	fp_out15 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p15.fasta','w')
	fp_out16 = open('./SILVA_123_SSUParc_tax_silva_trunc.2dna.p16.fasta','w')

	cut_thr1 = 311557
	cut_thr2 = cut_thr1 *2
	cut_thr3 = cut_thr1 *3
	cut_thr4 = cut_thr1 *4
	cut_thr5 = cut_thr1 *5
	cut_thr6 = cut_thr1 *6
	cut_thr7 = cut_thr1 *7
	cut_thr8 = cut_thr1 *8
	cut_thr9 = cut_thr1 *9
	cut_thr10 = cut_thr1 *10
	cut_thr11 = cut_thr1 *11
	cut_thr12 = cut_thr1 *12
	cut_thr13 = cut_thr1 *13
	cut_thr14 = cut_thr1 *14
	cut_thr15 = cut_thr1 *15
	count = 0
	while True:
		line = fp_in.readline()
		if not line:
			break

		if line.startswith('>'):
			count += 1

		if count < cut_thr1:
			fp_out1.write(line)
		elif count < cut_thr2:
			fp_out2.write(line)
		elif count < cut_thr3:
			fp_out3.write(line)
		elif count < cut_thr4:
			fp_out4.write(line)
		elif count < cut_thr5:
			fp_out5.write(line)
		elif count < cut_thr6:
			fp_out6.write(line)
		elif count < cut_thr7:
			fp_out7.write(line)
		elif count < cut_thr8:
			fp_out8.write(line)
		elif count < cut_thr9:
			fp_out9.write(line)
		elif count < cut_thr10:
			fp_out10.write(line)
		elif count < cut_thr11:
			fp_out11.write(line)
		elif count < cut_thr12:
			fp_out12.write(line)
		elif count < cut_thr13:
			fp_out13.write(line)
		elif count < cut_thr14:
			fp_out14.write(line)
		elif count < cut_thr15:
			fp_out15.write(line)
		else:
			fp_out16.write(line)

	fp_out1.close()
	fp_out2.close()
	fp_out3.close()
	fp_out4.close()
	fp_out5.close()
	fp_out6.close()
	fp_out7.close()
	fp_out8.close()

	fp_out9.close()
	fp_out10.close()
	fp_out11.close()
	fp_out12.close()
	fp_out13.close()
	fp_out14.close()
	fp_out15.close()
	fp_out16.close()
		