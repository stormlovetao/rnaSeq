import string
import cPickle as pickle
import matplotlib.pyplot as plt
# fp = open('./unmapped.filterLC.filterPhiX.hg38_snap.sam')
# mapq_list = []
# while True:
# 	line = fp.readline()
# 	if not line:
# 		break
# 	linesplit = line.split('\t')
# 	mapq = linesplit[4]
# 	mapq = string.atoi(mapq)
# 	mapq_list.append(mapq)
# pickle.dump(mapq_list, open('./test.pickle', 'wb'), True)

mapq_list = pickle.load(open('./test.pickle', 'rb'))

figure()

plt.hist(mapq_list)

savefig('test.pdf')

