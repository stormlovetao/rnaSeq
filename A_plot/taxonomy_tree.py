#export PATH=~/anaconda_ete/bin:$PATH
from ete3 import NCBITaxa
ncbi = NCBITaxa()

####### BRAINCODE viruses taxonomy tree ########
fp_in = open("/PHShome/tw786/localView/overview/Tree/BRAINCODE_viruses.txt")
viruses1 = fp_in.readlines()
viruses1 = [x.strip() for x in viruses1]
viruses_taxid = ncbi.get_name_translator(viruses1)
viruses_taxid = [x[0] for x in viruses_taxid.values()]
tree = ncbi.get_topology(viruses_taxid)
file_path = "/PHShome/tw786/localView/overview/Tree/BRAINCODE_viruses_tree.txt"
fp = open(file_path, 'w')
print>>fp, tree.get_ascii(attributes=['sci_name','rank'])
fp_in.close()
fp.close()
####### GTEx viruses taxonomy tree ########
fp_in = open("/PHShome/tw786/localView/overview/Tree/GTEx_viruses.txt")
viruses2 = fp_in.readlines()
viruses2 = [x.strip() for x in viruses2]
viruses_taxid = ncbi.get_name_translator(viruses2)
viruses_taxid = [x[0] for x in viruses_taxid.values()]
tree = ncbi.get_topology(viruses_taxid)
file_path = "/PHShome/tw786/localView/overview/Tree/GTEx_viruses_tree.txt"
fp = open(file_path, 'w')
print>>fp, tree.get_ascii(attributes=['sci_name','rank'])
fp_in.close()
fp.close()
####### BRAINCODE + GTEx viruses taxonomy tree ########
viruses_merge = viruses1 + viruses2
viruses_merge = list(set(viruses_merge))
viruses_taxid = ncbi.get_name_translator(viruses_merge)
viruses_taxid = [x[0] for x in viruses_taxid.values()]
tree = ncbi.get_topology(viruses_taxid)
file_path = "/PHShome/tw786/localView/overview/Tree/BRAINCODE+GTEx_viruses_tree.txt"
fp = open(file_path, 'w')
print>>fp, tree.get_ascii(attributes=['sci_name','rank'])
fp_in.close()
fp.close()