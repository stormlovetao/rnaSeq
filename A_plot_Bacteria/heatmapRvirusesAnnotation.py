#export PATH=~/anaconda_ete/bin:$PATH
from ete3 import NCBITaxa
ncbi = NCBITaxa()
# viruses_out_path = "/Users/Tao/BrainMicrobiome/Tables/gtex_table_subjects_reads10_viruses.txt"
viruses_out_path = '/PHShome/tw786/localView/GTExTables/gtex_table_subjects_reads10_viruses.txt'

fp_in = open(viruses_out_path)
viruses = fp_in.readlines()
viruses = [x.strip() for x in viruses]
viruses_taxid = ncbi.get_name_translator(viruses)
viruses_taxid = [x[0] for x in viruses_taxid.values()]
tree = ncbi.get_topology(viruses_taxid)
file_path = "/PHShome/tw786/localView/overview/Tree/GTEx_viruses_tree.txt"
fp = open(file_path, 'w')
print>>fp, tree.get_ascii(attributes=['sci_name','rank'])
fp_in.close()
fp.close()

fp = open(file_path)