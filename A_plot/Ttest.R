# Read the table into dataframe
table = read.csv(file = '/Users/Tao/BrainMicrobiome/Tables/braincode_table_samples_reads_1.xls', sep = '\t', check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
# Transform the data frame into a matrix. 
# The colnames are sample names
# The rownames are viruse names
table_mat = data.matrix(table[1:nrow(table),2:ncol(table)])
rownames(table_mat) = table[1:nrow(table),1]
colnames(table_mat) = colnames(table)[-1]