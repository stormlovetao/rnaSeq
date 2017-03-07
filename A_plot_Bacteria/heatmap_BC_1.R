# Based on heatmap_3.R
# For BrainCode data

library(pheatmap)
# Read the table into dataframe
table = read.csv(file = '/Users/Tao/BrainMicrobiome/Tables/braincode_table_samples_reads_1.xls', sep = '\t', check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
# Transform the data frame into a matrix. 
# The colnames are composed by "subjectID brainRegion"
# The rownames are viruses names
table_mat = data.matrix(table[1:nrow(table),2:ncol(table)])
rownames(table_mat) = table[1:nrow(table),1]
colnames(table_mat) = colnames(table)[-1]

# Shrink the source matrix to make sure each row and column contain at least one number greater than threshold
table_mat_selected = table_mat[,apply(table_mat,MARGIN = 2, function(x) any(x>=10))]
table_mat_selected = table_mat_selected[apply(table_mat_selected, MARGIN = 1, function(x) any(x>=10)),]
table_mat_selected[table_mat_selected>50] = 60 # for purpose of drawing legend

# Sort columns by batch ID
table_mat_selected_colnames_splitmat = do.call(rbind, strsplit(colnames(table_mat_selected),'_'))
table_mat_selected = table_mat_selected[,order(table_mat_selected_colnames_splitmat[,4])]
# Set column annotation by batch ID
annotation_col = data.frame(BATCH_ID = do.call(rbind, strsplit(colnames(table_mat_selected),'_'))[,4])
annotation_col$TYPE = do.call(rbind, strsplit(colnames(table_mat_selected),'_'))[,3]
rownames(annotation_col) = colnames(table_mat_selected)


# Set annotation for rows
annotation_table = read.table(file = '/Users/Tao/BrainMicrobiome/Tables/BRAINCODE_virusesAnno.txt',
                              sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(annotation_table) = annotation_table[,1]
annotation_table = annotation_table[,-1]
annotation_table[annotation_table==""] = 'NA'
annotation_row = annotation_table
# Sort rows of table_mat_selected
table_mat_selected = table_mat_selected[rownames(annotation_row),] # sort by solid order

# Set legend breaks and colors
my_breaks = c(0,10,20,30,40,50,max(table_mat_selected))
my_palette = colorRampPalette(c('white', "cyan","green","orange","yellow","red"))(6) 

pheatmap(table_mat_selected, color = my_palette, breaks = my_breaks,
         fontsize_row = 5, fontsize_col = 5, border_color = 'white',
         legend_breaks = c(10,20,30,40,50,60), legend_labels = c(10,20,30,40,50,'50+'),
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = annotation_col, annotation_row = annotation_row,
         fontsize = 5, main = "BRAINCODE(reads>=10)")