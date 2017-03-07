# Based on heatmap.R
# columns are clustered by subject ID


library(pheatmap)
# Read the table into dataframe
table = read.csv(file = '/Users/Tao/BrainMicrobiome/Tables/gtex_table_subjects_BacteriaFamilyLevel_reads_1_withSampleID.xls', sep = '\t', header = FALSE, stringsAsFactors = FALSE)

# Transform the data frame into a matrix. 
# The colnames are composed by "subjectID brainRegion"
# The rownames are viruses names
table_mat = data.matrix(table[3:nrow(table),2:ncol(table)])
rownames(table_mat) = table[3:nrow(table),1]
cnames = paste(table[1,],table[2,])
colnames(table_mat) = cnames[-1]

# Shrink the source matrix to make sure each row and column contain at least one number greater than threshold
table_mat_selected = table_mat[,apply(table_mat,MARGIN = 2, function(x) any(x>=1000))]
table_mat_selected = table_mat_selected[apply(table_mat_selected, MARGIN = 1, function(x) any(x>=1000)),]
table_mat_selected[table_mat_selected>5000] = 6000 # for purpose of drawing legend

# Sort columns of table_mat_selected by its subjects' ID
table_mat_selected_colnames = colnames(table_mat_selected)
table_mat_selected = table_mat_selected[,order(table_mat_selected_colnames)]

# Set annotation for columns
table_mat_selected_colnames = colnames(table_mat_selected)
table_mat_selected_colnames_mat = do.call(rbind, strsplit(table_mat_selected_colnames,' '))
annotation_col = data.frame(Subjects = table_mat_selected_colnames_mat[,1])
rownames(annotation_col) = table_mat_selected_colnames


# Set legend breaks and colors
my_breaks = c(0,1000,2000,3000,4000,5000,max(table_mat_selected))
my_palette = colorRampPalette(c('white', "cyan","green","orange","yellow","red"))(6) 

# Draw Heatmap
pheatmap(table_mat_selected,color = my_palette, breaks = my_breaks, 
         fontsize_row = 5,fontsize_col = 5, border_color = 'white', 
         legend_breaks = c(1000,2000,3000,4000,5000,6000), legend_labels = c(1000,2000,3000,4000,5000,'5000+'),
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = annotation_col, fontsize = 5,
         main = "GTEx(reads>=1000)")

