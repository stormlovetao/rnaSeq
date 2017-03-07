# Based on heatmap_1.R
# columns are sorted by subject ID, and each subject contains all its samples
# columns are annotated by subject ID, SMGEBTCH(Genotype or Expression Batch ID, Batch when DNA/RNA from a sample was analyzed)
library(pheatmap)
# Read the table into dataframe
table = read.csv(file = '/Users/Tao/BrainMicrobiome/Tables/gtex_table_subjects_BacteriaFamilyLevel_reads_1_withSampleID.xls', sep = '\t', header = FALSE, stringsAsFactors = FALSE)

# Transform the data frame into a matrix. 
# The colnames are composed by "subjectID brainRegion"
# The rownames are viruses names
table_mat = data.matrix(table[4:nrow(table),2:ncol(table)])
rownames(table_mat) = table[4:nrow(table),1]
cnames = table[2,]
colnames(table_mat) = cnames[-1]

# Shrink the source matrix to make sure each row and column contain at least one number greater than threshold
table_mat_selected = table_mat[,apply(table_mat,MARGIN = 2, function(x) any(x>=1000))]
table_mat_selected = table_mat_selected[apply(table_mat_selected, MARGIN = 1, function(x) any(x>=1000)),]
table_mat_selected[table_mat_selected>5000] = 6000 # for purpose of drawing legend

# Sort columns of table_mat_selected by its subjects' ID
table_mat_selected_colnames = colnames(table_mat_selected)
table_mat_selected = table_mat_selected[,order(table_mat_selected_colnames)]

# Set annotation for columns, annotation = subject_id
table_mat_selected_colnames = colnames(table_mat_selected)
table_mat_selected_colnames_mat = do.call(rbind, strsplit(table_mat_selected_colnames,'-'))
table_mat_selected_colnames_mat = table_mat_selected_colnames_mat[,1:2]
table_mat_selected_colnames_anno = paste(table_mat_selected_colnames_mat[,1], table_mat_selected_colnames_mat[,2], sep = '-')
annotation_col = data.frame(Subjects = table_mat_selected_colnames_anno)
rownames(annotation_col) = table_mat_selected_colnames


# Read GTEx samples' attributes
gtex_table = read.csv(file = '/Users/Tao/BrainMicrobiome/Tables/GTEx_Data_V6_Annotations_SampleAttributesDS.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(gtex_table) = gtex_table[,1]
gtex_table = gtex_table[colnames(table_mat_selected),] # Only keep those we use

# Set column annotation = batch id
gtex_table_smgebtch = gtex_table$SMGEBTCH
gtex_table_smgebtch[gtex_table_smgebtch==""] = "NoValue"
annotation_col$SMGEBTCH = gtex_table_smgebtch
#annotation_col = data.frame(SMCENTER = gtex_table_smgebtch)
rownames(annotation_col) = colnames(table_mat_selected)

# sort table_mat_selected by the order of SMGEBTCH/SMNABTCH
# table_mat_selected = table_mat_selected[,order(gtex_table_smgebtch)]

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
