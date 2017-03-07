# Based on heatmap.R
# columns are clustered by subject ID


library(pheatmap)
# Read the table into dataframe
table = read.csv(file = '/Users/Tao/BrainMicrobiome/Tables/gtex_table_subjects_reads_1.xls', sep = '\t', header = FALSE, stringsAsFactors = FALSE)

# Transform the data frame into a matrix. 
# The colnames are composed by "subjectID brainRegion"
# The rownames are viruses names
table_mat = data.matrix(table[3:nrow(table),2:ncol(table)])
rownames(table_mat) = table[3:nrow(table),1]
cnames = paste(table[1,],table[2,])
colnames(table_mat) = cnames[-1]

# Shrink the source matrix to make sure each row and column contain at least one number greater than threshold
table_mat_selected = table_mat[,apply(table_mat,MARGIN = 2, function(x) any(x>=10))]
table_mat_selected = table_mat_selected[apply(table_mat_selected, MARGIN = 1, function(x) any(x>=10)),]
table_mat_selected[table_mat_selected>50] = 60 # for purpose of drawing legend

# Sort columns of table_mat_selected by its subjects' ID
table_mat_selected_colnames = colnames(table_mat_selected)
table_mat_selected = table_mat_selected[,order(table_mat_selected_colnames)]

# Set annotation for columns
table_mat_selected_colnames = colnames(table_mat_selected)
table_mat_selected_colnames_mat = do.call(rbind, strsplit(table_mat_selected_colnames,' '))
annotation_col = data.frame(Subjects = table_mat_selected_colnames_mat[,1])
rownames(annotation_col) = table_mat_selected_colnames

# Set annotation for rows
annotation_table = read.table(file = '/Users/Tao/BrainMicrobiome/Tables/virusesAnno.txt',
                              sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(annotation_table) = annotation_table[,1]
annotation_table = annotation_table[,-1]
annotation_table[annotation_table==""] = 'NA'
annotation_row = annotation_table

# Sort rows of table_mat_selected
table_mat_selected = table_mat_selected[rownames(annotation_row),]
# Set annotation colors
ann_colors = list(
  Taxonomy1 = c('Herpesviridae, Family' = 'forestgreen', 'Papillomaviridae, Family' = "lightpink", 
                'Gammaretrovirus, Genus' = "dodgerblue4", 'Coronavirinae, Subfamily'= "burlywood1", 
                'Bunyaviridae, Family' = "burlywood3", 'Bromoviridae, Family' = "lightpink3", 
                'Baculoviridae, Family' ="lightpink2", 'Tobamovirus, Genus' = "deeppink", 
                'Potexvirus, Genus' = "goldenrod4", 'Picornavirales, Order' = "dodgerblue1",
                'Parvoviridae, Family' = "deepskyblue", 'NA' = "white"),
  Taxonomy2 = c('Roseolovirus, Genus' = "red",  'Simplexvirus, Genus' = 'blue', 
                'Varicellovirus, Genus' = 'yellow','NA' = "white")
)

# Set legend breaks and colors
my_breaks = c(0,10,20,30,40,50,max(table_mat_selected))
my_palette = colorRampPalette(c('white', "cyan","green","orange","yellow","red"))(6) 
gtex_table = read.csv(file = '/Users/Tao/GTEx_Data_V6_Annotations_SampleAttributesDS.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Draw Heatmap
pheatmap(table_mat_selected,color = my_palette, breaks = my_breaks, 
         fontsize_row = 5,fontsize_col = 5, border_color = 'white', 
         legend_breaks = c(10,20,30,40,50,60), legend_labels = c(10,20,30,40,50,'50+'),
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = annotation_row,  annotation_col = annotation_col, fontsize = 5,
         annotation_colors = ann_colors, main = "GTEx(reads>=10)")

