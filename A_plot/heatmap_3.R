# Based on heatmap_2.R
# columns are sorted by SMGEBTCH(Genotype or Expression Batch ID, Batch when DNA/RNA from a sample was analyzed)
# columns are annotated by SMGEBTCH
# replace SMGEBTCH to SMNABTCH,SMNABTCHD
library(pheatmap)
# Read the table into dataframe
table = read.csv(file = '/Users/Tao/BrainMicrobiome/Tables/gtex_table_subjects_reads_1_withSampleID.xls', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
# Transform the data frame into a matrix. 
# The colnames are composed by "subjectID brainRegion"
# The rownames are viruses names
table_mat = data.matrix(table[4:nrow(table),2:ncol(table)])
rownames(table_mat) = table[4:nrow(table),1]
cnames = table[2,]
colnames(table_mat) = cnames[-1]

# Shrink the source matrix to make sure each row and column contain at least one number greater than threshold
table_mat_selected = table_mat[,apply(table_mat,MARGIN = 2, function(x) any(x>10))]
table_mat_selected = table_mat_selected[apply(table_mat_selected, MARGIN = 1, function(x) any(x>10)),]
table_mat_selected[table_mat_selected>50] = 60 # for purpose of drawing legend

# filter phages Geobacillus virus E3, Lactococcus Phage ASCC191
table_mat_selected =  table_mat_selected[rownames(table_mat_selected) != "Geobacillus virus E3",]
table_mat_selected =  table_mat_selected[rownames(table_mat_selected) != "Lactococcus Phage ASCC191",]


# Set annotation for rows
# annotation_table = read.table(file = '/Users/Tao/BrainMicrobiome/Tables/gtex_virusesAnno.txt',
#                               sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# rownames(annotation_table) = annotation_table[,1]
# annotation_table = annotation_table[,-1]
# annotation_table[annotation_table==""] = 'NA'
# annotation_row = annotation_table
# # Sort rows of table_mat_selected
# table_mat_selected = table_mat_selected[rownames(annotation_row),] # sort by solid order


# Read GTEx samples' attributes
gtex_table = read.csv(file = '/Users/Tao/BrainMicrobiome/Tables/GTEx_Data_V6_Annotations_SampleAttributesDS.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(gtex_table) = gtex_table[,1]
gtex_table = gtex_table[colnames(table_mat_selected),] # Only keep those we use. rownames(gtex_table)==colnames(table_mat_selected)

# Set column annotation = batch id
gtex_table_smgebtch = gtex_table$SMGEBTCH
#gtex_table_smgebtch = gtex_table$SMCENTER
gtex_table_smgebtch[gtex_table_smgebtch==""] = "NoValue"
annotation_col = data.frame(SMGEBTCH = gtex_table_smgebtch)
# annotation_col = data.frame(SMCENTER = gtex_table_smgebtch)
rownames(annotation_col) = rownames(gtex_table) # or = colnames(table_mat_selected)

# sort table_mat_selected by the order of SMGEBTCH/SMNABTCH
table_mat_selected = table_mat_selected[,order(gtex_table_smgebtch)]




# Set annotation colors
# ann_colors = list(
#   Taxonomy1 = c('Herpesviridae, Family' = 'forestgreen', 'Papillomaviridae, Family' = "lightpink", 
#                 'Gammaretrovirus, Genus' = "dodgerblue4", 'Coronavirinae, Subfamily'= "burlywood1", 
#                 'Bunyaviridae, Family' = "burlywood3", 'Bromoviridae, Family' = "lightpink3", 
#                 'Baculoviridae, Family' ="lightpink2", 'Tobamovirus, Genus' = "deeppink", 
#                 'Potexvirus, Genus' = "goldenrod4", 'Picornavirales, Order' = "dodgerblue1",
#                 'Parvoviridae, Family' = "deepskyblue", 'NA' = "white"),
#   Taxonomy2 = c('Roseolovirus, Genus' = "red",  'Simplexvirus, Genus' = 'blue', 
#                 'Varicellovirus, Genus' = 'yellow','NA' = "white")
# )

# Set legend breaks and colors
my_breaks = c(0,10,20,30,40,50,max(table_mat_selected))
my_palette = colorRampPalette(c('white', "cyan","green","orange","yellow","red"))(6) 

table_mat_selected = table_mat_selected[,apply(table_mat_selected,MARGIN = 2, function(x) any(x>10))]
table_mat_selected = table_mat_selected[apply(table_mat_selected,MARGIN = 1, function(x) any(x>10)),]

# Draw Heatmap
pheatmap(table_mat_selected,color = my_palette, breaks = my_breaks, 
         fontsize_row = 5,fontsize_col = 5, border_color = 'white', 
         legend_breaks = c(10,20,30,40,50,60), legend_labels = c(10,20,30,40,50,'50+'),
         cluster_rows = TRUE, cluster_cols = FALSE,
         annotation_col = annotation_col,  fontsize = 5,
         annotation_colors = ann_colors, main = "GTEx(reads>10)")

