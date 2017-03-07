fpkm=read.table("~/RnaSeq/ERV/ERV_id_fpkm.txt", sep = "\t", stringsAsFactors = FALSE)
fpkm_mat = as.matrix(fpkm[,2:ncol(fpkm)])
 
header=read.table("~/RnaSeq/ERV/ERV_id.txt", sep = "\t",stringsAsFactors = FALSE)
rownames(fpkm_mat) = fpkm[,1]
colnames(fpkm_mat) = header[1,2:ncol(header)]
test = list(HC_ERV3_1=fpkm_mat[2,2:91], PD_ERV3_1=fpkm_mat[2,120:140])

test = list(HC_ERVK13_1=fpkm_mat[9,2:91], PD_ERVK13_1=fpkm_mat[9,120:140])
boxplot(test,outline=FALSE)

t.test(fpkm_mat[9,2:91], fpkm_mat[9,120:140])
t.test(fpkm_mat[9,2:91])