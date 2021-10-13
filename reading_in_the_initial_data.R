rna_express <- readRDS('E:/scrnaseqscript/foxa1_data/foxa1_expressionmatrix_SE.rds')
peak_matrix <- readRDS('E:/scrnaseqscript/foxa1_data/foxa1_peakmatrix_RSE.rds')


library(rtracklayer)

#rtracklayer::export(con = 'foxa1_data/bedpeaksfoxa1_new.bed', object = peak_matrix@rowRanges, ignore.strand = TRUE, format = 'bed')
sparse_mat <- peak_matrix@assays@data$PeakMatrix

library(Matrix)

library(Seurat)


#manage the peak matrix
idx <- which((gsub("#.*","",colnames(sparse_mat)) == 'H2_multiome') == TRUE)
celltypes <- rna_express$luminal_basal[idx]
write.table(celltypes,file = "E:/scrnaseqscript/foxa1_data/h2_data/H2luminalvsbasal.tsv", sep = '\t')
sparse_mat <- sparse_mat[,idx]
#rowidx <- sparse_mat[!apply(sparse_mat == 0, 1, all), , drop = FALSE]
rowidx <- which(rowSums(sparse_mat) > 0)
toexportpeaks <- peak_matrix@rowRanges[rowidx]
rtracklayer::export(con = 'E:/scrnaseqscript/foxa1_data/h2_data/H2_bedpeaksfoxa1_new.bed', object = toexportpeaks, ignore.strand = TRUE, format = 'bed')
sparse_mat <- sparse_mat[rowidx,]
write.table(colnames(sparse_mat), file = 'E:/scrnaseqscript/foxa1_data/h2_data/H2peakcelltypes.txt', sep = '\t')
writeMM(sparse_mat, 'E:/scrnaseqscript/foxa1_data/h2_data/foxA1_peak_counts_H2.mtx')
#only normalize if not already
#sparse_mat <- LogNormalize(sparse_mat)
#writeMM(sparse_mat,'foxa1_data/ev_data/normalizedpeaks.mtx')

#manage the rna matrix-looks like this one is already log normalized
#if the matrix is not log normalized, do that before subsetting.
rna_matrix <- rna_express@assays@data$GeneExpressionMatrix
idx_rna <- which((gsub("#.*","",colnames(rna_matrix)) == 'H2_multiome') == TRUE)

rna_matrix <- rna_matrix[,idx_rna]
rowidx <- which(rowSums(rna_matrix) > 0)
rna_matrix <- rna_matrix[rowidx, ]
length(colnames(rna_matrix))

write.table(rna_express@elementMetadata$name[rowidx], file = 'E:/scrnaseqscript/foxa1_data/h2_data/H2_rna_features.tsv', sep = '\t')
writeMM(rna_matrix, 'E:/scrnaseqscript/foxa1_data/h2_data/foxA1_rna_counts_H2.mtx')
write.table(colnames(rna_matrix), file = 'E:/scrnaseqscript/foxa1_data/h2_data/H2rnacelltypes.txt', sep = '\t')


