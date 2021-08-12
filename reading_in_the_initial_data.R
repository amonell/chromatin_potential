rna_express <- readRDS('foxa1_data/foxa1_expressionmatrix_SE.rds')
peak_matrix <- readRDS('foxa1_data/foxa1_peakmatrix_RSE.rds')
umap_coords <- read.table('foxa1_data/foxa1_archR_Umap.txt', sep = ' ')


library(rtracklayer)

#rtracklayer::export(con = 'foxa1_data/bedpeaksfoxa1_new.bed', object = peak_matrix@rowRanges, ignore.strand = TRUE, format = 'bed')
sparse_mat <- peak_matrix@assays@data$PeakMatrix

library(Matrix)

library(Seurat)


#manage the peak matrix
idx <- which((gsub("#.*","",colnames(sparse_mat)) == 'EV_multiome') == TRUE)
celltypes <- rna_express$luminal_basal[idx]
write.table(celltypes,file = "foxa1_data/ev_data/EVluminalvsbasal.tsv", sep = '\t')
sparse_mat <- sparse_mat[,idx]
#rowidx <- sparse_mat[!apply(sparse_mat == 0, 1, all), , drop = FALSE]
rowidx <- which(rowSums(sparse_mat) > 0)
toexportpeaks <- peak_matrix@rowRanges[rowidx]
rtracklayer::export(con = 'foxa1_data/ev_data/EV_bedpeaksfoxa1_new.bed', object = toexportpeaks, ignore.strand = TRUE, format = 'bed')
sparse_mat <- sparse_mat[rowidx,]
write.table(colnames(sparse_mat), file = 'foxa1_data/ev_data/EVpeakcelltypes.txt', sep = '\t')
writeMM(sparse_mat, 'foxa1_data/ev_data/foxA1_peak_counts_G2.mtx')
sparse_mat <- LogNormalize(sparse_mat)
writeMM(sparse_mat,'foxa1_data/ev_data/normalizedpeaks.mtx')

#manage the rna matrix-looks like this one is already log normalized
#if the matrix is not log normalized, do that before subsetting.
rna_matrix <- rna_express@assays@data$GeneExpressionMatrix
idx_rna <- which((gsub("#.*","",colnames(rna_matrix)) == 'EV_multiome') == TRUE)
rna_matrix <- rna_matrix[,idx_rna]
length(colnames(rna_matrix))

write.table(rna_express@elementMetadata$name, file = 'foxa1_data/ev_data/EV_rna_features.tsv', sep = '\t')
writeMM(rna_matrix, 'foxa1_data/ev_data/foxA1_rna_counts_EV.mtx')
write.table(colnames(rna_matrix), file = 'foxa1_data/ev_data/EVrnacelltypes.txt', sep = '\t')


