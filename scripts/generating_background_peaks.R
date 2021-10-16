
peakssp <- read.table("E:/scrnaseqscript/foxa1_data/h2_data/H2_peaks_used_in_background.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
peaksz <- read.table("E:/scrnaseqscript/foxa1_data/h2_data/H2_bedpeaksfoxa1_new.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

colnames(peaksz) <- c("chr", "start", "end")

peaksz$region <- paste(peaksz$start, peaksz$end, sep="_")
peaksz$total <- paste(peaksz$chr, peaksz$region, sep="_")


library(Seurat)
library(Matrix)
library(readr)

counts = readMM("E:/scrnaseqscript/foxa1_data/h2_data/foxA1_peak_counts_H2.mtx")
#cell_ids should be the names of the atac cell barcodes
cell_ids <- read.table("E:/scrnaseqscript/foxa1_data/h2_data/H2peakcelltypes.txt", sep = '\t')
#celltypes <- cell_ids$celltype
cell_ids <- cell_ids$x

colnames(counts) <- cell_ids
rownames(counts) <- peaksz$total

idx <- which(rownames(counts) %in% peakssp$peak)
counts <- counts[idx,]
peaks <- data.frame('chr' = peaksz$chr[idx], 'start' = peaksz$start[idx], 'end' = peaksz$end[idx])

peaks$region <- paste(peaks$start, peaks$end, sep="-")
peaks$total <- paste(peaks$chr, peaks$region, sep="_")

library(SummarizedExperiment)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
length(rownames(peaks))

peaks <- makeGRangesFromDataFrame(peaks)
fragment_counts <- SummarizedExperiment(assays = list(counts = counts), rowRanges = peaks)#, colRanges = cell_ids)
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Mmusculus.UCSC.mm10)

fragment_counts <- getBackgroundPeaks(fragment_counts, bias = rowData(fragment_counts)$bias, niterations = 100, w = 0.5, bs = 100)

#row names are the second column of the peak names from peakssp
write.table(fragment_counts, file = "E:/scrnaseqscript/foxa1_data/h2_data/backgroundpeaks_H2.txt", row.names = peakssp$peak)
