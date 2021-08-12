library(Seurat)
library(Matrix)
library(readr)
counts_peaks = readMM("foxa1_data/foxA1_peak_counts_G2.mtx")
cell_ids <- read.table("foxa1_data/G2peakcelltypes.txt", sep = "\t")
celltypes <- cell_ids$x
cell_ids <- cell_ids$x

colnames(counts_peaks) <- cell_ids

peaksz <- read.table("foxa1_data/bedpeaksfoxa1_new.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

colnames(peaksz) <- c("chr", "start", "end")

peaksz$region <- paste(peaksz$start, peaksz$end, sep="-")
peaksz$total <- paste(peaksz$chr, peaksz$region, sep="_")
rownames(counts_peaks) <- peaksz$total



library(devtools)
#devtools::install_github("aertslab/cisTopic")
library(cisTopic)

hold2 <- gsub("_", ":", rownames(counts_peaks))
rownames(counts_peaks) <- hold2

cisTopicObject <- createcisTopicObject(counts_peaks, project.name = 'cistopic')
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(50), seed=42, nCores=2, burnin = 120, iterations = 150, addModels=F)

saveRDS(cisTopicObject, file = "foxa1_data/cistopobj4")

#cisTopicObject <- readRDS('foxa1_data/cistopobj4')

cisTopicObject <- selectModel(cisTopicObject, type='maximum', select = 10)

cisTopicObject <- runUmap(cisTopicObject, target='cell', seed = 42, min_dist = 0.3)
plot(cisTopicObject@dr$cell$Umap, cex=0.1)

par(mfrow=c(2,5))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr = 'Probability', colorBy = NULL, cex.legend = 0.8, factor.max = .75, dim = 2, legend = TRUE, col.low = 'darkgreen', col.mid = 'yellow', col.high = 'brown1', intervals=20 )
cellTopicHeatmap(cisTopicObject, method='Probability')

cis2 <- modelMatSelection(cisTopicObject, target='cell',method='Probability')

write.table(cis2, file = "foxa1_data/cisTopic10Topicsoutput.tsv", sep= "\t")

umapcoords <- cisTopicObject@dr$cell$Umap

write.table(umapcoords, file = "foxa1_data/cisTopic10TopicsUmap.tsv", sep= "\t")
