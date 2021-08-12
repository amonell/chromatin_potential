library(ChIPpeakAnno)
library(chromVAR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
mm10genome <- load('mm10_refSeq.RData')
names(mm10TSSRanges) <- mm10TSSRanges$gene_name
mart = useMart(biomart = "ensembl", dataset= "mmusculus_gene_ensembl")

peaks2 <- getPeaks("foxa1_data/ev_data/EV_bedpeaksfoxa1_new.bed", sort_peaks = TRUE)
annotations = annotatePeakInBatch(peaks2, AnnotationData = mm10TSSRanges, output = c("both"), multiple = TRUE, maxgap = 50000L, PeakLocForDistance = c("middle"), FeatureLocForDistance = c("middle"))
hold <- strsplit(names(annotations), "\\.")
hold2 <- sapply(hold, "[", 2)
names(annotations) = hold2

#gene <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'ensembl_gene_id', values = names(annotations), mart = mart)

library(dplyr)
library(rtracklayer)
library(Repitools)


annotations_2 <- annotations[(abs(elementMetadata(annotations)[,7]) <50000)]
length(annotations_2)
rtracklayer::export.bed(annotations_2, "foxa1_data/ev_data/EV_foxa1_annotations.bed")

