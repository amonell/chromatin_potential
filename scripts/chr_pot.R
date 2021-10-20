# chromatin potential -----------------------------------------------------
SEsum_smooth <- as.matrix(read.table('foxa1_data/ev_data/smoothed_dorc_table_EV_40.csv', sep = ',', header = TRUE))
forcols <- read.table('foxa1_data/ev_data/smoothed_dorc_columns_EV.csv', sep = ',')
temp_cell_index <- read.table('foxa1_data/ev_data/EVpeakcelltypes.txt', sep = "\t", header = TRUE)
rownames_SEsum_smooth = SEsum_smooth[,1]
cellindex <- temp_cell_index$x
colnames <- cellindex
SEsum_smooth <- SEsum_smooth[,-c(1)]


RNA_smooth <- as.matrix(read.table('foxa1_data/ev_data/smoothed_rna_table_EV_40.csv', sep = ',', header = TRUE))
rownames <- RNA_smooth[,'gene']
RNA_smooth <- RNA_smooth[,-c(1)]

SEsum_smooth <- matrix(as.numeric(SEsum_smooth), ncol = ncol(SEsum_smooth))
RNA_smooth <- matrix(as.numeric(RNA_smooth), ncol = ncol(RNA_smooth))
rownames(RNA_smooth) <- rownames
rownames(SEsum_smooth) <- rownames_SEsum_smooth
colnames(RNA_smooth) <- colnames
colnames(SEsum_smooth) <- colnames

celltype <- read.table('foxa1_data/ev_data/EVluminalvsbasal.tsv', sep = '\t')
celltype <- celltype$x

library(data.table)

## find high var genes
x <- log10(rowMeans(SEsum_smooth[, cellindex]))
std <- apply(SEsum_smooth[, cellindex], 1, sd)
y = std^2/rowMeans(SEsum_smooth[, cellindex]) 
plot(x, y)
#genes <- names(y[y > 0.04675]) 
length(y)
genes <- names(y[y > 0.5]) 
length(genes) # try to get top half of genes

library(readr)
mapcoords <- read.table('foxa1_data/ev_data/cisTopic40EV_Umap', sep = '\t', header = TRUE)
#uncomment if using archr umap coords instead of cisTopic
##########################################################
#rownames(mapcoords) <- cellindex
#mapcoords <- subset(mapcoords, select = -c(X))
#########################################################
mapcoords <- subset(mapcoords, rownames(mapcoords) %in% cellindex)
mapcoords$celltype <- celltype

## knn.index.dist
library(KernelKnn)
# this step is very slow
ATAC.RNA.KNN3 <- knn.index.dist(t(SEsum_smooth[genes,cellindex]),t(RNA_smooth[genes,cellindex]), method="pearson_correlation", k=10, threads = 3) # pearson_correlation
dim(ATAC.RNA.KNN3$test_knn_idx)
# saveRDS(ATAC.RNA.KNN3, "ATAC.RNA.KNN3.rds")
# ATAC.RNA.KNN3 <- readRDS('ATAC.RNA.KNN3.rds')

# smooth in ATAC umap with 10 KNNs
umap.old <-  data.frame(umap1=mapcoords$UMAP1, umap2=mapcoords$UMAP2)
umap.new <- umap.old
umap.new$mean.dist <- 0
for (i in 1:nrow(umap.new)){
  umap.new[i, 1:2] <- colMeans(umap.old[ATAC.RNA.KNN3$test_knn_idx[i, ], ])
  umap.new$mean.dist[i] <- max(dist(umap.old[ATAC.RNA.KNN3$test_knn_idx[i, ], ],method = "euclidean"))
}


plot(umap.old[, 1:2])
plot(umap.new[, 1:2])

# remove ORS cells for visulization
temp <- data.frame(umap.old, umap.new)
head(temp)
rownames(temp) <- rownames(mapcoords)
#temp2 <- intersect(rownames(temp), cellindex)
#temp <- temp[temp2, ]
#temp <- subset(temp, temp$umap1 < -3 & temp$umap2 < 0)

# filter out top 5% long arrow
temp$arrow.legnth <- ((temp$umap1.1- temp$umap1)^2 + (temp$umap2.1- temp$umap2)^2)^0.5
hist(temp$arrow.legnth, breaks = 100)
temp$umap1.1[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)] <- temp$umap1[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)]
temp$umap2.1[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)] <- temp$umap2[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)]


library(FNN)
## smooth in umap space
umap.raw <- temp; head(umap.raw)
pcs <- data.frame(umap1=temp$umap1, umap2=temp$umap2)
k = 15 
knn.norm = get.knn(pcs, k = k)
knn.norm <- knn.norm$nn.index

umap.smooth <- umap.raw
for (i in 1:nrow(umap.smooth)){
  umap.smooth[i,3] <- mean(umap.raw[knn.norm[i, ],3])
  umap.smooth[i,4] <- mean(umap.raw[knn.norm[i, ],4])
}


umap.smooth$arrow.legnth <- ((umap.smooth$umap1.1- umap.smooth$umap1)^2 + (umap.smooth$umap2.1- umap.smooth$umap2)^2)^0.5

chromatin.potential <- umap.smooth 
# saveRDS(chromatin.potential, "chromatin.potential.rds")

scale.factor = 0.1
library(ggplot2)
library(FNN)
library(BuenColors)
library(RColorBrewer)
library(metR)

umap.smooth$celltype <- celltype
library(dplyr)
umap.smooth <- umap.smooth %>% mutate(Group =
                     case_when(celltype == 'Basal' ~ "Red", 
                               celltype == 'Luminal' ~ "Blue",
                               celltype == '' ~ "Gray")
)

randocells <- subset(umap.smooth, umap.smooth$celltype == '')

par(mfrow=c(1,1))

plot(umap.smooth$umap1, umap.smooth$umap2, col=umap.smooth$Group, xlim = c(-10, 10), ylim = c(-10, 10))# + legend(7,4.3,unique(umap.smooth$Group),col=1:length(umap.smooth$Group),pch=1) 

#plot(randocells$umap1, randocells$umap2, col = umap.smooth$Group)

ggplot(umap.smooth, aes(x=umap1, y=umap2)) +
  geom_point(aes(color=arrow.legnth), size=0.3, stroke=0, alpha = 1) +
  
  geom_arrow(aes(dx = (umap1.1-umap1),dy = (umap2.1-umap2)), skip = 2,
             color = "black",lineend="round", size=0.2, arrow.length = 0.5) +
  scale_mag(max_size = 3, guide = "none") +
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))) +
  pretty_plot(fontsize = 5) + L_border() + xlim(-10, 10) + ylim(-10, 10)
  guides(colour = FALSE)
ggsave("my.chromatin.ptime.png", width=3.8, height=3, units = "cm", dpi=1200)



