# chromatin potential -----------------------------------------------------
SEsum_smooth <- as.matrix(read.table('E:/scrnaseqscript/foxa1_data/m2_data/smoothed_dorc_table_M2_700.csv', sep = ',', header = TRUE))
forcols <- read.table('E:/scrnaseqscript/foxa1_data/m2_data/smoothed_dorc_columns_M2_700.csv', sep = ',')
temp_cell_index <- read.table('E:/scrnaseqscript/foxa1_data/m2_data/M2peakcelltypes.txt', sep = "\t", header = TRUE)
rownames_SEsum_smooth = SEsum_smooth[,1]
cellindex <- temp_cell_index$x
colnames <- cellindex
SEsum_smooth <- SEsum_smooth[,-c(1)]


RNA_smooth <- as.matrix(read.table('E:/scrnaseqscript/foxa1_data/m2_data/smoothed_rna_table_M2_700.csv', sep = ',', header = TRUE))
rownames <- RNA_smooth[,'gene']
RNA_smooth <- RNA_smooth[,-c(1)]

SEsum_smooth <- matrix(as.numeric(SEsum_smooth), ncol = ncol(SEsum_smooth))
RNA_smooth <- matrix(as.numeric(RNA_smooth), ncol = ncol(RNA_smooth))
rownames(RNA_smooth) <- rownames
rownames(SEsum_smooth) <- rownames_SEsum_smooth
colnames(RNA_smooth) <- colnames
colnames(SEsum_smooth) <- colnames

celltype <- read.table('E:/scrnaseqscript/foxa1_data/m2_data/M2luminalvsbasal.tsv', sep = '\t')
celltype <- celltype$x

library(data.table)

## find high var genes
x <- log10(rowMeans(SEsum_smooth[, cellindex]))
std <- apply(SEsum_smooth[, cellindex], 1, sd)
y = std^2/rowMeans(SEsum_smooth[, cellindex]) 
plot(x, y)
#genes <- names(y[y > 0.04675])
length(y)
genes <- names(y[y > 0.3]) 
length(genes) # try to get top half of genes

library(readr)
mapcoords <- read.table('E:/scrnaseqscript/foxa1_data/m2_data/cisTopic40M2_Umap', sep = '\t', header = TRUE)
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
newarrows <- read.table('E:/scrnaseqscript/foxa1_data/m2_data/cisTopic40M2.tsv', sep = '\t')

newarrowscolumns <- colnames(newarrows)
#colnames(newarrows) <- newarrowscolumns
newarrowscolumns <- gsub("e.", "e#", newarrowscolumns, fixed = TRUE)
colnames(newarrows) <- newarrowscolumns <- gsub(".", "-", newarrowscolumns, fixed = TRUE)
#newarrowscolumns <- subset(newarrowscolumns, newarrowscolumns %in% cellindex)
newarrows <- newarrows[,cellindex]
newends <- newarrows
for (i in 1:nrow(ATAC.RNA.KNN3$test_knn_idx)){
  startpoint <- newarrows[,i]
  endpoint <- rowMeans(newarrows[,ATAC.RNA.KNN3$test_knn_idx[i,]])
  newends[i] <- endpoint
}

colnames(newends) <- paste("end", colnames(newends), sep = '_')


#remotes::install_github("jlmelville/uwot")
library(uwot)
library(RcppAnnoy)
neural_embedding <- uwot::umap(t(newarrows), ret_model = TRUE)
#names(neural_embedding) <- c("layout", "data", "nn_index", "config")
ending_umaps <- umap_transform(t(newends), neural_embedding)
starting_umaps <- umap_transform(t(newarrows), neural_embedding)
#-------------------------------------------------------------------------------------------------------
newarrows <- cbind(newarrows, newends)


#if you want to smooth in high dimension, smooth this umap.  THis will take a long time however.  Cannot do with UMAP neural network
#library(FNN)
#k = 15 
#knn.norm = get.knn(t(newarrows), k = k)
#knn.norm <- knn.norm$nn.index

#umap.smooth <- newarrows
#print(nrow(umap.smooth))
#for (i in 1:nrow(umap.smooth)){
#  print(i)
#  for (j in 1:ncol(umap.smooth)){
    
#    umap.smooth[i,j] <- mean(newarrows[j,knn.norm[i, ]])
    #umap.smooth[i,j] <- mean(umap.raw[knn.norm[i, ],j])
#  }
#}
####################
library(cisTopic)
set.seed(42)
library(scater)
#dont run this line if high dimension smoothing
umap.smooth <- newarrows
umap.smooth <- calculateUMAP(umap.smooth)
colnames(umap.smooth) <- c('UMAP1', 'UMAP2')
rownames(umap.smooth) <- colnames(newarrows)

tomerge <- umap.smooth[colnames(newends),]
umap.smooth <- umap.smooth[newarrowscolumns, ]
colnames(tomerge) <- c('UMAP1.1', 'UMAP2.1')
umap.smooth <- cbind(umap.smooth, tomerge)
#-----------------------------------------------------------------------------------------------------
umap.smooth <- cbind(starting_umaps, ending_umaps)
colnames(umap.smooth) <- c('UMAP1', 'UMAP2', 'UMAP1.1', 'UMAP2.1')

temp <- as.data.frame(umap.smooth)

temp$celltype <- celltype

# remove ORS cells for visulization
#temp2 <- intersect(rownames(temp), colnames(atac.se[ ,atac.se$L1 | atac.se$L2 | atac.se$L3]))
#temp2 <- intersect(rownames(temp), cellindex)
#temp <- temp[temp2, ]
#temp <- subset(temp, temp$umap1 < -3 & temp$umap2 < 0)
# filter out top 5% long arrow
temp$arrow.legnth <- ((temp$UMAP1.1- temp$UMAP1)^2 + (temp$UMAP2.1- temp$UMAP2)^2)^0.5
hist(temp$arrow.legnth, breaks = 100)
temp$umap1.1[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)] <- temp$umap1[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)]
temp$umap2.1[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)] <- temp$umap2[temp$arrow.legnth > quantile(temp$arrow.legnth, probs=0.95)]



## smooth in umap space
umap.raw <- temp; head(umap.raw)
pcs <- data.frame(umap1=temp$UMAP1, umap2=temp$UMAP2)
library(KernelKnn)
k = 15 
knn.norm = get.knn(pcs, k = k)
knn.norm <- knn.norm$nn.index

umap.smooth <- umap.raw
for (i in 1:nrow(umap.smooth)){
  umap.smooth[i,3] <- mean(umap.raw[knn.norm[i, ],3])
  umap.smooth[i,4] <- mean(umap.raw[knn.norm[i, ],4])
}


umap.smooth$arrow.legnth <- ((umap.smooth$UMAP1.1- umap.smooth$UMAP1)^2 + (umap.smooth$UMAP2.1- umap.smooth$UMAP2)^2)^0.5

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

umap.smooth <- umap.smooth %>% mutate(GroupCol =
                                        case_when(celltype == 'Basal' ~ 2, 
                                                  celltype == 'Luminal' ~ 0,
                                                  celltype == '' ~ 1)
)

randocells <- subset(umap.smooth, umap.smooth$celltype == '')



par(mfrow=c(1,1))

plot(umap.smooth$UMAP1, umap.smooth$UMAP2, col=umap.smooth$Group, xlim=c(-8, 9) , ylim=c(-10, 5)) 
#legend(7,4.3,unique(umap.smooth$umap1),col=1:length(umap.smooth$umap2),pch=1)

#plot(randocells$umap1, randocells$umap2, col = umap.smooth$Group)

ggplot(umap.smooth, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(color=arrow.legnth), size=0.3, stroke=0, alpha = 1) +
  
  geom_arrow(aes(dx = (UMAP1.1-UMAP1),dy = (UMAP2.1-UMAP2)), skip = 2,
             color = "black",lineend="round", size=0.2, arrow.length = 0.6) +
  scale_mag(max_size = 2, guide = "none") +
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))) +
  pretty_plot(fontsize = 5) + L_border() + xlim(-8, 9) + ylim(-10, 5)
guides(colour = FALSE)

ggplot(umap.smooth, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(color=GroupCol), size=2, stroke=0, alpha = 0.1) +
  
  geom_arrow(aes(dx = (UMAP1.1-UMAP1),dy = (UMAP2.1-UMAP2)), skip = 1.5,
             color = "black",lineend="round", size=0.2, arrow.length = 0.6) +
  scale_mag(max_size = 2, guide = "none") +
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))) +
  pretty_plot(fontsize = 5) + L_border() + xlim(-8, 9) + ylim(-10, 5)

ggsave("my.chromatin.ptime.png", width=3.8, height=3, units = "cm", dpi=1200)


