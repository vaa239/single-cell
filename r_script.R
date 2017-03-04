setwd("/home/yaroslav/Documents/biohack/matrices_mex/hg19")
library(data.table)
library(Seurat)
library(dplyr)
library(Matrix)
source("./functions.R")
# matrix <- fread("matrix.mtx")
# barcodes <- fread("barcodes.tsv")
# genes <- fread("genes.tsv")
# load("jurkat.Robj")

jurkat.data <- Read10X("/home/yaroslav/Documents/biohack/matrices_mex/hg19")
jurkat <- new("seurat", raw.data = jurkat.data)
jurkat <- Setup(jurkat, min.cells = 3, min.genes = 200, do.logNormalize = T, 
              total.expr = 1e4, project = "10X_PBMC")

mito.genes <- grep("^MT-", rownames(jurkat@data), value = T)
percent.mito <- colSums(expm1(jurkat@data[mito.genes, ]))/colSums(expm1(jurkat@data))

jurkat <- AddMetaData(jurkat, percent.mito, "percent.mito")
# VlnPlot(jurkat, c("nGene", "nUMI", "percent.mito"), nCol = 3)

# par(mfrow = c(1, 2))
# GenePlot(jurkat, "nUMI", "percent.mito")
# GenePlot(jurkat, "nUMI", "nGene")

jurkat <- SubsetData(jurkat, subset.name = "nGene", accept.high = 3200)
jurkat <- SubsetData(jurkat, subset.name = "percent.mito", accept.high = 0.05)

jurkat <- RegressOut(jurkat, latent.vars = c("nUMI", "percent.mito"))

jurkat <- MeanVarPlot(jurkat ,fxn.x = expMean, fxn.y = logVarDivMean, 
                    x.low.cutoff = 0.0125, x.high.cutoff = 3, 
                    y.cutoff = 0.5, do.contour = F)

# MeanVarPlot(jurkat ,fxn.x = expMean, fxn.y = logVarDivMean, 
#             x.low.cutoff = 0.0125, x.high.cutoff = 3, 
#             y.cutoff = 0.5, do.contour = F)

jurkat <- PCA(jurkat, pc.genes = jurkat@var.genes, do.print = TRUE, 
            pcs.print = 5, genes.print = 5)

jurkat <- ProjectPCA(jurkat)

PrintPCA(jurkat, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
VizPCA(jurkat, 1:2)
PCAPlot(jurkat, 1, 2)
PCHeatmap(jurkat, pc.use = 1, cells.use = 100, do.balanced = TRUE)
PCHeatmap(jurkat, pc.use = 1:12, cells.use = 500, do.balanced = TRUE,
#           label.columns = FALSE, use.full = FALSE)

jurkat <- JackStraw(jurkat, num.replicate = 100, do.print = FALSE)

JackStrawPlot(jurkat, PCs = 1:12)

PCElbowPlot(jurkat)

jurkat <- FindClusters(jurkat, pc.use = 1:10 , resolution = 0.6, 
                       print.output = 0, save.SNN = T)

jurkat <- RunTSNE(jurkat, dims.use = 1:10, do.fast = T)

TSNEPlot(jurkat)
save(jurkat, 
     file = "/home/yaroslav/Documents/biohack/matrices_mex/hg19/jurkat.Robj")

jurkat.markers <- FindAllMarkers(jurkat, only.pos = TRUE, 
                               min.pct = 0.25, thresh.use = 0.25)
jurkat.markers %>% group_by(cluster) %>% top_n(2, avg_diff)
current.cluster.ids <- c(0, 1, 2, 3, 4)
new.cluster.ids <- c(0, 1, 2, 3, 4)
jurkat@ident <- plyr::mapvalues(jurkat@ident, from = current.cluster.ids, 
                              to = new.cluster.ids)
TSNEPlot(jurkat, do.label = T, pt.size = 0.5)


tsne <- jurkat@tsne.rot
tsne$cluster_id <- jurkat@ident
tsne$id <- paste(jurkat@ident)
write.tsv(tsne, dir=".")


