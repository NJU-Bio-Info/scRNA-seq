library(Seurat)
library(SeuratData)
library(DT)

#加载内置数据集
data("pbmc3k")
pbmc3k <- pbmc3k.final

invisible(gc())

#查看细胞类群
celltype <- table(pbmc3k@meta.data$seurat_annotations)
datatable(as.data.frame(celltype))

#经典的FindAllMarkers()函数
starttime <- Sys.time()
Seurat_markers <- FindAllMarkers(pbmc3k, only.pos = TRUE)
endtime <- Sys.time()
timecost <- endtime - starttime
print(timecost)


library(COSG)
library(dplyr)
Seurat_markers %>% count(cluster)
starttime <- Sys.time()
COSG_markers <- cosg(
  pbmc3k,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=900)
endtime <- Sys.time()
timecost <- endtime - starttime
print(timecost)

library(patchwork)
library(ggplot2)
p1 <- DotPlot(pbmc3k, 
              assay = 'RNA',
              features = subset(Seurat_markers, cluster == 'B')$gene[1:10]) + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45)) + 
  NoLegend()
p2 <- DotPlot(pbmc3k, 
              assay = 'RNA',
              features = COSG_markers$names$B[1:10]) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
p1 + p2
