#install packages
devtools::install_github("mahmoudibrahim/genesorteR") 
remotes::install_github(repo = 'genecell/COSGR')

library(Seurat)
library(SeuratData)
library(genesorteR)
library(COSG)

data("pbmc3k.final")

#identify gene marker with COSG
starttime <- Sys.time()
COSG_markers <- cosg(
  pbmc3k.final,
  groups = 'all',
  assay = 'RNA',
  slot = 'data',
  mu = 1)
endtime <- Sys.time()
timecost <- endtime - starttime
print(timecost)



#identify gene marker with genesorteR
starttime <- Sys.time()
gs <- sortGenes(x = pbmc3k.final@assays$RNA@data, 
                classLabels = Idents(pbmc3k.final))
head(gs$specScore)
genesorteR_markers <- getMarkers(gs = gs, quant = 0.99)
pp <- plotMarkerHeat(gs$inputMat, 
                     as.vector(gs$inputClass), 
                     genesorteR_markers$markers, 
                     clusterGenes = T, 
                     outs = T)
pp$gene_class_info
endtime <- Sys.time()
timecost <- endtime - starttime
print(timecost)

library(patchwork)
library(ggplot2)
#for cosg
p1 <- DotPlot(pbmc3k.final, 
              assay = 'RNA',
              features = COSG_markers$names$B[1:10]) +
  ggtitle(label = 'B cell Markers with COSG') + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45))

p2 <- DotPlot(pbmc3k.final, 
              assay = 'RNA',
              features = which(pp$gene_class_info == 5) %>% names()) +
  ggtitle(label = 'B cell Markers with genesorteR') + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
p1+p2
