# rm(list=ls())
# options(stringsAsFactors = F)
# library(Seurat)
# library(ggplot2)
# library(clustree)
# library(cowplot)
# library(dplyr)
# getwd()
# setwd('3-cell/')
# sce.all=readRDS( "../2-harmony/sce.all_int.rds")
# sce.all

sce.all
sce=sce.all 
library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',#T细胞#
                   'CD19', 'CD79A', 'MS4A1' ,#B细胞#
                   'IGHG1', 'MZB1', 'SDC1',#浆细胞#
                   'CD68', 'CD163', 'CD14', #巨噬细胞#
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB',  # mac
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'FCGR3A',
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'KLRB1','NCR1', # NK 
                   'FGF7','MME', 'ACTA2', ## fibo 
                   'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo 
                   'MKI67' , 'TOP2A', 
                   'PECAM1', 'VWF',  ## endo 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)#函数 str_to_upper 和 str_to_lower 是大小写转换#
genes_to_check
p <- DotPlot(sce.all, features = unique(genes_to_check),
             assay='RNA'  )  + coord_flip()

p 
ggsave('check_last_markers.pdf',height = 11,width = 11)#什么是last#

DimPlot(sce.all, reduction = "umap", 
        group.by = "RNA_snn_res.0.01",label = T) 
ggsave('umap_by_RNA_snn_res.0.01.pdf')

cg=c('NUPR1','RARRES2','MLEC','TMSB10',
     'HES1','SNHG5','TFF3','WFDC2','MUC2','ITLN1','SPINK4',
     'MGST1','ADH1C','UGT2B17','CRYBA2','SCGN','PCSK1N','PTMS',
     'TUBA1A','MT2A','MT1G','KRT19','PHGR1','S100A6','SELENBP1','CA1',
     'CEACAM1','CEACAM7','GUCA2A','SLC26A3','AQP8','CA7','BEST4',
     'OTOP2','MT1H','LYPD8','CA4','HLA−DRA','HLA−DQB1','HLA−DPB1',
     'CD74','CD83','EEF1A1','CREM','CCL5','TMSB4X',
     'B2M','SRGN','LAPTM5','VIM','CD44','TPSAB1','TPSB2')
p <- DotPlot(sce.all, features = unique(cg),
             assay='RNA'  )  + coord_flip()


p 
ggsave('check_paper_markers.pdf',height = 11,width = 11)

library(ggplot2) 
genes_to_check = c(  'PROM1', 'CD44' , 'THY1','ACTA2','MKI67', 'ESR1',
                     'EPCAM', 'KRT19', 
                     'KRT8',  'KRT18', 'KRT5',  'KRT14', 
                     'ALDH1A1', 'CD24' )
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip()

p  
ggsave(filename = 'DotPlot_epi_recluster.pdf')


# 需要自行看图，定细胞亚群：

celltype=data.frame(ClusterID=0:12,
                    celltype= 0:12) 
#定义细胞亚群
celltype[celltype$ClusterID %in% c( 0,3,4,5,7,10 ),2]='Mac' 
celltype[celltype$ClusterID %in% c( 1,2),2]='Tcells' 
celltype[celltype$ClusterID %in% c( 6),2]='PKIB'
celltype[celltype$ClusterID %in% c( 9 ),2]='Bcells'  
celltype[celltype$ClusterID %in% c( 8 ),2]='cycling'  
celltype[celltype$ClusterID %in% c( 12 ),2]='IDO1-DC'  
celltype[celltype$ClusterID %in% c( 11  ),2]='KRT19-epi'   

head(celltype)
celltype
table(celltype$celltype)
sce.all@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)


th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
 
p <- DotPlot(sce.all, features = unique(genes_to_check),
             assay='RNA' ,group.by = 'celltype' )  + coord_flip()  +th

p
ggsave(plot=p, filename="check_marker_by_celltype.pdf",#不懂为啥要checkmarker#
       width = 7 ,height = 8)

DotPlot(sce.all, features = genes_to_check,
        assay='RNA' ,group.by = 'celltype' ) + 
  coord_flip()+ scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
  NULL

table(sce.all@meta.data$celltype,sce.all@meta.data$RNA_snn_res.0.8)
#  

library(patchwork)#接下来这一串应该就是三种算法分群的比较？#
p_all_markers=DotPlot(sce.all, features = genes_to_check,
                      assay='RNA' ,group.by = 'celltype' )  + coord_flip()+th
p_umap=DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T)
p_all_markers+p_umap
ggsave('markers_umap_by_celltype.pdf',width = 12,height = 8)
p_harmony=DimPlot(sce.all, reduction = "harmony", group.by = "celltype",label = T)
p_all_markers+p_harmony
ggsave('markers_harmony_by_celltype.pdf',width = 12,height = 8)

sce.all=RunTSNE(sce.all,  dims = 1:15, 
        reduction = "harmony")
p_tsne=DimPlot(sce.all, reduction = "tsne", group.by = "celltype",label = T)
p_all_markers+p_tsne
ggsave('markers_tsne_by_celltype.pdf',width = 12,height = 8)

phe=sce.all@meta.data
save(phe,file = 'phe-by-markers.Rdata')


sce.all
table(Idents(sce.all))  
Idents(sce.all)=sce.all$celltype
table(Idents(sce.all))  

library(future)
# check the current active plan
plan()
plan("multiprocess", workers = 8)
plan()

sce.markers <- FindAllMarkers(object = sce.all, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers)
pro='markers'
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce.all,top10$gene,size=3)
ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'),height = 15)

library(dplyr) 
top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
DoHeatmap(sce.all,top3$gene,size=3)
ggsave(paste0(pro,'DoHeatmap_check_top3_markers_by_clusters.pdf'))
p <- DotPlot(sce.all, features = unique(top3$gene),
             assay='RNA'  )  + coord_flip()+th

p
ggsave(paste0(pro,'DotPlot_check_top3_markers_by_clusters.pdf'))
save(sce.markers,file = paste0(pro, 'sce.markers.Rdata'))


# 
#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1
library(ggplot2) 
ggsave(filename="Vlnplot1.pdf",plot=p1)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all,  features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2	
ggsave(filename="Vlnplot2.pdf",plot=p2)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", 
                  pt.size = 0.5)
ggsave(filename="Scatterplot.pdf",plot=p3)
