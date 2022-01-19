dir.create("2-harmony")
getwd()
setwd("2-harmony")

# sce.all=readRDS("../1-QC/sce.all_qc.rds")
sce=sce.all.filt 
sce#过滤过的#
sce <- NormalizeData(sce, #NormalizeData()这个函数是首先对基因的reads数进行了同一文库大小的校正#
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4) 
sce <- FindVariableFeatures(sce)#鉴定在细胞间表达高度变化的基因，后续研究需要集中于这部分基因。#
sce <- ScaleData(sce)#做标准化#
sce <- RunPCA(sce, features = VariableFeatures(object = sce))#主成分分析，抓出主要矛盾，找出对细胞异质性影响比较大的基因#

library(harmony)
seuratObj <- RunHarmony(sce, "orig.ident")#Harmony是在去除批次效应#
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                     reduction = "harmony")#跑一次UMAP可以整合十五个主成分的信息，PCA一次最多只能展示两个#
DimPlot(seuratObj,reduction = "umap",label=T ) #主成分的展示#

sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",
                     dims = 1:15)#Which dimensions to use as input features, used only if features is NULL#
#首先基于pca维度中（先前计算的pca数据）计算欧式距离（the euclidean distance），然后根据两个细胞在局部的重合情况（Jaccard 相似系数）优化两个细胞之间的边缘权值。#
sce.all=sce
#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all=FindClusters(sce.all, #graph.name = "CCA_snn", 
                       resolution = res, algorithm = 1)
}
colnames(sce.all@meta.data)
apply(sce.all@meta.data[,grep("RNA_snn",colnames(sce.all@meta.data))],2,table)#就是把table这个函数用到前面取出来的那个矩阵上#
p1_dim=plot_grid(ncol = 3, DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 14)#观察了分辨率为0.01，0.1，0.2的低分辨率分群效果#

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 18)#观察了高分辨率分群效果#


p2_tree=clustree(sce.all@meta.data, prefix = "RNA_snn_res.")#根据这个也看不出来为啥要0.8讷#
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf")


#接下来分析，按照分辨率为0.8进行 
sel.clust = "RNA_snn_res.0.8"
sce.all <- SetIdent(sce.all, value = sel.clust)#Value:要从对象元数据或身份本身中提取的身份名称#
table(sce.all@active.ident) 
saveRDS(sce.all, "sce.all_int.rds")#保存了这个对象#


###### step5:检查常见分群情况  ######

setwd('../')
dir.create("3-cell")
setwd("3-cell")  

DimPlot(sce.all, reduction = "umap", group.by = "seurat_clusters",label = T) 
DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.8",label = T) 
ggsave('umap_by_RNA_snn_res.0.8.pdf')


# Tumor cells were identified using MLANA, MITF, and DCT.
# Tumor cells were further divided into subgroups 
# by expression of PRAME and GEP genes

library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E',   # T cells
                   'KLRD1', # natural killer (NK) cells
                   'FCGR3B', #neutrophils
                   'MS4A1' ,'IGHG4', # b cells 
                   'RPE65','RCVRN',
                   'CD68', 'CD163', 'CD14',   'MKI67' ,'TOP2A',
                   'MLANA', 'MITF',  'DCT','PRAME' , 'GEP' )
library(stringr)  
p_paper_markers <- DotPlot(sce.all, features = genes_to_check,
                           assay='RNA'  )  + coord_flip()

p_paper_markers#这些markers为啥图里的跟我们检查的不对应呀，是因为有些没检测出来吗#
ggsave(plot=p_paper_markers,
       filename="check_paper_marker_by_seurat_cluster.pdf",width = 12)



# T Cells (CD3D, CD3E, CD8A), 
# B cells (CD19, CD79A, MS4A1 [CD20]), 
# Plasma cells (IGHG1, MZB1, SDC1, CD79A), 
# Monocytes and macrophages (CD68, CD163, CD14),
# NK Cells (FGFBP2, FCG3RA, CX3CR1),  
# Photoreceptor cells (RCVRN), 
# Fibroblasts (FGF7, MME), 
# Endothelial cells (PECAM1, VWF). 
# epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
#   immune (CD45+,PTPRC), epithelial/cancer (EpCAM+,EPCAM), 
# stromal (CD10+,MME,fibo or CD31+,PECAM1,endo) 

library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',#CD4是哪里来的##Tcells组#
                   'CD19', 'CD79A', 'MS4A1' ,#Bcells组#
                   'IGHG1', 'MZB1', 'SDC1',#Plasma cells浆细胞组#
                   'CD68', 'CD163', 'CD14', #巨噬细胞组#
                   'TPSAB1' , 'TPSB2',  # mast cells,#肥大细胞组#
                   'RCVRN','FPR1' , 'ITGAM' ,#Photoreceptor cells感光细胞组 后面俩是哪来的#
                   'FGF7','MME', 'ACTA2',#成纤维细胞组 后面1是哪来的#
                   'PECAM1', 'VWF', #血管内皮细胞组#
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )#肿瘤细胞组#
library(stringr)  
p_all_markers <- DotPlot(sce.all, features = genes_to_check,
                         assay='RNA'  )  + coord_flip()

p_all_markers
ggsave(plot=p_all_markers, filename="check_all_marker_by_seurat_cluster.pdf")



genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                   'CCR7', 'SELL' , 'TCF7','CXCR6' , 'ITGA1',
                   'FOXP3', 'IL2RA',  'CTLA4','GZMB', 'GZMK','CCL5',
                   'IFNG', 'CCL4', 'CCL3' ,
                   'PRF1' , 'NKG7') #T细胞组#
library(stringr)  
p <- DotPlot(sce.all, features = genes_to_check,
             assay='RNA'  )  + coord_flip()

p
ggsave(plot=p, filename="check_Tcells_marker_by_seurat_cluster.pdf")

# mast cells, TPSAB1 and TPSB2 
# B cell,  CD79A  and MS4A1 (CD20) 
# naive B cells, such as MS4A1 (CD20), CD19, CD22, TCL1A, and CD83, 
# plasma B cells, such as CD38, TNFRSF17 (BCMA), and IGHG1/IGHG4
genes_to_check = c('CD3D','MS4A1','CD79A',
                   'CD19', 'CD22', 'TCL1A',  'CD83', #  naive B cells
                   'CD38','TNFRSF17','IGHG1','IGHG4', # plasma B cells,
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'PTPRC' ) 
p <- DotPlot(sce.all, features = genes_to_check,
             assay='RNA'  )  + coord_flip()#B细胞组#

p
ggsave(plot=p, filename="check_Bcells_marker_by_seurat_cluster.pdf")


genes_to_check = c('CD68', 'CD163', 'CD14',  'CD86', 'LAMP3', ## DC 
                   'CD68',  'CD163','MRC1','MSR1','ITGAE','ITGAM','ITGAX','SIGLEC7', 
                   'MAF','APOE','FOLR2','RELB','BST2','BATF3')
p <- DotPlot(sce.all, features = unique(genes_to_check),
             assay='RNA'  )  + coord_flip()#髓系细胞#

p
ggsave(plot=p, filename="check_myeloids_marker_by_seurat_cluster.pdf")


# epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
# - alveolar type I cell (AT1; AGER+)肺泡Ⅰ型细胞
# - alveolar type II cell (AT2; SFTPA1)肺泡II型细胞
# - secretory club cell (Club; SCGB1A1+)分泌俱乐部细胞
# - basal airway epithelial cells (Basal; KRT17+)基底气道上皮细胞
# - ciliated airway epithelial cells (Ciliated; TPPP3+) 纤毛气道上皮细胞
#来自于正常的肺的上皮细胞，可以分成如上所示的5个亚群#

genes_to_check = c(  'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' ,
                     'AGER',
                     'SFTPA1',
                     'SCGB1A1','
                     KRT17',
                     'TPPP3',
                     'KRT4','KRT14','KRT8','KRT18',
                     'CD3D','PTPRC' ) #这些是？#
p <- DotPlot(sce.all, features = unique(genes_to_check),
             assay='RNA'  )  + coord_flip()

p
ggsave(plot=p, filename="check_epi_marker_by_seurat_cluster.pdf")


genes_to_check = c('TEK',"PTPRC","EPCAM","PDPN","PECAM1",'PDGFRB',
                   'CSPG4','GJB2', 'RGS5','ITGA7',
                   'ACTA2','RBP1','CD36', 'ADGRE5','COL11A1','FGF7', 'MME')
p <- DotPlot(sce.all, features = unique(genes_to_check),#骨髓基质细胞#
             assay='RNA'  )  + coord_flip()

p
ggsave(plot=p, filename="check_stromal_marker_by_seurat_cluster.pdf")


p_all_markers
p_umap=DimPlot(sce.all, reduction = "umap",
               group.by = "RNA_snn_res.0.8",label = T) 
library(patchwork)
p_all_markers+p_umap
ggsave('markers_umap.pdf',width = 15)
DimPlot(sce.all, reduction = "umap",split.by = 'orig.ident',
        group.by = "RNA_snn_res.0.8",label = T) 
ggsave('orig.ident_umap.pdf',width = 15)



sce.all  
sce=sce.all
table(Idents(sce))  

library(future)
# check the current active plan
plan()
plan("multiprocess", workers = 8)#不太懂#
#要访问Seurat中的函数的并行版本，需要加载future的包并设置plan。该plan将指定如何执行该函数。默认行为是以非并行的方式(顺序地)计算的。#
plan()


sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)#only.pos = TRUE只返回positive的基因##min.pct = 0.25旨在通过不测试很少表达的基因来加速功能。 默认值为 0.1#
DT::datatable(sce.markers)
 

pro='RNA_snn_res.0.8'
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)#找出每类中表达量最高的前十个基因#
DoHeatmap(sce,top10$gene,size=3)
ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'))

library(dplyr) 
top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
DoHeatmap(sce,top3$gene,size=3)
ggsave(paste0(pro,'DoHeatmap_check_top3_markers_by_clusters.pdf'))
p <- DotPlot(sce, features = unique(top3$gene),
             assay='RNA'  )  + coord_flip()

p
ggsave(paste0(pro,'DotPlot_check_top3_markers_by_clusters.pdf'),
       height = 15,width = 10)

save(sce.markers,file = 'cca-sce.markers.Rdata')
