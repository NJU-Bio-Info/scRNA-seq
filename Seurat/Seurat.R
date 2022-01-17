##Here I use Seurat built-in dataset pbmc3k to perform analysis
library(SeuratData)
library(Seurat)
library(clustree)
library(patchwork)
library(dplyr)

#load data
data("pbmc3k")
pbmc3k = force(pbmc3k)
#check pbmc3k project
pbmc3k

#Here, pbmc3k is already a seurat object
#how many cells and how many genes
dim(pbmc3k)

########Step1: Quality control#########
#calculated the mt genes percentage for every cell
pbmc3k[['percent.mt']] <- PercentageFeatureSet(object = pbmc3k,
                                               pattern = '^MT-')
#or we can:
if(F){
  pbmc3k <- PercentageFeatureSet(object = pbmc3k,
                                 pattern = '^MT-',
                                 col.name = 'percent.mt')
}
#check the mt genes percentage:
head(pbmc3k@meta.data)
#why we need to calculate mt genes percentage:
#Low-quality / dying cells often exhibit extensive mitochondrial contamination.
#check 'nCount_RNA', 'nFeature_RNA' and 'percent.mt' distribution
features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt')
VlnPlot(object = pbmc3k,
        features = features)
#FeatureScatter can be typically used to visualize feature-feature relationships
FeatureScatter(object = pbmc3k,
               feature1 = 'nCount_RNA',
               feature2 = 'nFeature_RNA')
#here, I chose the cell with >200 and <2500 genes and percent.mt<5% to perform downstream analysis
pbmc = subset(x = pbmc3k,
              subset = nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5) #not 0.05!

########Step2: Normalizing the data#########
#score: remove the effect of different sequencing depth for each cell
pbmc <- NormalizeData(object = pbmc,
                      scale.factor = 1e4,
                      normalization.method = 'LogNormalize')
#the results stored in pbmc[["RNA"]]@data

########Step3: Identify variable features (genes)#########
#score: calculate a subset of features that exhibit high cell-to-cell variation in the dataset
pbmc <- FindVariableFeatures(object = pbmc,
                             nfeatures = 2000,
                             selection.method = 'vst')#2000 features to be variable features
#check the result
sum(pbmc@assays[["RNA"]]@meta.features[["vst.variable"]])
(top10 = head(VariableFeatures(object = pbmc), n = 10))
#plot the variable features
VariableFeaturePlot(object = pbmc,
                    cols = c('black', 'red'))#define the color of different type of features
#we can find that the y-axis label is standardized variance
#add features' label to this plot
LabelPoints(plot = VariableFeaturePlot(pbmc), points = top10, repel = T)

########Step4: Scaling the data#########
#score:
#(1) Shifts the expression of each gene, so that the mean expression across cells is 0
#(2) Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#actually, the max value will be set as 10 to remove the effect of outlier
pbmc <- ScaleData(object = pbmc,
                  features = rownames(pbmc))#variable features by default
#The results of this are stored in pbmc[["RNA"]]@scale.data
#check the result:
if(F){
  max.value = apply(X = pbmc@assays$RNA@scale.data, MARGIN = 1, FUN = max)
  head(max.value, 10)
  #here select NOC2L to test
  mean(pbmc@assays$RNA@scale.data[rownames(pbmc@assays$RNA@scale.data)=='NOC2L', ])
  #approximately zero
  sd(pbmc@assays$RNA@scale.data[rownames(pbmc@assays$RNA@scale.data)=='NOC2L', ])
  #standard variance = 1
}



########Step5: Perform linear dimensional reduction##########
pbmc <- RunPCA(object = pbmc,
               features = VariableFeatures(pbmc))#use variable features by default
#pca analysis use scale.data as input!
#so you must perform NormalizeData() and ScaleData() at first
#the main target for pca analysis is to find the variable features that mainly contribute
#to the cell-to-cell heterogeneity, so that we can be computing time and computing space friendly
#some function to visualize the pca analysis output:
print(pbmc@reductions$pca, dims = 1:5, nfeatures = 5)
#here, all of the gene list out is variable features!
VizDimLoadings(object = pbmc, 
               dims = 1:3,#the first 3 pcs
               nfeatures = 10,#the 10 most variable features
               reduction = 'pca')
DimPlot(object = pbmc, dims = c(1,2), reduction = 'pca')#cells' position in pc1 and pc2
DimHeatmap(object = pbmc,
           cells = 500,#500 cells
           nfeatures = 30,#30 variable features
           dims = 1:6)#visualize the first 6 pcs
#DimHeatmap() can help us find the pcs to perform downstream analysis (more heterogeneity)
#you can find that 30th pc may not be too heterogeneous
DimHeatmap(object = pbmc, dims = 30, cells = 500, nfeatures = 30)

########Step6: Determine the ‘dimensionality’ of the dataset#########
#simply speaking, this part mainly want to determine how many pcs for us yo perform downstream analysis
#method 1: from the Step5's DimHeatmap() function's result
#method 2: use JackStraw()
pbmc <- JackStraw(object = pbmc,
                  reduction = 'pca',
                  dims = 20,#estimate the first 20 pCs
                  num.replicate = 100,#perform 100 times
                  prop.freq = 0.01)#subset 1% cells of the dataset
pbmc <- ScoreJackStraw(object = pbmc,
                       dims = 1:20)
JackStrawPlot(object = pbmc,
              dims = 1:20)
#here we can find that, especially the 19/20 PCs is not important
#the reason and some more details can be seen in source.txt
#method3: use ElbowPlot() function
ElbowPlot(object = pbmc,
          ndims = 20)
#the y-axis is the percentage of variance explained by each one PC 
#choose the first 10 PCs to perform down-stream analysis

########Step7: Cluster the cells###########
pbmc <- FindNeighbors(object = pbmc,
                      dims = 1:10)
pbmc <- FindClusters(object = pbmc,
                     resolution = 0.5)#how to choose a suitable resolution value?
#view the pbmc object
View(pbmc)
#check the cluster
table(pbmc@active.ident)
head(Idents(pbmc), 10)

########Step8: Run non-linear dimensional reduction (UMAP/tSNE)#########
#UMAP reduction
pbmc <- RunUMAP(object = pbmc,
                dims = 1:10)
#tSNE reduction
pbmc <- RunTSNE(object = pbmc,
                dims = 1:10)
#visualize
DimPlot(object = pbmc,
        reduction = 'umap')
DimPlot(object = pbmc,
        reduction = 'tsne')
#all the results are stored in pbmc@reductions
########Step9: Identify differentially expressed features (cluster biomarkers)#########
#find all markers for cluster2
cluster2.markers <- FindMarkers(object = pbmc,
                                ident.1 = 2,
                                min.pct = 0.25,
                                #a feature to be detected at a minimum percentage in either of the two groups of cells
                                logfc.threshold = 0.25
                                #max.cells.per.ident = 0.3 #only sample 30% cells to test
                                )
head(cluster2.markers, 10)
#plot a volcano plot
ggplot(data = cluster2.markers, aes(x = avg_log2FC, y = p_val_adj)) +
  geom_point() + 
  theme_test()
#find markers distinguish cluster 5 and cluster 0/1
cluster5.markers <- FindMarkers(object = pbmc,
                                ident.1 = 5,
                                ident.2 = c(0, 1),
                                min.pct = 0.25)
head(cluster5.markers, 10)
#find markers for every cluster compared to all remaining cells, report only the positive
markers <- FindAllMarkers(object = pbmc,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
markers %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 2)
#visualize the markers
#type I: VlnPlot()
VlnPlot(object = pbmc,
        features = c('CCR7', 'S100A9'))
#type II: FeaturePlot()
FeaturePlot(object = pbmc,
            features = c('CCR7', 'S100A9'))
#type III: DoHeatMap()
markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> top10_markers
DoHeatmap(object = pbmc,
          features = top10_markers$gene)#important!


########Step10: Assigning cell type identity to clusters###########
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", 
                     "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")#according to celltype markers
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
#check result
levels(pbmc)
head(pbmc@active.ident)
#visualize
DimPlot(object = pbmc, reduction = 'umap')

########Step11: Perform subset analysis#############
#here, I want to furtherly analysis Naive CD4 T cell cluster
#extract subset
naive_CD4_cluster = subset(x = pbmc,
                           idents = 'Naive CD4 T')
#here we can find that the results produced by previous analysis are still in this object
#how can we clean this object
#Method 1: create new Seurat Object
new_naive_CD4_cluster = CreateSeuratObject(counts = naive_CD4_cluster@assays$RNA@counts)
#then you can rerun the previous pipeline
#Method 2: use DietSeurat() function
clean_naive_CD4_cluster = DietSeurat(object = naive_CD4_cluster,
                                     counts = T,#remain counts data
                                     data = T)#remain data
#recommend Method1