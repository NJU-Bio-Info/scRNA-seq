rm(list=ls()) #清空历史变量#
options(stringsAsFactors = F) #在读入数据时，遇到字符串之后，不将其转换为factors，仍然保留为字符串格式#
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)

###### step1:导入数据 ######     
library(data.table)
dir='E-MTAB-9265/outputs' 
samples=list.files(dir)
samples 
# matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv  
sceList = lapply(samples,function(pro){ #function用来自己定义函数#
  # pro=samples[1]
  folder=file.path( dir ,pro) #file.path（a, b）把a与b用“/”连接起来形成一个路径#
  print(pro)
  print(folder)
  print(list.files(folder))
  sce=CreateSeuratObject(counts = Read10X(folder),
                         project =  pro ,
                         min.cells = 5,
                         min.features = 300)#nFeature_RNA # 每个细胞所检测到的基因数目#
  
  return(sce)
})
names(sceList)  
samples  
names(sceList)  = samples 
sce.all=merge(x=sceList[[1]],#直接提取第1个子表中的所有元素#
              y=sceList[ -1 ],#表示输入剩下的多个#
              add.cell.ids = samples)#用add.cell.id参数将样品特异的前缀添加到每个细胞ID，用以区分不同样本中的同个细胞ID#

as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])#counts看的是每个基因在不同细胞里的表达量#
head(sce.all@meta.data, 10)#meta.data可以看每个细胞中基因和RNA的表达量#
table(sce.all$orig.ident)#orig.ident：通常包含所知的样品名，默认为我们赋给project的值##table函数可用于统计因子数量#

if(F){
  
  library(clusterProfiler)
  # BiocManager::install('org.Cf.eg.db')
  library(org.Cf.eg.db)
  head(rownames(sce.all))
  tmp=as.data.frame(rownames(sce.all))
  ids=bitr(rownames(sce.all),'ENSEMBL','SYMBOL',org.Cf.eg.db)
  ids=ids[!duplicated(ids$SYMBOL),]
  head(ids)
  pos=match(ids$ENSEMBL,rl$V1)#match函数获取括号中左边的在右边的中的位置信息
  ct <- mtx[pos,] 
  rownames(ct) <- ids$SYMBOL
  ct[1:4,1:4]
  colnames(ct)=cl$V1
  sce.all=CreateSeuratObject(counts = ct,
                             project = 'E-MTAB-7427')
  
  as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
  head(sce.all@meta.data, 10)
  table(sce.all@meta.data$orig.ident) 
}

  

#rm(list=ls())
###### step2:QC质控 ######
dir.create("./1-QC")
setwd("./1-QC")
# sce.all=readRDS("../sce.all_raw.rds")
#计算线粒体基因比例
# 人和鼠的基因名字稍微不一样 
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all))] #rownames就是可以提取数据框的第一列,grep返回坐标#
mito_genes #13个线粒体基因
sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)

#计算核糖体基因比例
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
ribo_genes
sce.all=PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")#col.name就是要分配的 meta.data 列中的名称。#
fivenum(sce.all@meta.data$percent_ribo)
#计算红血细胞基因比例
rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
sce.all=PercentageFeatureSet(sce.all, "^HB[^(P)]", col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1
library(ggplot2) 
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2	
ggsave(filename="Vlnplot2.pdf",plot=p2)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
ggsave(filename="Scatterplot.pdf",plot=p3)
#根据上述指标，过滤低质量细胞/基因#咋根据？#
#过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 300)
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0 ) > 3]
#sce.all@assays$RNA@counts > 0 返回的结果仍然是一个matrix，只不过这个时候的matrix元素只有TRUE和FALSE两个然后再对这个matrix进行按行求和（Matrix::rowSums()）就能知道一个基因在多少个细胞中有表达#

sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)#一个过滤函数#
dim(sce.all) #R语言中的dim()函数用于获取或设置指定矩阵、数组或 DataFrame 的维数##比如这里结果是10813 4124 就说明有10813行4124列#
dim(sce.all.filt) 
#  可以看到，主要是过滤了基因，其次才是细胞#嗯嗯是这样吗？#

# par(mar = c(4, 8, 2, 1))
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
# 这里的C 这个矩阵，有一点大，可以考虑随抽样
C=C[,sample(1:ncol(C),1000)]
most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
pdf("TOP50_most_expressed_gene.pdf",width=14)
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(50)[50:1], 
        horizontal = TRUE)
se()
rm(C)#rm () R语言中的函数用于从内存中删除对象#

#过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 50)#whichcells函数返回向量，相当于记下有哪些细胞#
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 1 )
length(selected_hb)#length函数用于获取其他向量的长度#
length(selected_ribo)
length(selected_mito)


sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)

table(sce.all.filt$orig.ident) 

#可视化过滤后的情况
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered)

#过滤指标3:过滤特定基因
# Filter MALAT1 管家基因#为什么要过滤管家基因，是因为管家基因哪哪都表达所以不重要是不#
sce.all.filt <- sce.all.filt[!grepl("MALAT1", rownames(sce.all.filt),ignore.case = T), ]
# Filter Mitocondrial 线粒体基因
sce.all.filt <- sce.all.filt[!grepl("^MT-", rownames(sce.all.filt),ignore.case = T), ]
# 当然，还可以过滤更多

dim(sce.all.filt) 

#细胞周期评分#基于G2/M和S期的经典marker基因的表达，计算每个细胞可能所处的细胞周期的分数#
sce.all.filt = NormalizeData(sce.all.filt)#这样就去除了测序深度的不同对结果造成的影响#
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce.all.filt=CellCycleScoring(object = sce.all.filt, 
                              s.features = s.genes, 
                              g2m.features = g2m.genes, 
                              set.ident = TRUE)
p4=VlnPlot(sce.all.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
           ncol = 2, pt.size = 0)
ggsave(filename="Vlnplot4_cycle.pdf",plot=p4)

sce.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
ggsave(filename="cycle_details.pdf" )
# S.Score较高的为S期，G2M.Score较高的为G2M期，都比较低的为G1期
 
dim(sce.all) #就是莫得过滤过的#

#saveRDS(sce.all.filt, "sce.all_qc.rds")
