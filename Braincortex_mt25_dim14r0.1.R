setwd("~/Brain_cortex_mt25")
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
Brain_cortex_mt25 <- readRDS("~/Brain_cortex_mt25/Brain_cortex_mt25.rds")
Braincortex_mt25 <- merge(Brain_cortex_mt25[[1]], y=c(Brain_cortex_mt25[[2]], Brain_cortex_mt25[[3]], 
                                            Brain_cortex_mt25[[4]], Brain_cortex_mt25[[5]], 
                                            Brain_cortex_mt25[[6]], Brain_cortex_mt25[[7]], 
                                            Brain_cortex_mt25[[8]], Brain_cortex_mt25[[9]], 
                                            Brain_cortex_mt25[[10]],Brain_cortex_mt25[[11]],
                                            Brain_cortex_mt25[[12]],Brain_cortex_mt25[[13]],
                                            Brain_cortex_mt25[[14]],Brain_cortex_mt25[[15]]))
dim(Braincortex_mt25)
table(Braincortex_mt25@meta.data$orig.ident) 
Braincortex_mt25 <- NormalizeData(Braincortex_mt25) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
DimPlot(Braincortex_mt25, reduction = "pca")
#harmony整合
Braincortex_mt25_Harmony <- RunHarmony(Braincortex_mt25, group.by.vars = "orig.ident")
ElbowPlot(Braincortex_mt25_Harmony, ndims = 50)
##################
#dim：作者提出一个确定PC阈值的三个标准：主成分累积贡献大于90%、PC本身对方差贡献小于5%、两个连续PCs之间差异小于0.1%
# Determine percent of variation associated with each PC
pct <- Braincortex_mt25_Harmony [["pca"]]@stdev / sum( Braincortex_mt25_Harmony [["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
# Elbow plot to visualize  #算出14
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
###################
Braincortex_mt25_Harmony <- RunUMAP(Braincortex_mt25_Harmony, reduction = "harmony", dims = 1:14,reduction.name = 'harmony_umap', reduction.key = 'harmony_umap')
Braincortex_mt25_Harmony <- FindNeighbors(Braincortex_mt25_Harmony, reduction = "harmony", dims = 1:14) %>% FindClusters(resolution = 0.1)
#group_by_cluster
plot1 = DimPlot(Braincortex_mt25_Harmony, reduction = "harmony_umap", label=T) 
#group_by_sample
plot2 = DimPlot(Braincortex_mt25_Harmony, reduction = "harmony_umap", group.by='orig.ident') 
#combinate
plotc <- plot1+plot2

###随着resolution的调高，具体是哪些群在不停地继续裂变成为多个群
#先运行Braincortex_mt25_Harmony <- FindNeighbors(Braincortex_mt25_Harmony, reduction = "harmony", dims = 1:14)
# Check clustering stability at given resolution  
# Set different resolutions 
res.used <- seq(0.1,0.5,by=0.1)
res.used
# Loop over and perform clustering of different resolutions 
for(i in res.used){
  sce <- FindClusters(object = Braincortex_mt25_Harmony, verbose = T, resolution = res.used)
}
# Make plot 
library(clustree)
clus.tree.out <- clustree(sce) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")

clus.tree.out
###

#识别保守的细胞类型标记
DefaultAssay(Braincortex_mt25_Harmony) <- "RNA"
cluster0_conserved_markers <- FindConservedMarkers(Braincortex_mt25_Harmony,
                                                   ident.1 = 0,
                                                   grouping.var = "group",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
get_conserved <- function(cluster){
  FindConservedMarkers(Braincortex_mt25_Harmony,
                       ident.1 = cluster,
                       grouping.var = "group",
                       only.pos = TRUE) %>%
    cbind(cluster_id = cluster, .)
}
conserved_markers <- map_dfr(0:9, get_conserved)
# 提取每个类群排名靠前的20个标记物
gene <- row.names(conserved_markers)
conserved_markers <- cbind(gene,conserved_markers)
top20 <- conserved_markers %>% 
  mutate(avg_fc = (Non_Viral_Control_avg_log2FC + Covid19_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)
top20_gene <- str_split(top20$gene,'[...]',simplify = T)[,1]
top20_1 <- cbind(top20_gene,top20[,-1])
names(top20_1)[names(top20_1) == '...1'] <- 'gene'
View(top20_1)

write.table(top20_1,file = "./top.markers_mt25_dim14r0.1.txt",sep = "\t")

features  <-  c("SLC17A7","ST18","GFAP","GAD1","VCAN","CD74","CLDN5")
FeaturePlot(Braincortex_mt25_Harmony,features = features,min.cutoff = "q9")
DotPlot(Braincortex_mt25_Harmony, features = features, col.min = 0)
DotPlot(Braincortex_mt25_Harmony, features = features,cols = c("blue", "red"),dot.scale = 8,split.by = "group") +RotatedAxis()
RidgePlot(Braincortex_mt25_Harmony, features,group.by = 'seurat_clusters', ncol = 2)
VlnPlot(Braincortex_mt25_Harmony, features = features,pt.size=0)
#热图，采取随机取样的方式将细胞数量downsample到30k以下
# 给所有基因都做一下scale
Braincortex <- ScaleData(Braincortex_mt25_Harmony)
# 每个cluster选取300个细胞来画热图
set.seed(42)
subobj <- subset(Braincortex, downsample = 300)
DoHeatmap(subobj, features = top20$gene,group.by = 'seurat_clusters')

# 注释
new.cluster.ids <- c(
  "Glutamatergic neuron", # 00 02 03
  "Oligodendrocyte", # 01
  "Glutamatergic neuron", #02
  "Glutamatergic neuron", #03
  "Astrocyte", # 04
  "GABAergic interneuron", # 05 06
  "GABAergic interneuron", # 06
  "OPC", # 07
  "Microglia", # 08
  "Endothelial" # 09
)
names(new.cluster.ids) <- levels(Braincortex_mt25_Harmony)
Braincortex_mt25_Harmony <- RenameIdents(Braincortex_mt25_Harmony, new.cluster.ids)
Braincortex_mt25_Harmony$cell_type <- Idents(Braincortex_mt25_Harmony)
DimPlot(Braincortex_mt25_Harmony, reduction = "harmony_umap", label=T) 

#细胞比例
table(Braincortex_mt25_Harmony$orig.ident)##查看各样本细胞数
prop.table(table(Idents(Braincortex_mt25_Harmony)))#各细胞百分比
table(Idents(Braincortex_mt25_Harmony), Braincortex_mt25_Harmony$orig.ident)#各样本不同细胞群细胞数
Cellratio <- prop.table(table(Idents(Braincortex_mt25_Harmony), Braincortex_mt25_Harmony$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.8,size = 0.2,colour = 'black')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  scale_fill_manual(values=c("#D4807C",
                             "#E1A7C9",
                             "#926CAC",
                             "#8C3B63",
                             "#C09343",
                             "#06A19F",
                             "#126CB5"))

saveRDS(Braincortex_mt25_Harmony, file = "./Braincortex_mt25_Harmony_dim14r0.1.rds")
#分组细胞比例
table(Braincortex_mt25_Harmony$group)##查看各组细胞数
table(Idents(Braincortex_mt25_Harmony), Braincortex_mt25_Harmony$group)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(Braincortex_mt25_Harmony), Braincortex_mt25_Harmony$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio1 <- prop.table(table(Idents(Braincortex_mt25_Harmony)))#各细胞百分比
Cellratio <- as.data.frame(Cellratio)
Cellratio1 <- as.data.frame(Cellratio1)
Cellratio1$Var2 <- rep("All",times=7)
Cellratio<-rbind(Cellratio,Cellratio1)#合并，总，control，case组
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.2,colour = 'white')+ 
  theme_classic() +
  labs(x='Group',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  scale_fill_manual(values=c("#D4807C",
                             "#E1A7C9",
                             "#926CAC",
                             "#8C3B63",
                             "#C09343",
                             "#06A19F",
                             "#126CB5"))
