###S1
setwd("~/202210_movAPA/S1")
Braincortex_mt25_Harmony <- readRDS("~/movAPA/Brain_cortex_mt25/Braincortex_mt25_Harmony_dim14r0.1.rds")
####S1A cellratio
Braincortex_mt25_Harmony@meta.data$stim <- Braincortex_mt25_Harmony@meta.data$group
Braincortex_control <- subset(Braincortex_mt25_Harmony,stim=="Non_Viral_Control")
table(Braincortex_control$orig.ident)
prop.table(table(Idents(Braincortex_control)))
table(Idents(Braincortex_control), Braincortex_control$orig.ident)
Cellratio <- prop.table(table(Idents(Braincortex_control), Braincortex_control$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
###
##S1B Heatmap
top20 <- read.table("top.markers_mt25_dim14r0.1.txt",head=TRUE)
top10 <- top20 %>% group_by(cluster_id) %>% top_n(n = 10, wt = avg_fc)
#
set.seed(42)
subobj <- subset(Braincortex_mt25_Harmony, downsample = 1000)
DoHeatmap(subobj,
          features = as.character(top10$gene),
          group.by = "cell_type",
          assay = 'RNA',
          label=F)+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))
###
##F2C marker
#小提琴图
library(Seurat)
library(ggplot2)
features  <-  c("SLC17A7","ST18","GFAP","GAD1","VCAN","CD74","CLDN5")
remotes::install_github("lyc-1995/MySeuratWrappers")#通过链接安装包  
library(MySeuratWrappers)  
my36colors <-c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87", "#E95C59", "#E59CC4", "#AB3282")#颜色设置  
my36colors <-c("#D4807C",
               "#E1A7C9",
               "#926CAC",
               "#8C3B63",
               "#C09343",
               "#06A19F",
               "#126CB5")
VlnPlot(Braincortex_mt25_Harmony, features = features, stacked=T,pt.size=0,
        cols = my36colors, direction = "horizontal", #水平作图
        x.lab = "", y.lab = "")+ #横纵轴不标记任何东西  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())#不显示坐标刻度
### F2D 
library(UCell)
library(irGSEA)
library(AUCell)
Braincortex_mt25_Harmony <- readRDS("~/movAPA/Brain_cortex_mt25/Braincortex_mt25_Harmony_dim14r0.1.rds")
#control
Braincortex_mt25_Harmony@meta.data$stim <- Braincortex_mt25_Harmony@meta.data$group
Braincortex_control <- subset(Braincortex_mt25_Harmony,stim=="Non_Viral_Control")
Braincortex_control <- irGSEA.score(object = Braincortex_control, assay = "RNA",
                                    slot = "data", seeds = 123, ncores = 1,
                                    min.cells = 3, min.feature = 0,
                                    custom = F, geneset = NULL, msigdb = T,
                                    species = "Homo sapiens", category = "H",  
                                    subcategory = NULL, geneid = "symbol",
                                    method = c("AUCell", "UCell", "singscore","ssgsea"),
                                    aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                    kcdf = 'Gaussian')
Seurat::Assays(Braincortex_control)
Braincortex_control.KEGG <- irGSEA.score(object = Braincortex_control, assay = "RNA",
                                         slot = "data", seeds = 123, ncores = 1,
                                         min.cells = 3, min.feature = 0,
                                         custom = F, geneset = NULL, msigdb = T,
                                         species = "Homo sapiens", category = "C2",
                                         subcategory = "CP:KEGG", geneid = "symbol",
                                         method = c("AUCell", "UCell", "singscore","ssgsea"),
                                         aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                         kcdf = 'Gaussian')
Seurat::Assays(Braincortex_control.KEGG)
control.dge.KEGG <- irGSEA.integrate(object = Braincortex_control.KEGG,
                                     group.by = "cell_type",
                                     metadata = NULL, col.name = NULL,
                                     method = c("AUCell","UCell","singscore",
                                                "ssgsea"))
irGSEA.heatmap.control.KEGG <- irGSEA.heatmap(object = control.dge.KEGG,
                                              method = "RRA",
                                              top = 50,
                                              show.geneset = NULL)
irGSEA.heatmap.control.KEGG
#saveRDS(Braincortex_control.KEGG, file = "./Braincortex_control.KEGG.rds")
control.dge <- irGSEA.integrate(object = Braincortex_control,
                                group.by = "cell_type",
                                metadata = NULL, col.name = NULL,
                                method = c("AUCell","UCell","singscore",
                                           "ssgsea"))
#全局展示，热图
irGSEA.heatmap.control <- irGSEA.heatmap(object = control.dge,
                                         method = "RRA",
                                         top = 50,
                                         show.geneset = NULL)
irGSEA.heatmap.control

##################################################################################




