###F1BC motif 
setwd("~/20221004_scAPAtrap/F1")
scPACds <- readRDS("~/20221004_scAPAtrap/20221005_all_scPACds.rds")
#control组(Non_Viral_Control)
scPACds <- subsetPACds(scPACds, group='group',conds="Non_Viral_Control")
table(scPACds@anno$ftr)

library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
##替换部分chr，如GL000194v1.1到chr14_GL000194v1_random
seqnames(BSgenome.Hsapiens.UCSC.hg38)
table(scPACds@anno$chr)
scPACds@anno$chr[which(scPACds@anno$chr=="GL000194.1") ]  <- "chr14_GL000194v1_random"
scPACds@anno$chr[which(scPACds@anno$chr=="GL000195.1") ]  <- "chrUn_GL000195v1"
scPACds@anno$chr[which(scPACds@anno$chr=="GL000218.1") ]  <- "chrUn_GL000218v1"
scPACds@anno$chr[which(scPACds@anno$chr=="GL000219.1") ]  <- "chrUn_GL000219v1"
scPACds@anno$chr[which(scPACds@anno$chr=="KI270711.1") ]  <- "chr1_KI270711v1_random"
scPACds@anno$chr[which(scPACds@anno$chr=="KI270721.1") ]  <- "chr11_KI270721v1_random"
scPACds@anno$chr[which(scPACds@anno$chr=="KI270728.1") ]  <- "chr16_KI270728v1_random"
scPACds@anno$chr[which(scPACds@anno$chr=="KI270734.1") ]  <- "chr22_KI270734v1_random"
table(scPACds@anno$chr)
scPACds@anno$coord <- as.numeric(scPACds@anno$coord)

#'mm' motif顺序与hg一致
v <- getVarGrams('mm')
priority <- c(1,2,rep(3,length(v)-2))
PACdsMM <- annotateByPAS(scPACds, bsgenome, grams=v,
                         priority=priority, 
                         from=-50, to=-1, label='mm')
table(PACdsMM@anno$mm_gram)

#Plot signal logos.
pas <- PACdsMM@anno$mm_gram[!is.na(PACdsMM@anno$mm_gram)]
plotSeqLogo(pas)

#Prepare the data to plot PAS distributions
pas <- as.data.frame(cbind(region=PACdsMM@anno$ftr, PAS=PACdsMM@anno$mm_gram))
pas$PAS[is.na(pas$PAS)] <- 'NOPAS'
pas$PAS[pas$PAS %in% v[-c(1:2)]] <- 'Variants'
n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
#> `summarise()` regrouping output by 'region' (override with `.groups` argument)
n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
#> `summarise()` ungrouping output (override with `.groups` argument)
n=merge(n, n2)
n$PAC=n$nPAC/n$nTot
n=n[n$PAS!='NOPAS', ]
n$PAS=factor(n$PAS, levels=rev(c('AATAAA', 'ATTAAA','Variants', 'NOPAS')))
n$region=factor(n$region, 
                levels=c('3UTR','Ext_3UTR', 'intergenic','intron','CDS','5UTR'))
#Plot PAS distributions.
ggplot(data=n, aes(x=region, y=PAC, fill=PAS)) + 
  geom_bar(stat="identity") + 
  ylab("PAC Fraction") + theme_bw()

################################################################################
###F1 F
##将pasusage导入metadata中
Braincortex_mt25_Harmony_dim14r0.1 <- readRDS("~/movAPA/Brain_cortex_mt25/Braincortex_mt25_Harmony_dim14r0.1.rds")
DimPlot(Braincortex_mt25_Harmony_dim14r0.1, reduction = "harmony_umap", label=T)

Braincortex_mt25_Harmony_dim14r0.1@meta.data$ident <- rownames(Braincortex_mt25_Harmony_dim14r0.1@meta.data) 
Braincortex_mt25_Harmony_dim14r0.1@meta.data$ident <- str_split(Braincortex_mt25_Harmony_dim14r0.1@meta.data$ident,"\\_",simplify = T)[,1]
Braincortex_mt25_Harmony_dim14r0.1@meta.data <- tidyr::unite(Braincortex_mt25_Harmony_dim14r0.1@meta.data, "Ident", orig.ident, ident,sep = "_")
Braincortex_mt25_Harmony_dim14r0.1@meta.data$rowname <- rownames(Braincortex_mt25_Harmony_dim14r0.1@meta.data)
#加入PAS信息
Brain_APA_MT_25_Harmony_dim14r0.1_Proximal <- readRDS("~/movAPA/Brain_APA_MT_25_Harmony_dim14r0.1_Proximal.rds")
Gene_Average_cell_PUI <- Brain_APA_MT_25_Harmony_dim14r0.1_Proximal@meta.data
Gene_Average_cell_PUI$Ident <- rownames(Gene_Average_cell_PUI)
colnames(Gene_Average_cell_PUI)[10:12]
Gene_Average_cell_PUI <- Gene_Average_cell_PUI[,10:12]
Braincortex_mt25_Harmony_dim14r0.1@meta.data <- left_join(Braincortex_mt25_Harmony_dim14r0.1@meta.data,Gene_Average_cell_PUI,by="Ident")#将pasusage导入
rownames(Braincortex_mt25_Harmony_dim14r0.1@meta.data) <- Braincortex_mt25_Harmony_dim14r0.1@meta.data$rowname
Braincortex_mt25_Harmony_dim14r0.1@meta.data <- Braincortex_mt25_Harmony_dim14r0.1@meta.data[,-10]

#boxplot
df <- data.frame(Braincortex_mt25_Harmony_dim14r0.1@meta.data)
df <- df[df$group=="Non_Viral_Control",]
ggplot(df, aes(x=cell_type, y=ProximalRateMean, fill=cell_type)) + 
  geom_boxplot()+
  theme_classic()+
  coord_flip()+
  scale_fill_manual(values=c("#D4807C",
                             "#E1A7C9",
                             "#926CAC",
                             "#8C3B63",
                             "#C09343",
                             "#06A19F",
                             "#126CB5"))
###F1G
##细胞特异性PAS柱状及比例图
cells <- Braincortex_mt25_Harmony_dim14r0.1@meta.data$Ident
Brain_APA <- subset(Brain_APA_MT_25_Harmony_dim14r0.1_Proximal,cells = cells)#根据单细胞细胞ID提取
Non_Viral_Control <- subset(Brain_APA,group=="Non-viral-control")
Non_Viral_Control$Ident <- rownames(Non_Viral_Control@meta.data)
Non_Viral_Control@meta.data <- Non_Viral_Control@meta.data[,-9]
Non_Viral_Control@meta.data <- left_join(Non_Viral_Control@meta.data,df[,c(1,9)],by="Ident")#将单细胞celltype导入
rownames(Non_Viral_Control@meta.data) <- Non_Viral_Control$Ident
#分细胞类型统计PAS表达
#Glutamatergic_neuron
Glutamatergic_neuron <- subset(Non_Viral_Control,cell_type == "Glutamatergic neuron")
object <- Glutamatergic_neuron
group.by <- "cell_type" 
mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
mat <- expm1(x = mat) 
mat <- aggregate(mat, by=list(object@meta.data[[group.by]]), FUN="mean") 
rownames(mat) <- mat$Group.1 
mat <- t(mat[,-1]) 
head(mat)
mat_Glu <- mat
table(mat_Glu==0)
#saveRDS(mat_Glu,file = "./mat_Glu.rds")
#GABAergic interneuron
GABAergic_interneuron <- subset(Non_Viral_Control,cell_type == "GABAergic interneuron")
object <- GABAergic_interneuron
group.by <- "cell_type" 
mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
mat <- expm1(x = mat) 
mat <- aggregate(mat, by=list(object@meta.data[[group.by]]), FUN="mean") 
rownames(mat) <- mat$Group.1 
mat <- t(mat[,-1]) 
head(mat)
mat_GABA <- mat
table(mat_GABA==0)
#saveRDS(mat_GABA,file = "./mat_GABA.rds")
#Oligodendrocyte
Oligodendrocyte <- subset(Non_Viral_Control,cell_type == "Oligodendrocyte")
object <- Oligodendrocyte
group.by <- "cell_type" 
mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
mat <- expm1(x = mat) 
mat <- aggregate(mat, by=list(object@meta.data[[group.by]]), FUN="mean") 
rownames(mat) <- mat$Group.1 
mat <- t(mat[,-1]) 
head(mat)
mat_Oligodendrocyte <- mat
table(mat_Oligodendrocyte==0)
#saveRDS(mat_Oligodendrocyte,file = "./mat_Oligodendrocyte.rds")
#OPC
OPC <- subset(Non_Viral_Control,cell_type == "OPC")
object <- OPC
group.by <- "cell_type" 
mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
mat <- expm1(x = mat) 
mat <- aggregate(mat, by=list(object@meta.data[[group.by]]), FUN="mean") 
rownames(mat) <- mat$Group.1 
mat <- t(mat[,-1]) 
head(mat)
mat_OPC <- mat
table(mat_OPC==0)
#saveRDS(mat_OPC,file = "./mat_OPC.rds")
#Astrocyte
Astrocyte <- subset(Non_Viral_Control,cell_type == "Astrocyte")
object <- Astrocyte
group.by <- "cell_type" 
mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
mat <- expm1(x = mat) 
mat <- aggregate(mat, by=list(object@meta.data[[group.by]]), FUN="mean") 
rownames(mat) <- mat$Group.1 
mat <- t(mat[,-1]) 
head(mat)
mat_Astrocyte <- mat
table(mat_Astrocyte==0)
#saveRDS(mat_Astrocyte,file = "./mat_Astrocyte.rds")
#Microglia
Microglia <- subset(Non_Viral_Control,cell_type == "Microglia")
object <- Microglia
group.by <- "cell_type" 
mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
mat <- expm1(x = mat) 
mat <- aggregate(mat, by=list(object@meta.data[[group.by]]), FUN="mean") 
rownames(mat) <- mat$Group.1 
mat <- t(mat[,-1]) 
head(mat)
mat_Microglia <- mat
table(mat_Microglia==0)
#saveRDS(mat_Microglia,file = "./mat_Microglia.rds")
#Endothelial
Endothelial <- subset(Non_Viral_Control,cell_type == "Endothelial")
object <- Endothelial
group.by <- "cell_type" 
mat <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "RNA", slot = "data")))) 
mat <- expm1(x = mat) 
mat <- aggregate(mat, by=list(object@meta.data[[group.by]]), FUN="mean") 
rownames(mat) <- mat$Group.1 
mat <- t(mat[,-1]) 
head(mat)
mat_Endothelial <- mat
table(mat_Endothelial==0)
#saveRDS(mat_Endothelial,file = "./mat_Endothelial.rds")

##画图
#mat_Microglia <- readRDS("./mat_Microglia.rds")
mat_Glu <- data.frame(mat_Glu)
mat_Glu$gene <- rownames(mat_Glu)

mat_GABA <- data.frame(mat_GABA)
mat_GABA$gene <- rownames(mat_GABA)

mat_Oligodendrocyte <- data.frame(mat_Oligodendrocyte)
mat_Oligodendrocyte$gene <- rownames(mat_Oligodendrocyte)

mat_OPC <- data.frame(mat_OPC)
mat_OPC$gene <- rownames(mat_OPC)

mat_Astrocyte <- data.frame(mat_Astrocyte)
mat_Astrocyte$gene <- rownames(mat_Astrocyte)

mat_Microglia <- data.frame(mat_Microglia)
mat_Microglia$gene <- rownames(mat_Microglia)

mat_Endothelial <- data.frame(mat_Endothelial)
mat_Endothelial$gene <- rownames(mat_Endothelial)
#合并
all <- merge(merge(merge(merge(merge(merge(mat_Glu,mat_GABA,all=T),mat_Oligodendrocyte,all=T),mat_OPC,all=T),
                         mat_Astrocyte,all=T),mat_Microglia,all=T),mat_Endothelial,all=T)
rownames(all) <- all$gene
all <- all[,-1]

#upset
library(UpSetR)
all[all>0] <- 1
upset(all,nset = 7,nintersects = 30)
##堆叠柱状图
##左
library(ggplot2)
library(cols4all)
library(patchwork)
library(ggbreak)
dt <- read.csv('specificpas.csv',header = T)
head(dt)
#转换为因子固定绘图顺序：
dt$Celltype <- factor(dt$Celltype, levels = unique(dt$Celltype))
#自定义色板：
c4a_gui() #挑选色板
mycol <- c4a('set1',2)
#数值型堆叠柱形图绘制：
p <- ggplot(dt,aes(x = Values, y = Celltype, fill = Group)) +
  geom_bar(position = "stack",
           stat="identity",
           alpha = 1) +
  theme_classic()
p
#图表美化：
p3 <- p +
  scale_fill_manual(values = rev(mycol)) + #更换为自定义配色
  theme(legend.position = "top", #图例放到顶部
        legend.title = element_blank(), #去掉图例标题
        axis.title.y = element_blank(),#去掉y轴标题
        axis.text.y = element_blank())+ #去掉y轴标签
  xlab("Number of PACs") #修改x轴标题
p3
#这里为x轴添加“break”；
#注意，这里scale_x_break的expand参数数值会影响“断点”的位置；
p4 <- p3 + scale_x_break(c(8000, 180000),scales = "fixed",
                         expand=expansion(add = c(0, 1000)))
p4
#令scales = "free"或1，使上下两半部分高度相同；
#使用scale_x_continuous定义“断点”之前的标签；
#使用ticklabels参数自定义“断点”之后的标签；
p5 <- p3 + 
  scale_x_break(c(6000, 180000),scales = "free",
                ticklabels=c(180000,240000,300000,360000),
                expand=expansion(add = c(0, 0)))+ 
  theme(legend.position="none")
#或者
p5 <- p3 + scale_x_continuous(limits = c(0, 360000),
                              breaks = c(0,2000,4000,6000),
                              label = c("0","2000","4000","6000"))+
  scale_x_break(c(6000, 180000),scales = "free",
                ticklabels=c(180000,240000,300000,360000),
                expand=expansion(add = c(0, 0)))+ 
  theme(legend.position="right")
p5
#计算特异性在总的中的占比：
dt1 <- read.csv('dt1.csv',header = T)
#dt$Fraction <- dt$Values/sum(dt$Values) 算的不是spefic的比例
head(dt1)
#转换为因子固定绘图顺序：
dt1$Celltype <- factor(dt1$Celltype, levels = unique(dt1$Celltype))
#绘制柱形图：
p6 <- ggplot(dt1,aes(x = Fraction, y = Celltype)) +
  geom_col(fill = c('#E41A1C')) +
  theme_classic()
p6
#图表美化/修改：
p7 <- p6 +
  scale_x_continuous(expand = c(0,0),limits = c(0.000,1.600)) +
  theme(axis.title.y = element_blank())+#去掉y轴标题 
  xlab("Fraction of cell type specific PACs (%)") #修改x轴标题
p7

p8 <- p7 + 
  scale_x_break(c(0.150, 1.200),scales = "free",
                ticklabels=c(1.200,1.600),
                expand=expansion(add = c(0, 0)))
#将Y轴放到右侧：
p9 <- p8 +
  scale_y_discrete(position = "right") + #y轴逆转
  scale_x_reverse(expand = c(0,0)) #数值也需要同时逆转
p9
p9+p5

###F1H S1
#神经元高表达PAS，功能分析
A <- Non_Viral_Control
A@meta.data$Celltype <- "Others"
A@meta.data[A@meta.data$cell_type=='Glutamatergic neuron'|A@meta.data$cell_type=='GABAergic interneuron',]$Celltype <- "Neuron"
Idents(A)="Celltype"
deg <- FindMarkers(A,ident.1 = "Neuron" ,ident.2 = "Others")
deg1 <- deg %>% subset(p_val<0.05) %>% subset(avg_log2FC>1)
deg1$PAS <- rownames(deg1)
deg1$gene <- str_split(deg1$PAS,"\\-",simplify = T)[,1]  
gene <- unique(deg1$gene)  #去重后153
#write.csv(gene,"degPAC.csv")

#MEGS PACs批量小提琴图  
MEG3_genes <- rownames(A@assays$RNA)[grep("^MEG3-",rownames(A@assays$RNA))]
P <- VlnPlot(A, features = MEG3_genes)
vln.df <- as.data.frame(A[["RNA"]]@data[MEG3_genes,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("Ident","exp")
head(vln.df)
anno=A@meta.data[,c("Ident","Celltype")]
head(anno)
vln.df=inner_join(vln.df,anno,by="Ident")
vln.df$gene=factor(vln.df$gene,levels = MEG3_genes) #为了控制画图的基因顺序
head(vln.df)
##小提琴图，横向
vln.df%>%ggplot(aes(Celltype,exp))+geom_violin(aes(fill=Celltype),scale = "width")+
  facet_grid(.~gene,scales = "free_y")+
  scale_fill_manual(values=c("#E41A1C","#377EB8"))+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )+
  ylab("Expression of PACs")

