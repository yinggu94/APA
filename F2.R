###F2
##F2E
library('ggalluvial')
sang_tot <- read.csv("./sang_tot_1201.csv")
sang_tot$regulation <- "N"
sang_tot[sang_tot$avg_log2FC  > 0.25,]$regulation <- "Up"
sang_tot[sang_tot$avg_log2FC  < (-0.25) ,]$regulation <- "Down"
name <- c("subtype","PAStype","cell","regulation")
sang_tot1 <- sang_tot[,name]
sang_tot2 <- sang_tot1 %>% group_by(subtype,PAStype,cell,regulation) %>% summarise(n = n())
df_lodes <- to_lodes_form(sang_tot2,key ="x", value = "stratum", id = "alluvium",axes =1:3)
head(df_lodes,12)
p1 <- ggplot(df_lodes,
             aes(x = x, y = n, stratum =stratum, alluvium = alluvium,fill = regulation,label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_flow(width = 0.2, knot.pos = 0.1) +
  geom_stratum(alpha = 1,color="white",width = 1/6,fill = "#DCDCDC") +
  geom_text(stat = "stratum", size =2,color="black") +
  scale_fill_manual(values = c("#377EB8","#E41A1C")) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),axis.ticks =element_blank(),axis.text.y =element_blank())
p1
##F2F
library(ggplot2)
sang_tot <- sang_tot[sang_tot$subtype=='3UTR',] # 1036
name <- c("cell","PAStype","regulation")
sang_tot1 <- sang_tot[,name]
sang_tot2 <- sang_tot1 %>% group_by(cell,PAStype,regulation) %>% summarise(n = n())
ggplot(sang_tot2,aes(x=PAStype,y=n,fill=regulation))+
  geom_bar(stat="identity",position = "stack")+
  facet_wrap(~cell, nrow = 1)+
  scale_fill_manual(values = c("#377EB8","#E41A1C")) +
  theme_bw()
##F2H
sang_tot <- sang_tot[sang_tot$subtype=='3UTR',] # 1036
r.deg1 <- sang_tot
data_up <- r.deg1[r.deg1$avg_log2FC > 0.25,] # 370
data_down <-r.deg1[r.deg1$avg_log2FC <= -0.25,] # 666
#up
up <- as.data.frame.array(table(data_up$pas,data_up$cell))
up$Microglia <- 0
up$Endothelial <- 0
up <- up[,c(3,4,1,2,5,6,7)]
up$sum <- apply(up,1,sum) # 298
uniq <- up[up$sum ==1,c(-8)] # 251
same <- up[up$sum !=1,c(-8)] # 47

#complexheatmap
library(pheatmap)
library(ComplexHeatmap)
p1 <- Heatmap(same,cluster_columns = F,show_row_names = F,border = T,col = c("#ffffff","#E41A1C"),show_row_dend = F,show_heatmap_legend = F)
p2 <- Heatmap(uniq,cluster_columns = F,show_row_names = F,border = T,col = c("#ffffff","#E41A1C"),show_row_dend = F,show_heatmap_legend = F)
p1_p2 = p1 %v% p2
draw(p1_p2)
#统计
uniq_up <- colSums(uniq) #上调的基因中每个细胞类型中特有的总数
total_uniq_up <- sum(unlist(uniq)) #上调的基因中每个细胞类型特有的基因数目
same_up <- colSums(same) #上调的基因中每个细胞类型中共有的总数
total_same_up <- sum(unlist(same)) #上调的基因中每个细胞类型中共有的基因数目
#down
down <- as.data.frame.array(table(data_down$pas,data_down$cell))
down <- down[,c(4,6,1,3,7,5,2)]
down$sum <- apply(down,1,sum)
uniq2 <- down[down$sum ==1,c(-8)] # 321
same2 <- down[down$sum !=1,c(-8)] # 156
p1 <- Heatmap(same2,cluster_columns = F,show_row_names = F,border = T,col = c("#ffffff","#377EB8"),show_row_dend = F,show_heatmap_legend = F)
p2 <- Heatmap(uniq2,cluster_columns = F,show_row_names = F,border = T,col = c("#ffffff","#377EB8"),show_row_dend = F,show_heatmap_legend = F)
p1_p2 = p1 %v% p2
draw(p1_p2)
uniq_down <- colSums(uniq2)
total_uniq_down <- sum(unlist(uniq2))
same_down <- colSums(same2)
total_same_down <- sum(unlist(same2))
##F2H
#distal
data_up_distal <- data_up[data_up$PAStype=="Distal",] # 310
data_up_proximal <- data_up[data_up$PAStype=="Proximal",] # 60
data_down_distal <- data_down[data_down$PAStype=="Distal",] # 580
data_down_proximal <- data_down[data_down$PAStype=="Proximal",] # 86
#up
pas_up <- unique(data_up_distal[,c("genename","cell")]) # 292
table(pas_up$cell)
pas_down <- unique(data_down_distal[,c("genename","cell")]) # 482 
table(pas_down$cell)
#BP
library("clusterProfiler")
library(org.Hs.eg.db)
library(stringr)
library(enrichplot)
type <- c("Astrocyte","GABAergic interneuron","Glutamatergic neuron","Oligodendrocyte","OPC")
for (i in 1:5) {
  input <- pas_up[pas_up$cell==type[i],]
  inputgene <- bitr(input$genename, fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  GO_BP <- enrichGO(inputgene$ENTREZID,
                    OrgDb = 'org.Hs.eg.db',
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.1,
                    readable = T)
  write.csv(GO_BP, file = paste0(type[i],'.bp.csv') )
}
#down
type <- c("Astrocyte","GABAergic interneuron","Glutamatergic neuron","Oligodendrocyte","OPC","Endothelial")
for (i in 1:6) {
  input <- pas_down[pas_up$cell==type[i],]
  inputgene <- bitr(input$genename, fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  GO_BP <- enrichGO(inputgene$ENTREZID,
                    OrgDb = 'org.Hs.eg.db',
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.1,
                    readable = T)
  write.csv(GO_BP, file = paste0(type[i],'.bp.csv') )
}
#proximal
#up
pas_up <- unique(data_up_proximal[,c("genename","cell")]) # 60
table(pas_up$cell)
pas_down <- unique(data_down_proximal[,c("genename","cell")]) # 86 
table(pas_down$cell)
#up
type <- c("Astrocyte","GABAergic interneuron","Glutamatergic neuron","Oligodendrocyte")
for (i in 1:4) {
  input <- pas_up[pas_up$cell==type[i],]
  inputgene <- bitr(input$genename, fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  GO_BP <- enrichGO(inputgene$ENTREZID,
                    OrgDb = 'org.Hs.eg.db',
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.5,
                    readable = T)
  write.csv(GO_BP, file = paste0(type[i],'.bp.csv') )
}

#down
type <- c("Astrocyte","GABAergic interneuron","Glutamatergic neuron","Oligodendrocyte")
for (i in 1:4) {
  input <- pas_down[pas_up$cell==type[i],]
  inputgene <- bitr(input$genename, fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  GO_BP <- enrichGO(inputgene$ENTREZID,
                    OrgDb = 'org.Hs.eg.db',
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.5,
                    readable = T)
  write.csv(GO_BP, file = paste0(type[i],'.bp.csv') )
}
#气泡矩阵图
#distal.up
library(reshape2)
library(ggplot2)
distal.up.heatmap <- read.csv("~/20221004_scAPAtrap/F2/distalupheatmap.csv")
distal.up.heatmap1 <- melt(distal.up.heatmap)
names(distal.up.heatmap1) = c('BP','Celltype','Value')
distal.up.pvalue <- read.csv("~/20221004_scAPAtrap/F2/distaluppvalue.csv")
distal.up.pvalue1 <- melt(distal.up.pvalue)
names(distal.up.pvalue1) = c('BP','Celltype','Pvalue')

A <- merge(distal.up.heatmap1,distal.up.pvalue1,sort = F)

A$BP <- factor(A$BP, levels = distal.up.heatmap$X)
A$Celltype <- factor(A$Celltype, levels = c("Glu","GABA","Ast","OPC"))
A[order(A$BP,A$Celltype),]
A[A == 0] <- NA

ggplot(A,aes(x = Celltype, y = BP))+
  geom_point(aes(color=Pvalue,size=Value))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "darkgrey",size=1.5,fill=NA),
        axis.text.x=element_text(angle = 45,hjust = 1))+
  scale_color_gradient(low="#CD2626",high="#D2752C",na.value = "white")+
  scale_size_continuous(range=c(3,7))+
  theme(text=element_text(face="bold",color = "black", size=10),
        plot.title = element_text(size=10, face="bold",hjust = 0.5),
        legend.position = "left")+
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Upregulated")

#distal.down
distal.down.heatmap <- read.csv("~/20221004_scAPAtrap/F2/distaldownheatmap.csv")
distal.down.heatmap1 <- melt(distal.down.heatmap)
names(distal.down.heatmap1) = c('BP','Celltype','Value')
distal.down.pvalue <- read.csv("~/20221004_scAPAtrap/F2/distaldownpvalue.csv")
distal.down.pvalue1 <- melt(distal.down.pvalue)
names(distal.down.pvalue1) = c('BP','Celltype','Pvalue')

B <- merge(distal.down.heatmap1,distal.down.pvalue1,sort = F)

B$BP <- factor(B$BP, levels = distal.down.heatmap$X)
B$Celltype <- factor(B$Celltype, levels = c("Glu","Ast","Olig"))
B[order(B$BP,B$Celltype),]
B[B == 0] <- NA

ggplot(B,aes(x = Celltype, y = BP))+
  geom_point(aes(color=Pvalue,size=Value))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "darkgrey",size=1.5,fill=NA),
        axis.text.x=element_text(angle = 45,hjust = 1))+
  scale_color_gradient(low="#4472C4",high="#42A4B7",na.value = "white")+
  scale_size_continuous(range=c(3,7))+
  scale_y_discrete(position = "right")+
  theme(text=element_text(face="bold",color = "black", size=10),
        plot.title = element_text(size=10, face="bold",hjust = 0.5),
        legend.position = "right")+
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Downregulated")
#proximal.up
proximal.up.heatmap <- read.csv("~/20221004_scAPAtrap/F2new2/proximalupheatmap.csv")
proximal.up.heatmap1 <- melt(proximal.up.heatmap)
names(proximal.up.heatmap1) = c('BP','Celltype','Value')
proximal.up.pvalue <- read.csv("~/20221004_scAPAtrap/F2new2/proximaluppvalue.csv")
proximal.up.pvalue1 <- melt(proximal.up.pvalue)
names(proximal.up.pvalue1) = c('BP','Celltype','Pvalue')
C <- merge(proximal.up.heatmap1,proximal.up.pvalue1,sort = F)
C$BP <- factor(C$BP, levels = proximal.up.heatmap$X)
C$Celltype <- factor(C$Celltype, levels = c("Glu","Ast","Olig"))
C[order(C$BP,C$Celltype),]
C[C == 0] <- NA
ggplot(C,aes(x = Celltype, y = BP))+
  geom_point(aes(color=Pvalue,size=Value))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "darkgrey",size=1.5,fill=NA),
        axis.text.x=element_text(angle = 45,hjust = 1))+
  scale_color_gradient(low="#CD2626",high="#D2752C",na.value = "white")+
  scale_size_continuous(range=c(5,9))+
  theme(text=element_text(face="bold",color = "black", size=12),
        plot.title = element_text(size=10, face="bold",hjust = 0.5),
        legend.position = "left")+
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Upregulated")

#proximal.down
proximal.down.heatmap <- read.csv("~/20221004_scAPAtrap/F2new2/proximaldownheatmap.csv")
proximal.down.heatmap1 <- melt(proximal.down.heatmap)
names(proximal.down.heatmap1) = c('BP','Celltype','Value')
proximal.down.pvalue <- read.csv("~/20221004_scAPAtrap/F2new2/proximaldownpvalue.csv")
proximal.down.pvalue1 <- melt(proximal.down.pvalue)
names(proximal.down.pvalue1) = c('BP','Celltype','Pvalue')
D <- merge(proximal.down.heatmap1,proximal.down.pvalue1,sort = F)
D$BP <- factor(D$BP, levels = proximal.down.heatmap$X)
D$Celltype <- factor(D$Celltype, levels = c("Glu","GABA","Ast","Olig"))
D[order(D$BP,D$Celltype),]
D[D == 0] <- NA
ggplot(D,aes(x = Celltype, y = BP))+
  geom_point(aes(color=Pvalue,size=Value))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "darkgrey",size=1.5,fill=NA),
        axis.text.x=element_text(angle = 45,hjust = 1))+
  scale_color_gradient(low="#4472C4",high="#42A4B7",na.value = "white")+
  scale_size_continuous(range=c(5,9))+
  scale_y_discrete(position = "right")+
  theme(text=element_text(face="bold",color = "black", size=12),
        plot.title = element_text(size=10, face="bold",hjust = 0.5),
        legend.position = "right")+
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Downregulated")