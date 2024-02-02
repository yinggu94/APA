###F4
#F4A 计每种细胞类型的miRNAs数目
Glutamatergic_neuron <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Glutamatergic neuron.csv")
GABAergic_interneuron <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/GABAergic interneuron.csv")
Astrocyte <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Astrocyte.csv")
Oligodendrocyte <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Oligodendrocyte.csv")
OPC <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/OPC.csv")
Microglia <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Microglia.csv")
Endothelial <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Endothelial.csv")
celltype <- list(Glutamatergic_neuron,GABAergic_interneuron,Astrocyte,
                 Oligodendrocyte,OPC,Microglia,Endothelial)
type <- c("Glutamatergic_neuron","GABAergic_interneuron","Astrocyte",
          "Oligodendrocyte","OPC","Microglia","Endothelial")
Line5 <- Line[,1:2]
all <- data.frame()
for (i in 1:7) {
  input <- celltype[[i]] %>% drop_na()
  A <- subset(Line5,Line5$Gene %in% unlist(input$hgnc_symbol))
  B <- celltype[[i]][,c(8,13)] %>% distinct()
  colnames(B) <- c("Gene","flag")
  degs.miR <- left_join(A,B,by='Gene')
  write.csv(degs.miR,file = paste0(type[i],'.Disease.csv') )
  C <- unique(degs.miR[,2:3])
  D <- data.frame(table(C$flag))  
  D$celltype <- type[i]
  all <- rbind(all,D)
}
#write.csv(all, file = "miRNA_duiji.csv")
all$celltype <- factor(all$celltype,
                       levels = c("Glutamatergic_neuron","GABAergic_interneuron","Astrocyte",
                                  "OPC","Oligodendrocyte","Microglia","Endothelial"))
ggplot(all,aes(x=celltype,y=Freq,fill=Var1))+
  geom_bar(stat = "identity",
           position = "stack",
           colour = 'white')+
  scale_fill_manual(values = alpha(c("1"="tomato","2"="steelblue1",
                                     "3"="steelblue","4"="red4"),0.8))+
  coord_polar()+
  ylim(-1000,NA)+
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

#F4B 统计每种细胞类型的miRNAs热图
i <- 1
input <- celltype[[i]] %>% drop_na()
A <- subset(Line5,Line5$Gene %in% unlist(input$hgnc_symbol))
B <- celltype[[i]][,c(8,13)] %>% distinct()
colnames(B) <- c("Gene","flag")
degs.miR <- left_join(A,B,by='Gene')
degs.miR[degs.miR$flag==2|degs.miR$flag==3,]$flag <- -1
degs.miR[degs.miR$flag==4,]$flag <- 1
C <- unique(degs.miR[,2:3])
all.heatmap1 <- C[C$flag==1,]
colnames(all.heatmap1)[2] <- paste0(type[i],'_distal')
all.heatmap2 <- C[C$flag== -1,]
colnames(all.heatmap2)[2] <- paste0(type[i],'_proximal')
all.heatmap <- full_join(all.heatmap1,all.heatmap2)
for (i in 2:7) {
  input <- celltype[[i]] %>% drop_na()
  A <- subset(Line5,Line5$Gene %in% unlist(input$hgnc_symbol))
  B <- celltype[[i]][,c(8,13)] %>% distinct()
  colnames(B) <- c("Gene","flag")
  degs.miR <- left_join(A,B,by='Gene')
  degs.miR[degs.miR$flag==2|degs.miR$flag==3,]$flag <- -1
  degs.miR[degs.miR$flag==4,]$flag <- 1
  C <- unique(degs.miR[,2:3])
  D1 <- C[C$flag==1,]
  colnames(D1)[2] <- paste0(type[i],'_distal')
  D2 <- C[C$flag== -1,]
  colnames(D2)[2] <- paste0(type[i],'_proximal')
  all.heatmap <- full_join(all.heatmap,D1)
  all.heatmap <- full_join(all.heatmap,D2)
}
#i为6，7,补上
i <- 6
i <- 7

rownames(all.heatmap) <- all.heatmap$Target
miR.heatmap <- all.heatmap[,-1]
miR.heatmap[is.na(miR.heatmap)] <- 0
miR.heatmap <- t(miR.heatmap)
miR.heatmap <- miR.heatmap[c(1,3,5,7,9,11,13,2,4,6,8,10),]
pheatmap(miR.heatmap,
         scale = "none",
         color = colorRampPalette(colors = c("#4472C4","white","#CD2626"))(3),
         use_raster = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,treeheight_col = 0,
         gaps_row = 7)
#CD
pie <- data.frame(t(miR.heatmap)) # 2403
pie$sum1 <- apply(pie[,1:7], 1, sum)
pie$sum2 <- apply(pie[,8:12], 1, sum)
pie <- pie[,13:14]
length(which(pie$sum1 == 0)) # 171 distal
length(which(pie$sum2 == 0)) # 260 proximal
#修改
pie <- all.heatmap[,-1]
pie[is.na(pie)] <- 0
Glu.pie <- pie[,1:2]
Glu.pie <- Glu.pie[Glu.pie$Glutamatergic_neuron_distal!=0|Glu.pie$Glutamatergic_neuron_proximal!=0,] # 1809
Glu.pie <- Glu.pie[Glu.pie$Glutamatergic_neuron_distal==0|Glu.pie$Glutamatergic_neuron_proximal==0,] # 1014
table(Glu.pie$Glutamatergic_neuron_distal) # 785
table(Glu.pie$Glutamatergic_neuron_proximal) # 229

Gaba.pie <- pie[,3:4]
Gaba.pie <- Gaba.pie[Gaba.pie$GABAergic_interneuron_distal!=0|Gaba.pie$GABAergic_interneuron_proximal!=0,] # 1949
Gaba.pie <- Gaba.pie[Gaba.pie$GABAergic_interneuron_distal==0|Gaba.pie$GABAergic_interneuron_proximal==0,] # 801
table(Gaba.pie$GABAergic_interneuron_distal) # 227
table(Gaba.pie$GABAergic_interneuron_proximal) # 574


Ast.pie <- pie[,5:6]
Ast.pie <- Ast.pie[Ast.pie$Astrocyte_distal!=0|Ast.pie$Astrocyte_proximal!=0,] # 1821
Ast.pie <- Ast.pie[Ast.pie$Astrocyte_distal==0|Ast.pie$Astrocyte_proximal==0,] # 886
table(Ast.pie$Astrocyte_distal) # 578
table(Ast.pie$Astrocyte_proximal) # 308

Olig.pie <- pie[,7:8]
Olig.pie <- Olig.pie[Olig.pie$Oligodendrocyte_distal!=0|Olig.pie$Oligodendrocyte_proximal!=0,] # 1942
Olig.pie <- Olig.pie[Olig.pie$Oligodendrocyte_distal==0|Olig.pie$Oligodendrocyte_proximal==0,] # 826
table(Olig.pie$Oligodendrocyte_distal) # 534
table(Olig.pie$Oligodendrocyte_proximal) # 292


OPC.pie <- pie[,9:10]
OPC.pie <- OPC.pie[OPC.pie$OPC_distal!=0|OPC.pie$OPC_proximal!=0,] # 326
OPC.pie <- OPC.pie[OPC.pie$OPC_distal==0|OPC.pie$OPC_proximal==0,] # 250
table(OPC.pie$OPC_distal) # 95
table(OPC.pie$OPC_proximal) # 155

Mic.pie <- pie[,11:12]
table(Mic.pie$Microglia_distal) # 118

Endo.pie <- pie[,13:14]
table(Endo.pie$Endothelial_distal) # 131