gc()
setwd("~/20221004_scAPAtrap/F5miRNAs")
###F5A
#根据四象图(diff_usage)做疾病数据库pie矩阵
#读取细胞亚型疾病风险基因
Glu.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Glutamatergic_neuron.Disease.csv")
Gaba.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/GABAergic_interneuron.Disease.csv")
Ast.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Astrocyte.Disease.csv")
Olig.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Oligodendrocyte.Disease.csv")
OPC.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/OPC.Disease.csv")
Mic.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Microglia.Disease.csv")
End.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Endothelial.Disease.csv")

Glu.dis$celltype <- "Glutamatergic_neuron"
Gaba.dis$celltype <- "GABAergic_interneuron"
Ast.dis$celltype <- "Astrocyte"
Olig.dis$celltype <- "Oligodendrocyte"
OPC.dis$celltype <- "OPC"
Mic.dis$celltype <- "Microglia"
End.dis$celltype <- "Endothelial"

All.dis <- rbind(rbind(rbind(rbind(rbind(rbind(Glu.dis,Gaba.dis),
                                         Ast.dis),Olig.dis),OPC.dis),Mic.dis),End.dis)
All.dis1 <- All.dis[,c(2,3,5,6)]
table(All.dis1$Disease,All.dis1$celltype)
All.dis1$Disease[All.dis1$Disease=="Alzheimer.s.disease"] <- "Alzheimer's Disease"
All.dis1$Disease[All.dis1$Disease=="Alzheimer's Disease"] <- "AD"
All.dis1$Disease[All.dis1$Disease=="Bipolar.disorder"] <- "Bipolar Disorder"
All.dis1$Disease[All.dis1$Disease=="Parkinson.s.disease"] <- "Parkinson Disease"
All.dis1$Disease[All.dis1$Disease=="Parkinson Disease"] <- "PD"
All.dis1$Disease[All.dis1$Disease=="Amyotrophic Lateral Sclerosis"] <- "ALS"
All.dis1$Disease[All.dis1$Disease=="Autism"] <- "Autism Spectrum Disorders"
All.dis1$Disease[All.dis1$Disease=="Depression"] <- "Depressive disorder"
All.dis1$Disease[All.dis1$Disease=="Attention deficit hyperactivity disorder"] <- "ADHD"
All.dis1$Disease[All.dis1$Disease=="Huntington Disease"] <- "HD"
All.dis1$Disease[All.dis1$Disease=="Multiple.sclerosis"] <- "MS"
All.dis1 <- All.dis1[All.dis1$Disease!='Psychosis',]
All.dis1 <- All.dis1[All.dis1$Disease!='Cognitive',]
All.dis1 <- All.dis1[All.dis1$Disease!='Dementia',]
All.dis1 <- All.dis1[All.dis1$Disease!='Insomnia',]
All.dis1 <- All.dis1[All.dis1$Disease!='Intelligence',]
All.dis1 <- All.dis1[All.dis1$Disease!='PTSD',]
All.dis1 <- All.dis1[All.dis1$Disease!='Risk.behavior',]
All.dis1 <- All.dis1[All.dis1$Disease!='Neuroticism',]
All.dis1 <- All.dis1[All.dis1$Disease!='MSA',]

table(All.dis1$Disease,All.dis1$celltype) # 1091

#去重
All.dis1 <- All.dis1[!duplicated(All.dis1),] # 1036
##write.csv(All.dis1,'All.dis1.csv')

#构建pie矩阵
pie <- data.frame(matrix(NA,105,12)) # 15*7

disease <- unique(All.dis1$Disease)
celltype <- unique(All.dis1$celltype)

colnames(pie) <- c("disease","celltype","group","r","s1","s2","s3","s4","s1.gene","s2.gene","s3.gene","s4.gene")
pie$disease <- rep(disease,each = 7)
pie$celltype <- rep(celltype,times = 15)
pie$group <- c(1:105)

#统计每种疾病里出现的四象限基因
for (j in 1:7){
  for (i in 1:15) {
    for (n in 1:4) {
      A <- All.dis1[All.dis1$celltype==celltype[j],]
      B <- A[A$Disease==disease[i],]
      C <- B[B$flag==n,]
      pie[(7*i-7+j),4+n] <- length(C$Gene)
      pie[(7*i-7+j),8+n] <- paste(as.character(C$Gene), collapse=",")
    }
  }
}

pie$r <- rowSums(pie[5:8]) # pie半径总count
##write.csv(pie,'pie.csv')

pie1 <- pie[,1:8]
#去除总数为0
pie2 <- pie1[pie1$r!=0,]
#缩放至5-10 (1-49)
k <- (10-5)/(49-1)  
pie2$r <- 5+k*(pie2$r-1)  

#pie矩阵作图
library(ggplot2)
library(jjPlot) #devtools::install_github('junjunlab/jjPlot')

View(pie1)
df.long <- reshape2::melt(pie2,id.vars = c("disease","celltype","group",'r'),variable.name = 'type',value.name = 'per')

df.long$disease <- factor(df.long$disease, levels = disease)
df.long$celltype <- factor(df.long$celltype, levels = c("Glutamatergic_neuron","GABAergic_interneuron","Astrocyte",
                                                        "Oligodendrocyte","OPC","Microglia","Endothelial"))
df.long[order(df.long$disease,df.long$celltype),]

ggplot(df.long,aes(x =disease,y =celltype ,group = group)) +
  geom_jjPointPie(aes(pievar = per,
                      fill = type,
                      filltype = type,
                      width = r/7),
                  color = 'white',line.size = 0.5) +
  coord_fixed(ratio = 0.8)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = c(s1="#EB5A40",s2="#6AADE0",s3="#4A84B5",s4="#870E07"))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,size = 10),
        axis.text.y = element_text(size = 10),
        text=element_text(family="sans"))+
  scale_x_discrete(labels=function(x) str_wrap(x, width=40)) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="dotted"))+
  geom_vline(xintercept = c(seq(from=1.5, to=14.5, by=1)),
             color="lightgrey", linewidth=0.7, linetype="dotted")

###F5B
library(multiMiR)
db.tables <- multimir_dbTables()                 
db.tables
db.info <- multimir_dbInfo()                      
db.info
predicted_tables()           
diseasedrug_tables()         
reverse_table_lookup("targetscan")                
db.count <- multimir_dbCount()                    
db.count
apply(db.count[,-1], 2, sum)                     
## 查看multiMiR收录的miRNA、gene、drug和disease
# 当前版本multiMiR收录来自人类、小鼠和大鼠的 5830 个 miRNA 和 97186 个靶基因，以及 64 种药物和 223 个疾病术语
miRNAs   <- list_multimir("mirna", limit <- 10)
genes    <- list_multimir("gene", limit <- 10)
drugs    <- list_multimir("drug", limit <- 10)
diseases <- list_multimir("disease", limit <- 10)
head(miRNAs)
head(genes)
head(drugs)
head(diseases)

#合并每个细胞类型四象限图中出现的基因
Glutamatergic_neuron <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Glutamatergic neuron.csv")
GABAergic_interneuron <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/GABAergic interneuron.csv")
Astrocyte <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Astrocyte.csv")
Oligodendrocyte <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Oligodendrocyte.csv")
OPC <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/OPC.csv")
Microglia <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Microglia.csv")
Endothelial <- read.csv(file = "~/202210_movAPA/sixiangxiantu/covid19_control_cellType/Endothelial.csv")
Glutamatergic_neuron$celltype <- "Glutamatergic_neuron"
GABAergic_interneuron$celltype <- "GABAergic_interneuron"
Astrocyte$celltype <- "Astrocyte"
Oligodendrocyte$celltype <- "Oligodendrocyte"
OPC$celltype <- "OPC"
Microglia$celltype <- "Microglia"
Endothelial$celltype <- "Endothelial"
All <- rbind(rbind(rbind(rbind(rbind(rbind(Glutamatergic_neuron,GABAergic_interneuron),
                                     Astrocyte),Oligodendrocyte),OPC),Microglia),Endothelial)
table(All$celltype)
all.gene <- unique(All$hgnc_symbol)
# 已知mRNA列表检索miRNA
example4_2 <- get_multimir(org = "hsa",
                           target = all.gene,
                           table = "validated",          #validated实验验证的三个数据库
                           summary = TRUE,
                           predicted.cutoff.type = "p",
                           predicted.cutoff = 10,
                           use.tibble = TRUE)

save(example4_2, file = "target_miRNAs_validated.Rdata")      

table(example4_2@data$type)                    
example4_result2 <- example4_2@data
head(example4_result2)
example4_sum2 <- example4_2@summary
head(example4_sum2)
apply(example4_2@summary[, 6:8], 2, sum)
a <- unique(example4_2@summary$mature_mirna_id)  # 33285条对应，2404个miRNAs

#构建cytoscape表格
Line <- example4_2@summary[,c(3,2)] # 33285
table(!duplicated(Line))
Line <- Line[!duplicated(Line),] # 32917
colnames(Line) <- c("Gene","Target")
Line$describe <- "miRNA"

#构建cytoscape中基因与疾病的关系
Glu.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Glutamatergic_neuron.Disease.csv")
Gaba.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/GABAergic_interneuron.Disease.csv")
Ast.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Astrocyte.Disease.csv")
Olig.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Oligodendrocyte.Disease.csv")
OPC.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/OPC.Disease.csv")
Mic.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Microglia.Disease.csv")
End.dis <- read.csv(file = "~/20221004_scAPAtrap/F4/Endothelial.Disease.csv")
All.dis <- rbind(rbind(rbind(rbind(rbind(rbind(Glu.dis,Gaba.dis),
                                         Ast.dis),Olig.dis),OPC.dis),Mic.dis),End.dis)
All.dis1 <- All.dis[,c(2,3,4)]
table(All.dis1$Disease,All.dis1$DB)
All.dis1$Disease[All.dis1$Disease=="Alzheimer.s.disease"] <- "Alzheimer's Disease"
All.dis1$Disease[All.dis1$Disease=="Alzheimer's Disease"] <- "AD"
All.dis1$Disease[All.dis1$Disease=="Bipolar.disorder"] <- "Bipolar Disorder"
All.dis1$Disease[All.dis1$Disease=="Parkinson.s.disease"] <- "Parkinson Disease"
All.dis1$Disease[All.dis1$Disease=="Parkinson Disease"] <- "PD"
All.dis1$Disease[All.dis1$Disease=="Amyotrophic Lateral Sclerosis"] <- "ALS"
All.dis1$Disease[All.dis1$Disease=="Autism"] <- "Autism Spectrum Disorders"
All.dis1$Disease[All.dis1$Disease=="Depression"] <- "Depressive disorder"
All.dis1$Disease[All.dis1$Disease=="Attention deficit hyperactivity disorder"] <- "ADHD"
All.dis1$Disease[All.dis1$Disease=="Huntington Disease"] <- "HD"
All.dis1$Disease[All.dis1$Disease=="Multiple.sclerosis"] <- "MS"
All.dis1 <- All.dis1[All.dis1$Disease!='Psychosis',]
All.dis1 <- All.dis1[All.dis1$Disease!='Cognitive',]
All.dis1 <- All.dis1[All.dis1$Disease!='Dementia',]
All.dis1 <- All.dis1[All.dis1$Disease!='Insomnia',]
All.dis1 <- All.dis1[All.dis1$Disease!='Intelligence',]
All.dis1 <- All.dis1[All.dis1$Disease!='PTSD',]
All.dis1 <- All.dis1[All.dis1$Disease!='Risk.behavior',]
All.dis1 <- All.dis1[All.dis1$Disease!='Neuroticism',]
All.dis1 <- All.dis1[All.dis1$Disease!='MSA',]

table(All.dis1$Disease,All.dis1$DB)
Line2 <- All.dis1[,c(1,2)] # 1091
Line2 <- Line2[!duplicated(Line2),] # 853
colnames(Line2) <- c("Gene","Target")
Line2$describe <- "disease"
A <- unique(Line2$Gene) #267
#Line
Line <- data.frame(Line)
Line3 <- rbind(Line,Line2) # 33770
write.csv(Line3, file = "Line.csv")

#node
Node1 <- data.frame(unique(Line3$Gene)) # 409
colnames(Node1) <- "Gene"
Node1$describe <- "DEG"
Node2 <- Line3[,c(2,3)]
Node2 <- Node2[!duplicated(Node2),] # 2419
colnames(Node2) <- c("Gene","describe")
Node <- rbind(Node1,Node2) # 2828
write.csv(Node, file = "Node.csv")
################################################################################


