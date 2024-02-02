###S7
#调取3个基因子网络
Line3 <- read.csv("~/F5/Line.csv")
CALM1 <- Line3[Line3$Gene=="CALM1",] # 207
table(CALM1$describe) # disease6 miRNA201
write.csv(CALM1, file = "CALM1_line.csv")
CALM1.Node <- CALM1[,c(2,3)] # 207
CALM1.Node[nrow(CALM1.Node) + 1,] = c("CALM1","gene")
write.csv(CALM1.Node, file = "CALM1_Node.csv")
#
APP <- Line3[Line3$Gene=="APP",] # 205
table(APP$describe) # disease13  miRNA 192
write.csv(APP, file = "APP_line.csv")
APP.Node <- APP[,c(2,3)] # 205
APP.Node[nrow(APP.Node) + 1,] = c("APP","gene")
write.csv(APP.Node, file = "APP_Node.csv")
#
ACTB <- Line3[Line3$Gene=="ACTB",] # 217
table(ACTB$describe) # disease 10  miRNA 207
write.csv(ACTB, file = "ACTB_line.csv")
ACTB.Node <- ACTB[,c(2,3)] # 217
ACTB.Node[nrow(ACTB.Node) + 1,] = c("ACTB","gene")
write.csv(ACTB.Node, file = "ACTB_Node.csv")
################################################################################