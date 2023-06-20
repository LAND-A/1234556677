library("RCircos")
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t", check.names=F)
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.7
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)
pdf(file="RCircos.pdf", width=7, height=7)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()
inputFile="normalize.txt" 
source("CIBERSORT.R")
outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=T)
outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)
library(ConsensusClusterPlus)
expFile="GeneExp.txt" 
workDir=""
setwd(workDir)
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)
group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="Treat"]
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")
calcICL(results, title="consensusScore", plot="png")
clusterNum=2
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cluster")
cluster$Cluster=paste0("C", cluster$Cluster)
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="cluster.txt", sep="\t", quote=F, col.names=F)
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
clusterFile="cluster.txt"
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[order(rt$Cluster),]
data=t(rt[,1:(ncol(rt)-1),drop=F])
Type=rt[,ncol(rt),drop=F]
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
crgCluCol=bioCol[1:length(levels(factor(Type$Cluster)))]
names(crgCluCol)=levels(factor(Type$Cluster))
ann_colors[["Cluster"]]=crgCluCol
pdf("heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()
data=melt(rt, id.vars=c("Cluster"))
colnames(data)=c("Cluster", "Gene", "Expression")
p=ggboxplot(data, x="Gene", y="Expression", color = "Cluster",
            xlab="",
            ylab="Gene expression",
            legend.title="Cluster",
            palette = crgCluCol,
            width=0.8,
            add="point")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Cluster),
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()
library(limma)
library(ggplot2)
clusterFile="cluster.txt" 
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Cluster=as.vector(rt[,ncol(rt)])
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
crgCluCol=bioCol[1:length(levels(factor(Cluster)))]
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Cluster))){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],
                                                   veganCovEllipse(cov.wt(cbind(PC1,PC2),
                                                                          wt=rep(1/length(PC1),length(PC1)))$cov,
                                                                   center=c(mean(PC1),mean(PC2))))), Cluster=g))
}
pdf(file="PCA.pdf", width=6.5, height=5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
  scale_colour_manual(name="Cluster", values =crgCluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=Cluster), size=1, linetype=2)+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
library(limma)
library(reshape2)
library(ggpubr)
clusterFile="cluster.txt"
immFile="CIBERSORT-Results.txt"
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
group=gsub("(.*)\\_(.*)", "\\2", row.names(immune))
data=immune[group=="treat",,drop=F]
Cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(Cluster))
rt=cbind(data[sameSample,,drop=F], Cluster[sameSample,"Cluster",drop=F])
rt=rt[order(rt$Cluster, decreasing=F),]
conNum=nrow(rt[rt$Cluster=="C1",])
treatNum=nrow(rt[rt$Cluster=="C2",])
data=t(rt[,-ncol(rt)])
pdf(file="barplot.pdf", width=14.5, height=8)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data, col=col, xaxt="n", yaxt="n", ylab="Relative Percent", cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"C1",cex=2)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5 , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"C2",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()
data=rt
data=melt(data, id.vars=c("Cluster"))
colnames(data)=c("Cluster", "Immune", "Expression")
group=levels(factor(data$Cluster))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", color="Cluster",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Cluster",
                  add="point",
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Cluster),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
expFile="normalize.txt"   
clusterFile="cluster.txt" 
gmtFile=""
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
ssgseaScore=gsva(data, geneSets, method='gsva')
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
nameC1=row.names(cluster[cluster$Cluster=="C1",,drop=F])
nameC2=row.names(cluster[cluster$Cluster=="C2",,drop=F])
dataC1=ssgseaScore[,nameC1,drop=F]
dataC2=ssgseaScore[,nameC2,drop=F]
conNum=ncol(dataC1)
treatNum=ncol(dataC2)
data=cbind(dataC1, dataC2)
Type=c(rep("C1",conNum), rep("C2",treatNum))
outTab=data.frame()
for(i in row.names(data)){
  test=t.test(data[i,] ~ Type)
  pvalue=test$p.value
  t=test$statistic
  if(pvalue<0.05){
    Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
    outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
  }
}
termNum=10
outTab=outTab[order(outTab$t),]
outTab=outTab[c(1:termNum,(nrow(outTab)-termNum):nrow(outTab)),]
pdf(file="barplot.pdf", width=9, height=6)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
              palette=c("blue3", "red3"), sort.val = "asc", sort.by.groups = T,
              rotate=TRUE, legend="right", title="",
              xlab="Term", ylab="t value of GSVA score, C2 vs C1",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()
library(glmnet)              
inputFile="normalize.txt"    
geneFile="Gene.txt"    
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=as.data.frame(data)
rt$Type=ifelse(group=="con", 0, 1)
x=as.matrix(rt[,1:(ncol(rt)-1)])
y=rt[,"Type"]
fit=glmnet(x, y, family = "binomial", alpha=1)
pdf(file="lasso.pdf",width=6,height=5.5)
plot(fit)
dev.off()
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="lasso.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)
outTab=rt[,c(lassoGene, "Type")]
outTab=cbind(id=row.names(outTab), outTab)
write.table(outTab, file="lasso.geneExp.txt", sep="\t", quote=F, row.names=F)
inputFile="lasso.geneExp.txt"
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
fit=glm(Type ~ ., family="binomial", data=rt)
summ=summary(fit)
newGene=row.names(summ$coefficients)[summ$coefficients[,"Pr(>|z|)"]<0.05]
rt2=rt[,c("Type", newGene)]
fit2=glm(Type ~ ., family="binomial", data=rt2)
fit2=step(fit2)
conf=confint(fit2, level=0.95)
summ2=summary(fit2)
gene=row.names(summ2$coefficients)
coef=summ2$coefficients[,"Estimate"]
OR=exp(summ2$coefficients[,"Estimate"])
OR.95L=exp(conf[,1])
OR.95H=exp(conf[,2])
pvalue=summ2$coefficients[,"Pr(>|z|)"]
geneCoef=cbind(gene,coef, OR, OR.95L, OR.95H, pvalue)
write.table(geneCoef[-1,], file="modelGene.txt", sep="\t", quote=F, row.names=F)
pred=predict(fit2, type="response")
outTab=cbind(rt[,c("Type", gene[-1])], riskScore=pred)
outTab=cbind(id=row.names(outTab), outTab)
write.table(outTab, file="risk.txt", sep="\t", quote=F, row.names=F)
library(limma)
library(pheatmap)
expFile="normalize.txt"  
geneFile="modelGene.txt" 
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
Type=gsub("(.*?)\\_(.*)", "\\2", colnames(data))
colnames(data)=gsub("(.*?)\\_(.*)", "\\1", colnames(data))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=8, height=5)
pheatmap(data, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
         cluster_cols=F,
         show_colnames=F,
         scale="row",
         fontsize = 7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()
library(reshape2)
library(ggpubr)
inputFile="CIBERSORT-Results.txt"
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
con=grepl("_con", rownames(rt), ignore.case=T)
treat=grepl("_treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData,treatData))
pdf(file="barplot.pdf", width=15, height=8)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="#008B45FF")
text(a1[conNum]/2,-0.035,"Control",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="#EE0000FF")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"Treat",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
group=levels(factor(data$Type))
bioCol=c("#008B45FF","#EE0000FF","#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Type",
                  notch=T,
                  #add="point",
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)
expFile="normalize.txt"
geneFile="modelGene.txt"
immFile="CIBERSORT-Results.txt"
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
data=t(data)
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]
outTab=data.frame()
for(cell in colnames(immune)){
  if(sd(immune[,cell])==0){next}
  for(gene in colnames(data)){
    x=as.numeric(immune[,cell])
    y=as.numeric(data[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
  }
}
outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=7, height=4.8)
ggplot(outTab, aes(Immune, Gene)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold")) +    
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) + 
  scale_x_discrete(position = "bottom")     
dev.off()
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
gene="ADCY7"
expFile="normalize.txt" 
gmtFile=""
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaScore=ssgseaScore[order(apply(ssgseaScore,1,sd),decreasing=T),]
ssgseaScore=ssgseaScore[1:50,]
lowName=colnames(data)[data[gene,]<median(data[gene,])] 
highName=colnames(data)[data[gene,]>=median(data[gene,])]
lowScore=ssgseaScore[,lowName]
highScore=ssgseaScore[,highName]
data=cbind(lowScore, highScore)
conNum=ncol(lowScore)
treatNum=ncol(highScore)
Type=c(rep("Control",conNum), rep("Treat",treatNum))
outTab=data.frame()
for(i in row.names(data)){
  test=t.test(data[i,] ~ Type)
  pvalue=test$p.value
  t=test$statistic
  Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
  outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
}
pdf(file="barplot.pdf", width=12, height=9)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
              palette=c("green3","grey","red3"), sort.val = "asc", sort.by.groups = T,
              rotate=TRUE, legend="right", title=gene,
              xlab="Term", ylab="t value of GSVA score",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()
library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)
expFile="expMatirx.txt" 
annFile="cellAnn.txt"
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
meta=read.table(annFile, header=T, sep="\t", check.names=F, row.names=1)
cellchat <- createCellChat(object = data, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use="labels")
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net=subsetCommunication(cellchat)
write.table(file="Comm.network.xls", df.net, sep="\t", row.names=F, quote=F)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
pdf(file="cellNetworkCount.pdf", width=7, height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf(file="cellNetworkWeight.pdf", width=7, height=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
dev.off()
pdf(file=paste0("COMM.", pathways.show , ".hierarchy.pdf"), width=12, height=6)
hierarchy=netVisual_aggregate(cellchat, signaling=pathways.show, layout="hierarchy",  vertex.receiver=seq(1,4), vertex.size = groupSize)
print(hierarchy)
dev.off()
pdf(file=paste0("COMM.", pathways.show , ".netAnalysis.pdf"), width=6, height=5)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis=netAnalysis_signalingRole_network(cellchat, signaling =pathways.show, width = 8, height = 5, font.size = 12)
print(netAnalysis)
dev.off()

