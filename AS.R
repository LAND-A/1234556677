ifyou need any information,plea
se connect me ,I'm sorry, I only p
rovided a portion of the code
iry(limma)
lry(ggplot2)
inputFile="geneMatrix.txt"  
conFile=nxujhdeuw
logFCfilter=nxnuwheuwbdaheader=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(nxsuhdeuwibcdiw
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5hxuhwwubduwygew(qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=nobxuhwyeudnuwenArrays(rt)
ple1=rd.table(conFile, header=F, sep="\t", ceck.nmes=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sahxuwhdwuewy=gsub("^ | $", "", as.vector(mple1[,1]))
sampleNahsuwqgdbq"^ | $", "", as.vector(sample2[,1]))
conDatdata[,sampleName1]
trnxswuebcdsbxcdsdata[,sampleName2]
data=cnd(conData,treatData)
conNum=ncol(conData)
treatm=ncol(treatData)
Type=c(rep("con",conNum), rep("treat",treatNum))
design hsuwiubcxzs("con","treat")
fit nxhbdsjbcontrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fhskawubdxbh
nxsdnncxhdfd
bzhasbwndsncxv
nxawjrisvcnkxbc= element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()
pdf(file=hsuwunckshe=5.5, height=5)
print(p)
dev.off()
libraoBandIdeogram
SHWUDNJSBA
nznfnsjcncmsecos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()
inputFile="normalize.txt" 
source("CIBERSORT.R")
outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=T)
outTa
setwd(workDir)
dashanbahwbbXZs.matrix(data)
group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="Treat"]
maxK=9
lot="png")
cnxhsuwhdjsbcdjs=2
)
crgCluCol=bHSRUZNDEJCNSSl)=levels(factor(Type$Cluster))
ann_colors[["CSnsjsfnuiewr width=7, height=4.5)
phenxhdsujhfnkczbcx
data=melt(rt, id.vars=c("Cluster"))
colnanxjsdhwuedsznxjbzxjshvadlab="Gene expression",
            legend.SHENZDNKCB ScrgCluCol,
            width=0.8,
            add="point")
p=p+rotate_x_text(60)
p1=p+stat_comparnxjshudweihcz  symnuSHUWDUEBXNints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
pdf(file="boxplot.pdf"SHIWHNUSDEBght=5)
ry(limma)
library(ggdgehcvdrgettttttrhfhdhkUBXSDr.txt" 
rt=read.table(ete45=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Cluster=as.vector(rt[,ncol(rt)])
data.pca=prc)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcat[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCamecrtggregate(PCAdger4ist(Cluster=PCA$Cluster), mean)
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
crgCluCol=bind(cos(theta), sin(theta))
  t(center +.frame()
for(g in levels(factor(PCA$Clusteryrrry5df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],
                                                   veganCovEllipse(covnxjsuheduwbxjabhwymean(PC1),mean(PC2))))), Cluster=g))
}
pdf(file="PCA.pdf", width=6.5, height=5)
ggale_colour_manual(name="Cluyryrues =crgCluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  geoll, aes(x=Psfsgdh2, colour=ClushftrPCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dry(limma)
librarle="cluster.txt"
immFile=64657-Results.txt"
immud.table4757ile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=i475.names(data), row.names(Cluster))
rt=cbind(da576ta[sameSample,,drop=F], Cluster[sameSample,"Cluster",drop=F])
rt=rt[order(rt$Cluster, decreasing=F),]
conNum=nrow(rt[rt$Cluster=="C1",])
treatNumpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/ryyyr",cex=2)
re('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()
data=rt
data=melt(data, id.vars=yth(group)]
boxplot=ggboxplotuster),symnum.argeyr 0.001, 0.01, 0.05, 1)ryry**", "**", "*", "")), label="p.signif")
pdf(firete.diff.pdf", width=8, height=6)
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
rowna53rt)=rt[,1]
exp=rt[etl(rt)]
dimnameeeeeeeeeeete(rownames(exp),colnames(exp))
data=matrix(axgeet5(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=sgfdddddddts=getGmt(gmtaliretn(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalixeC2,drop=F]
conNum=ncol(dataC1)
treatNum=ncol(dataC2)
data=cbind(dataC1, dataC2)
Type=c(rep("C1",conNum), rep("C2",treatNum))
outTab=6sfea.frame()
for(i in row.names(data)){
  test=t.test(data[i,] ~ Type)
  pvalue=t666786ue
  t=test$statistic
  if(pvalue<0.05){
    Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
    outTab=rbind(outTabttu76y=i, t, pvalue, Sig))
  }
}
termNum=10
outTab=outTab[order(outTab$t),]
outTab=outTab[c(1:termNum,(nrow(outTab)-termNum):nrow(outTab)),]
pdf(file="barplot.pdf", width=9, height=6)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
              palette=c("bluavftthufjrttuu3", "val = "asc", sort.by.groups = T,
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
rt$Type=ifelse(g"con", 0, 1)
x=as.matrix(rt[,1:(ncol(rt)-1)])
y=rt[,"Type"]
fit=glmnet(x, y, family = "binomial", alpha=1)
pdf(file="lasso.pth=6,height=5.5)
plot(fit)
dev.off()
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lae[-1]
write.table(lassoGene, file="lasso.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)
outTab=rt[,c(lassoGened=row.names(outTab), outTab)
write.table(outTab, file="lasso.geneExp.txt", sep="\t", quote=F, row.names=F)
input.names(summ$coefficients)[summ$coefficients[,"Pr(>|z|)"]<0.05]
rt2=rt[,c("Type", newGene)]
fit2=glm(Type ~ ., family="binomial", data=rt2)
fit2=step(fit2)
conf=confint(fit2, level=0.95)
summ2=summary(fit2)
gene=rw.names(summoefficients[,"Estimate"]
OR=exp(summ2$coefficients[,"Estimate"])
OR.9[,1])
OR.95H=exp(conf[,2])
pvalue=summ2$coefficients[,"Pr(>|z|)"]
genf=cbind(gene,coef, OR, f[-1,], file="modelGene.txt", sep="\t", quote=F, row.names=F)
pred, type="response")
bind(rt[,c("Type", gene[-1])], riskScore=pred)
outTa(id=row.names(outTab), outTab)
writab, file="risk.txt", sep="\t", quote=F, row.names=F)
librmma)
library(pheatmap)
ize.txt"  
geneFile="modelGene.txt" 
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=ax(rt)
rownames(rt)=rt[,1]
exp:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
Type=gsub("(.*?)\\_(.*)", "\\2", colnames(data))
colnames(data)=gsub("(.*?)\\_
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=8, height=5)
pheatmap(data, 
         annotation=Type, 
         colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
         cluster_cols=F,
         show_colnames=F,
      "row",
         fontsize = 20,
     ze_row=7,
         fontsize=7)
dev.off()
library(reshape2)
libubr)
inputFile="CIBERSORT-Results.txt"
rt=reacon=grepl("_con", rownames(rt), ignore.case=T)
treat=grepl("_treat", rownames(rt), ignore.case=T)
c=rt[con,]
treatData=rt[treat,]
cm=nrow(conData)
treatNum=nrow(treatData)
dat(conData,treatData))
pdf(file="barplot.pdf", width=15, height=8)
cainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barploa2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
t=0,xpd=T)
rect(xletext(a1[conNum]/2,-0.035,"Control",cex=1.8)
reca1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="#EE0000FF")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"Treat",cex=1.8)
ytick2 = cumsuta)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend'usr')[2]*0.98usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(a.frame(t(data)), Type)
data=melt(data,(data)=c("Type", "Immune", "Expression")
group=levels(factor(data$Type))
bioCol="#EE0000FF","#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=biolength(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
                  xlab="",
                  ylab=egend.title="Type",
                  notch=T,
                  #add="point",
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
pdf(file="immun", width=8, height=6)
print(boxplot)
dev.off()
library(limma)
library(reshape2)
library(se)
library(ggplot2)
expFile="normalize.txt"
geneFile="modelGene.txt"
immFile="CIBERSORT-Results.txt"
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(r
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(aereps(data)
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=daas.vector(geneRT[,1]),]
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
data=t(data)
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=indataple,,drop=F]
immune=immune[sameSample,,drop=F]
outTarame()
for colnames(immune)){
  if(sd(immune[,cell])==0){next}
  for(gene mes(dats.numeric(immune[,cell])
    y=as.numeric(data[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
  }
}
utTab$cor=as.numeric(outTab$cor)
pf(file="cor.pdf", width=7, height=4.8)
gplot(outTab, aes(Immune, Gene)) + 
om_tile(aes(fill = cor), colour = "grey", size = 1)+
  s_fill_gradi2(low = "#", mid = " high D") + 
  geomeabel=text),col ="bsize = 3) +
 me_minimal() + 
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        at.text.y = element_text(size = 8, face = "bold")) +    
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) + 
  scale_x_discrete(position = "bottom")     
dev.off()
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(
gene="ADCY7"
expFile="normalize.txt" 
gmtFile=""
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
group=gsub(",group=="Treat",drop=F]
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(xeaScore=normalize(ssgseaScore)
ssgseaScore=ssgseaScore[order(apply(ssgseaScore,1,sd),decreasing=T),]
ssgseaScore=ssgseaScore[1:50,]
lowName=ghName=colnames(data)[data[gene,]>=median(data[gene,])]
lowScore=ssgseae[,highName]
data=cbind(lowScore, highScore)
conNum=ncol(lowScore)
treatNum=ncol(hcore)
Type=c(rep("Control",conNum), rep("Treat",treatNum))
outTab=data.frame()
for(i in row.names(data)){
  test=t.test(da=test$p.value
  t=test$statistic
  Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
  outTabnd(outTab, cbind(Pathway=i, t, pvalue, Sig))
}
pdf(file="barplot.pdf", width=12, height=9)
outTab$t=americ(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
              palen3","grey","red3"), sort.val = "asc", sort.by.groups = T,
              rotate=TRUE, legend="right", title=gene,
              xlrm", ylab="t value of GSVA score",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.ff()
library(limma)
library(NMF)
libraggplot2)
library(ggalluvial)
library
library(CellChat)
expFile="" 
annFile=""
rttable(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=datannFile, header=T, sep="\t", check.names=F, row.names=1)
qishit <- createCellChat(object = data, meta = meta, group.by = "labels")
qishit <- setIdent(cellchat, ident.use="labels")
grou CellChatDB.human
tDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")


