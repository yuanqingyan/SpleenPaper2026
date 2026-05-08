#############################################################################################################
############                                                                                       ##########
############                      Blood CM vs Day7 after splenectomy                               ##########
############                                                                                       ##########
#############################################################################################################

inputFolder<-"./SpBlood/fastq/RawCount"
outputFolder<-"./DEG"
SaveFolder<-"./SpBlood/RDS"
source("~/code/sourceRcode/bulkRNAseq_sourceCode.R")
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
list.files(inputFolder)
dat<-read.delim(sprintf("%s/AllGeneExpCount.txt",inputFolder),header=T,sep="\t",row.names=1,stringsAsFactors = F)
colnames(dat)<-sapply(strsplit(colnames(dat),split="\\."),function(x) x[8])
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
anoFile<-read.delim("~/code/sourceRcode/sourceData/mart_export_grcm38.p6.txt",header=T,sep="\t",stringsAsFactors=FALSE);head(anoFile)
anoFile$ID<-paste(anoFile[,1],anoFile[,3],sep="_")
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]#

dat_coding<-dat[sapply(strsplit(row.names(dat),split="\\_"),function(x) x[1]) %in% codingGene$Gene.stable.ID,];dim(dat_coding)
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
pdata<-data.frame(ID=colnames(dat_coding),
                  trt=rep(c('Bld','D1','D7'), c(7,5,6)));pdata
pdata<-pdata[pdata$trt %in% c("Bld",'D7'),]
pdata$trt <- factor(pdata$trt, levels=c("Bld",'D7'))
pdata$col<-c("green",'red')[factor(pdata$trt)]
saveRDS(pdata, file=sprintf("%s/pataUsed.RDS",inputFolder))

#############################################################################################################
############                                                                                       ##########
#############################################################################################################

dat_input<-dat_coding[,paste(pdata$ID)];dim(dat_input) #
dat_input<-dat_input[rowSums(dat_input)>0,];dim(dat_input)#
sum(colnames(dat_input)==pdata$ID)
boxplot(log1p(dat_input),las=2)
##############
library(edgeR)
edgeInput1<-DGEList(dat_input,group=as.factor(pdata$trt))
edgeInput2 <- edgeInput1[rowSums(1e+06*edgeInput1$counts/expandAsMatrix(edgeInput1$samples$lib.size,dim(edgeInput1))>=1)>=ncol(edgeInput1)*0.1, ]
edgeInput2 <- calcNormFactors(edgeInput2)
edgeInput2$samples

#############################################################################################################
############                                                                                       ##########
#############################################################################################################
design.mat <- model.matrix(~pdata$trt)
colnames(design.mat) <- sapply(strsplit(colnames(design.mat),split="trt"),function(x) x[2])
colnames(design.mat)[c(1)]<- c("(Intercept)")
rownames(design.mat) <- colnames(edgeInput2)
design.mat

#############################################################################################################
############                                                                                       ##########
#############################################################################################################
cpm_data<-cpm(edgeInput2,prior.count=1,log=TRUE)

plotMDS(cpm_data, 
        col=pdata$col,
        pch=19,
        cex=2,
        xlab="Dim1",
        ylab="Dim2")


cpmDat<-cpm_data
colnames(cpmDat)<-paste(rep(c("Bld",'Splx_D7'),c(7,6)), c(1:7, 1:6), sep="_")

library(dendextend)
hc_list<-lapply(list(dat=cpmDat),function(x){
  scaldata=scale(x,center=TRUE,scale=FALSE)
  d<-dist(as.matrix(t(scaldata)))
  hc<-hclust(d,method="ward.D2")
  return(hc)
})
plot(hc_list[[1]],main="",xlab="",ylab="")

typeN="Clustering_SpD7vsBlood"
pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=8,height=6)
  plot(hc_list[[1]],main="",xlab="",ylab="");dev.off()
  
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
edgeInput3 <- estimateDisp(edgeInput2, design.mat, robust=TRUE) 
edgeInput3$common.dispersion
plotBCV(edgeInput3)
fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE)
plotQLDisp(fit)

#############################################################################################################
############                                                                                       ##########
#############################################################################################################
D7vsBld<-glmQLFTest(fit, coef=2)

DE.ALL<-list(D7vsBld=D7vsBld)

De.table<-lapply(DE.ALL,function(xy) {
  temp<-xy$table
  temp$FDR<-p.adjust(temp$PValue,method="BH")
  temp$stat_gsea<-zscoreT(sign(temp$logFC)*sqrt(temp$F),df=xy$df.total);return(temp)})

lapply(De.table,function(xy) min(xy$FDR))

De.table[[1]]["ENSMUSG00000027398_Il1b",]

pdf(sprintf("%s/MAplot_SplxD7_Blood.pdf",outputFolder),width=8,height=6)
par(mfrow=c(1,1))
compName<-names(De.table)[1]
temp<-De.table[[1]]

temp$Gene<-sapply(strsplit(row.names(temp),split="\\_"),function(ix) ix[2])
temp$col<-ifelse(temp$FDR<0.05, 'red','black')

plot(temp$logCPM,temp$logFC,xlab="LogCPM",ylab="logFC",col=temp$col,main='SplxVsBlood')
abline(h=0,col="grey",lty=2)

dev.off()



datIL1b<-pdata
datIL1b$IL1B<-as.numeric(cpmDat["ENSMUSG00000027398_Il1b",])
  
boxplot(IL1B~trt, data=datIL1b)


library(ggplot2)

pIl1b<-ggplot(datIL1b, aes(x = factor(trt), y = IL1B, fill = trt)) +
  geom_boxplot(alpha = 0.80) +
  geom_point(aes(fill = trt), size = 5, shape = 21, position = position_jitterdodge()) +
  ggtitle("Il1b")+
  scale_fill_manual(values = c("green", "purple"))+
  theme_bw()+
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())


pdf(sprintf("%s/Fig3_SpD7_Blood.pdf",outputFolder),width=6,height=6)
print(pIl1b)
dev.off()


########################################################################################
#             GSEA analysis
########################################################################################
library(msigdbr)
library(fgsea)
library(data.table)

source("/home/yyw9094/code/sourceRcode/sourceFgsea.R")
##################################################################
unique(as.data.frame(msigdbr(species = "Mus musculus", category ="C5"))$gs_subcat)
unique(as.data.frame(msigdbr(species = "Mus musculus", category ="H"))$gs_subcat)
head(as.data.frame(msigdbr(species = "Mus musculus", category ="H")))
##################################################################
H_list<-lapply(De.table,function(x){
  H_gsea<-gsea_Ana_edgeR(species = "Mus musculus",
                          category = "H",
                          subcategory = NULL,
                          DEinput=x,
                          cut_padj=0.05)
  return(H_gsea)
})


lapply(H_list,function(x) {sum(x$eaRes$padj<0.05)}) #

H_list$D7vsBld$eaRes
H_list$D7vsBld$eaRes[H_list$D7vsBld$eaRes$pathway=="HALLMARK_INFLAMMATORY_RESPONSE",]

library(ggplot2)
pdf(sprintf("%s/D7vsBlood_SpTxD7vsBld.pdf",outputFolder), width=6,height=4)
p1<-plotEnrichment(H_list$D7vsBld$m_list[['HALLMARK_INFLAMMATORY_RESPONSE']], H_list$D7vsBld$stats) + labs(title='Inflammatory Response')
print(p1)
dev.off()



