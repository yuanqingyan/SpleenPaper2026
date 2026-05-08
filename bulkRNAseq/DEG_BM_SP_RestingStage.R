
#############################################################################################################
############                                                                                       ##########
############            This is for DEG between Bone marrow vs Spleen at resting stage             ##########
############                                                                                       ##########
#############################################################################################################

inputFolder<-"./Apr2023/RawCount"
outputFolder<-"./DEG"
SaveFolder<-"./Apr2023/RDS"
source("~/code/sourceRcode/bulkRNAseq_sourceCode.R")
#############################################################################################################
############                                                                                       ##########
############                           read the data                                               ##########
############                                                                                       ##########
#############################################################################################################
list.files(inputFolder)
dat<-read.delim(sprintf("%s/AllGeneExpCount.txt",inputFolder),header=T,sep="\t",row.names=1,stringsAsFactors = F)
colnames(dat)<-sapply(strsplit(colnames(dat),split="\\."),function(x) x[9])
#############################################################################################################
############                                                                                       ##########
############                           protein coding only                                         ##########
############                                                                                       ##########
#############################################################################################################
anoFile<-read.delim("~/code/sourceRcode/sourceData/mart_export_grcm38.p6.txt",header=T,sep="\t",stringsAsFactors=FALSE);head(anoFile)
anoFile$ID<-paste(anoFile[,1],anoFile[,3],sep="_")
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]

dat_coding<-dat[sapply(strsplit(row.names(dat),split="\\_"),function(x) x[1]) %in% codingGene$Gene.stable.ID,];dim(dat_coding)#
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
pdata<-data.frame(ID=colnames(dat_coding),
                  trt=sapply(strsplit(colnames(dat_coding),split="\\_"),function(x) substr(x[1],1,2)));pdata
pdata$trt <- factor(pdata$trt, levels=c("BM","Sp"))
pdata$col<-c("green","purple")[factor(pdata$trt)]

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
#############################################################################################################
cpm_data<-cpm(edgeInput2,prior.count=1,log=TRUE)
listCPM<-list(cpm=cpm_data,edeeInput2=edgeInput2)
saveRDS(listCPM,file=sprintf("%s/listCPM_forBMSP_Resting.rds",SaveFolder))
boxplot(cpm_data,col=pdata$col)
library(dendextend)
hc_list<-lapply(list(dat=cpm_data),function(x){
  scaldata=scale(x,center=TRUE,scale=FALSE)
  d<-dist(as.matrix(t(scaldata)))
  hc<-hclust(d,method="complete")
  return(hc)
})
plot(hc_list[[1]],main="",xlab="",ylab="")


pdf(sprintf("%s/MDS_BmvsSp_resting.pdf",outputFolder),width=6,height=6)
plotMDS(cpm(edgeInput2,prior.count=1,log=TRUE), 
        col=pdata$col,
        pch=19,
        cex=2,
        xlab="Dim1",
        ylab="Dim2")
dev.off()

#############################################################################################################
############                                                                                       ##########
#############################################################################################################
edgeInput3 <- estimateDisp(edgeInput2, design.mat, robust=TRUE) ##similar to estimateCommonDisp and estimateTagwiseDisp. 
edgeInput3$common.dispersion
plotBCV(edgeInput3)
fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE)
plotQLDisp(fit)

#############################################################################################################
############                                                                                       ##########
#############################################################################################################
SpvsBM<-glmQLFTest(fit, coef=2)

DE.ALL<-list(SpvsBM1=SpvsBM)
names(DE.ALL);lapply(DE.ALL,function(x) summary(decideTests(x)))
saveRDS(DE.ALL,file=sprintf("%s/DE.ALL_resting.RDS",SaveFolder))

DE.ALL<-readRDS(file=sprintf("%s/DE.ALL_resting.RDS",SaveFolder))
De.table<-lapply(DE.ALL,function(x) {
  temp<-x$table
  temp$FDR<-p.adjust(temp$PValue,method="BH")
  temp$Bon<-p.adjust(temp$PValue,method="bon")
  temp$stat_gsea<-zscoreT(sign(temp$logFC)*sqrt(temp$F),df=x$df.total);return(temp)})
lapply(De.table,function(x) head(x))
saveRDS(De.table,file=sprintf("%s/DeTable_resting.rds",outputFolder))

FDRCut=0.05;logCut=log2(2)
singCount<-do.call("rbind",lapply(De.table,function(x) {
  tSig<-sum(x$FDR<FDRCut & abs(x$logFC)>logCut)
  uSig<-sum(x$FDR<FDRCut & x$logFC>logCut)
  dSig<-sum(x$FDR<FDRCut & x$logFC<(-logCut))
  out<-c(tSig,uSig,dSig)}))
colnames(singCount)<-c("TotalSig","UpSig","DnSig");singCount
write.table(data.frame(Comparison=row.names(singCount),singCount),file=sprintf("%s/singCount_resting.xls",outputFolder),col.names=T,row.names=F,sep="\t")

###write stat out
lapply(1:length(De.table),function(x){
  compName<-names(De.table)[x]
  temp<-De.table[[x]]
  write.table(data.frame(Gene=row.names(temp),temp),file=sprintf("%s/StatOut_resting.%s.xls",outputFolder,compName),row.names=F,col.names=T,sep="\t")
})

#############################################################################################################
############                                                                                       ##########
#############################################################################################################
length(De.table)

upgene<-c("Top2a",'Stmn1','Retnlg','Spp1','Mmp9','Cit','Lcn2','Kif11')
dngene<-c("Ppp1r15a",'Ccl3','hist1h1c','Pim3','Zfp36','Ccrl2','Btg2','Gm37352','Rgs1','Nr4a1')

proliGene<-c("Cdk1", "Ccna2", "Ccnb1", "Ccnb2", "Wee1", "Mki67", "Top2a")
Il1bGene<-c("Il1b")

pdf(sprintf("%s/Volcano.pdf",outputFolder),width=8,height=8)
par(mfrow=c(1,1))
compName<-names(De.table)[1]
temp<-De.table[[1]]
temp$Gene<-sapply(strsplit(row.names(temp),split="\\_"),function(ix) ix[2])
plot(temp$logFC,-log10(temp$FDR),xlab="Log2 FC",ylab="-log10 FDR",main=compName)
labGene<-temp[temp$Gene %in% proliGene,]
points(labGene$logFC, -log10(labGene$FDR), pch=20,col='red')
labGene2<-temp[temp$Gene %in% Il1bGene,]
points(labGene2$logFC, -log10(labGene2$FDR), pch=20,col='green')
abline(h=-log10(0.05),col="red",lty=2)
dev.off()



  
pdf(sprintf("%s/Fig1A_MAplot.pdf",outputFolder),width=10,height=8)
par(mfrow=c(1,1))
compName<-names(De.table)[1]
temp<-De.table[[1]]

proliGene<-c("Cdk1", "Ccna2", "Ccnb1", "Ccnb2", "Wee1", "Mki67", "Top2a")
Il1bGene<-c("Il1b","Nlrp3",'Traf1', 'Ccrl2')
target_genes <- c(proliGene, Il1bGene)
temp$Gene<-sapply(strsplit(row.names(temp),split="\\_"),function(ix) ix[2])
temp$col<-ifelse(temp$Gene %in% proliGene, 'red','grey50')
temp$col<-ifelse(temp$Gene %in% Il1bGene, 'green',temp$col)

plot(temp$logCPM,temp$logFC,xlab="LogCPM",ylab="logFC",col=temp$col,main=compName)
abline(h=0,col="grey",lty=2)
targets_df <- temp[temp$Gene %in% target_genes, ]
if(nrow(targets_df) > 0){
  text(targets_df$logCPM, targets_df$logFC, 
       labels = targets_df$Gene, 
       pos = 3,      
       cex = 0.8,    
       col = "red") 
}

dev.off()


  
#############################################################################################################
############                                                                                       ##########
############                                                                                       ##########
#############################################################################################################
sigGene<-lapply(De.table,function(x){x[x$FDR<0.05 & abs(x$logFC)>1,]})[[1]];head(sigGene)
dim(sigGene)
nrow(sigGene[sigGene$logFC>0,])
nrow(sigGene[sigGene$logFC<0,])

sigGene$Gene<-sapply(strsplit(row.names(sigGene),split="\\_"),function(y) y[2])

stat_spl<-De.table[[1]]
stat_spl$Gene<-sapply(strsplit(row.names(stat_spl),split="\\_"),function(y) y[2])
stat_spl[stat_spl$Gene %in% c("Il1b","Nlrp3"),]

sum(pdata$ID==colnames(cpm_data))
cpm_2gene<-data.frame(pdata,t(cpm_data[row.names(stat_spl[stat_spl$Gene %in% c("Il1b","Nlrp3"),]),]));head(cpm_2gene)
colnames(cpm_2gene)[4:5]<-c("Il1b","Nlrp3")


sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/Il1bNlrp3.tiff",outputFolder), res=600, width=8,height=4,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/Il1bNlrp3.pdf",outputFolder),width=8,height=4)}
  par(mfrow=c(1,2))
  boxplot(Il1b~trt,cpm_2gene,main="Il1b")
  boxplot(Nlrp3~trt,cpm_2gene,main="Nlrp3")
  dev.off()
})



########################################################################################
#             heat map
########################################################################################
library(gplots)
head(cpm_data)
fdrCut=0.05
logFCcut=1
lapply(De.table,function(x) sum(x$FDR<fdrCut & abs(x$logFC)>logFCcut))
inputGeneList<-lapply(De.table,function(x) row.names(x[(x$FDR<fdrCut & abs(x$logFC)>logFCcut),]))
lapply(inputGeneList,function(x) length(x))

unique(pdata$trt)
sum(colnames(cpm_data)==pdata$ID)
sigPlt<-list(Sp.vs.BM=cpm_data[inputGeneList[[1]],pdata[,"ID"]])
lapply(sigPlt,function(x) dim(x))

###########heatmap####################
RNAseq_heatmap_update<-function(plotDat,outN=outputFolder){
  row.names(plotDat)<-sapply(strsplit(row.names(plotDat),split="\\_"),function(x) x[2])
  color.map<-rep(c("green","magenta"),each=5)
  mini<-min(plotDat);max<-max(plotDat)
  len_bk=100
  bk<-seq(mini,max,by=(max-mini)/len_bk)
  mycols<-c(colorRampPalette(colors = c("darkblue","blue","white"))(length(bk)),
            colorRampPalette(colors = c("white","red","darkred"))(length(bk)))
  dist.my <- function(x) as.dist(1-cor(t(x)))
  hclust.my <- function(x) hclust(x, method="complete")
  for(ip in 1:2){
    if(ip==1){tiff (sprintf("%s.tiff",outN), res=600, width=10,height=10,units="in",compression="lzw")}
    if(ip==2){pdf(sprintf("%s.pdf",outN),width=10,height=10)}
    sidebarcolors <- color.map
    heatmap.2(as.matrix(plotDat),Rowv=TRUE,Colv=TRUE,cexRow=0.75,cexCol=1.2,labRow=FALSE,
              distfun=dist.my,hclustfun = hclust.my,ColSideColors=sidebarcolors,
              dendrogram=c("both"),col=mycols,key=TRUE,keysize=0.6,symkey=TRUE,scale="row",
              density.info="none",trace="none",margins=c(8,6),main="")
    dev.off()
  }
}

RNAseq_heatmap_update(plotDat=sigPlt[[1]],outN=sprintf("%s/SpvsBM",outputFolder))

########################################################################################
#             GSEA analysis
########################################################################################
library(msigdbr)
library(fgsea)
library(data.table)

source("/home/yyw9094/code/sourceRcode/sourceFgsea.R")
##################################################################
unique(as.data.frame(msigdbr(species = "Homo sapiens", category ="C5"))$gs_subcat)
unique(as.data.frame(msigdbr(species = "Mus musculus", category ="C5"))$gs_subcat)
##################################################################
bp_list<-lapply(De.table,function(x){
  bp_gsea<-gsea_Ana_edgeR(species = "Mus musculus",
                          category = "C5",
                          subcategory = "GO:BP",
                          DEinput=x,
                          cut_padj=0.05)
  return(bp_gsea)
})
saveRDS(bp_list, file=sprintf("%s/BP_list.ALL_resting.RDS",SaveFolder))

lapply(bp_list,function(x) {sum(x$eaRes$padj<0.05)}) #


####################
bp_SpVsBm<-bp_list[[1]]$collDat;head(bp_SpVsBm)
plot(bp_SpVsBm$NES, -log10(bp_SpVsBm$padj))

bmPath<-c(bp_SpVsBm$pathway[grep('CYCLE', bp_SpVsBm$pathway)], 
  bp_SpVsBm$pathway[grep('CHROMOSOME', bp_SpVsBm$pathway)], 
  bp_SpVsBm$pathway[grep('REPLICATION', bp_SpVsBm$pathway)]
)

(bmPath<-bmPath[!bmPath %in% c("GOBP_NEGATIVE_REGULATION_OF_CHROMOSOME_ORGANIZATION",
                               "GOBP_MEIOTIC_CELL_CYCLE",
                              "GOBP_DNA_STRAND_ELONGATION_INVOLVED_IN_DNA_REPLICATION",
                              "GOBP_DNA_DEPENDENT_DNA_REPLICATION",
                              "GOBP_REGULATION_OF_CELL_CYCLE_CHECKPOINT")])


spPath<-c(bp_SpVsBm$pathway[grep('TOLL_LIKE', bp_SpVsBm$pathway)],
  bp_SpVsBm$pathway[grep('INTERLEUKIN', bp_SpVsBm$pathway)],
  bp_SpVsBm$pathway[grep('CYTOKINE', bp_SpVsBm$pathway)])

bp_SpVsBm[bp_SpVsBm$pathway %in% bmPath,]
bp_SpVsBm[bp_SpVsBm$pathway %in% spPath,]

(plt_BP<-bp_SpVsBm[bp_SpVsBm$pathway %in% c(bmPath,spPath),])


plt_BP$label <- ifelse(plt_BP$padj<0.001, 'FDR<0.001', paste("FDR=",round(plt_BP$padj,3)))

library(ggplot2)
ggplot(plt_BP, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = label, 
                hjust = ifelse(NES > 0, -0.1, 1.1)), 
            size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) + 
  labs(x = "Pathway", title = "GSEA Results")



bpPlot<-ggplot(plt_BP, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.3) +
  geom_hline(yintercept = 0, color = "black") +
  geom_text(aes(label = label, 
                vjust = ifelse(NES > 0, -0.5, 1.5)), 
            size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(colour = "black"), # Keep the vertical axis line
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) + 
  labs(x = "Pathway", title = "GSEA Results")


pdf(sprintf("%s/Fig1B_BP.Pathway.pdf",outputFolder),width=12,height=12)
print(bpPlot)
dev.off()


library(ggplot2)
lapply(names(bp_list),function(x){
  print(x)
  temp<-bp_list[[x]]
  fwrite(temp$eaRes, file=sprintf("%s/BP.%s.xls",outputFolder,x), sep="\t", sep2=c("", " ", ""))
  fwrite(temp$collDat, file=sprintf("%s/BP_collDat.%s.xls",outputFolder,x), sep="\t", sep2=c("", " ", ""))
})


#######
bp_list<-lapply(De.table,function(x){
  bp_gsea<-gsea_Ana_edgeR(species = "Mus musculus",
                          category = "C5",
                          subcategory = "GO:BP",
                          DEinput=x,
                          cut_padj=0.05)
  return(bp_gsea)
})

lapply(bp_list,function(x) {sum(x$eaRes$padj<0.05)})

library(ggplot2)
lapply(names(bp_list),function(x){
  print(x)
  temp<-bp_list[[x]]
  fwrite(temp$eaRes, file=sprintf("%s/BP_resting.%s.xls",outputFolder,x), sep="\t", sep2=c("", " ", ""))
  fwrite(temp$collDat, file=sprintf("%s/BP_collDat_resting.%s.xls",outputFolder,x), sep="\t", sep2=c("", " ", ""))
})


#######Reactome
Rec_list<-lapply(De.table,function(x){
  out_gsea<-gsea_Ana_edgeR(species = "Mus musculus",
                          category = "C2",
                          subcategory = "CP:REACTOME",
                          DEinput=x,
                          cut_padj=0.05)
  return(out_gsea)
})
lapply(Rec_list,function(x) {sum(x$eaRes$padj<0.05)})
library(ggplot2)
lapply(names(Rec_list),function(x){
  print(x)
  temp<-Rec_list[[x]]
  fwrite(temp$eaRes, file=sprintf("%s/SpvsBM_REACTOME.%s.xls",outputFolder,x), sep="\t", sep2=c("", " ", ""))
})


########################################################################################
#             gene plot
########################################################################################
geneToPlot<-"Tgf"
head(cpm_data)
tgfGene<-row.names(cpm_data)[grepl("Tgf",row.names(cpm_data))]
tgfData<-t(cpm_data[tgfGene[1:4],]);head(tgfData)
row.names(tgfData)
df_tgf<-data.frame(Organ=rep(c("BM","Sp"),each=5),tgfData);head(df_tgf)
sapply(1:length(tgfGene[1:4]),function(x){
  temp<-df_tgf[,c("Organ",tgfGene[x])]
  colnames(temp)[2]<-"Gene"
  
  pTemp<-ggplot(temp, aes(x = factor(Organ), y = Gene, fill = Organ)) +
    geom_boxplot(alpha = 0.80) +
    geom_point(aes(fill = Organ), size = 5, shape = 21, position = position_jitterdodge()) +
    ggtitle(tgfGene[x])+
    theme(text = element_text(size = 18),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  pdf(sprintf("%s/%s.pdf",outputFolder,tgfGene[x]), width=6,height=4)
  print(pTemp)
  dev.off()
})




targetGene<-"Csf1r"
geneToPlot<-"Csf1r"
head(cpm_data)
tgfGene<-row.names(cpm_data)[grepl("Csf1r",row.names(cpm_data))]
tgfData<-t(cpm_data[tgfGene[1],,drop=F]);head(tgfData)
row.names(tgfData)
df_tgf<-data.frame(Organ=rep(c("BM","Sp"),each=5),tgfData);head(df_tgf)

  temp<-df_tgf[,c("Organ",tgfGene[1])]
  colnames(temp)[2]<-"Gene"
  
  pTemp<-ggplot(temp, aes(x = factor(Organ), y = Gene, fill = Organ)) +
    geom_boxplot(alpha = 0.80) +
    geom_point(aes(fill = Organ), size = 5, shape = 21, position = position_jitterdodge()) +
    ggtitle(tgfGene[1])+
    theme(text = element_text(size = 18),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  pdf(sprintf("%s/%s.pdf",outputFolder,tgfGene[1]), width=6,height=4)
  print(pTemp)
  dev.off()








