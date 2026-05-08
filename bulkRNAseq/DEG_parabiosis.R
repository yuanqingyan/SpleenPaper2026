
inputFolder<-"./March2025/mouse/RawCount"
outputFolder<-"./DEG"
SaveFolder<-"./March2025/mouse/RDS"
if (!dir.exists(SaveFolder)) {dir.create(SaveFolder, recursive = TRUE)}
source("~/code/sourceRcode/bulkRNAseq_sourceCode.R")


#############################################################################################################
############                                                                                       ##########
#############################################################################################################
list.files(inputFolder)
dat<-read.delim(sprintf("%s/AllGeneExpCount.txt",inputFolder),header=T,sep="\t",row.names=1,stringsAsFactors = F)
colnames(dat)<-sapply(strsplit(colnames(dat),split="\\DEout."),function(x) x[2])
colnames(dat)<-sapply(strsplit(colnames(dat),split="\\_S"),function(x) x[1])
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
anoFile<-read.delim("~/code/sourceRcode/sourceData/mart_export_grcm38.p6.txt",header=T,sep="\t",stringsAsFactors=FALSE);head(anoFile)
anoFile$ID<-paste(anoFile[,1],anoFile[,3],sep="_")
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]#;

dat_coding<-dat[sapply(strsplit(row.names(dat),split="\\_"),function(x) x[1]) %in% codingGene$Gene.stable.ID,];dim(dat_coding)
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
pdata<-data.frame(ID=colnames(dat_coding),
                  trt=sapply(strsplit(colnames(dat_coding),split="\\_"),function(x) x[1]));pdata
unique(pdata$trt)
pdata$col<-c("green","purple","blue",'red')[factor(pdata$trt)]
pdata$Grp<-c("Grp4","Grp2","Grp3",'Grp1')[factor(pdata$trt)]

#############################################################################################################
############                                                                                       ##########
#############################################################################################################
saveRDS(pdata, file=sprintf("%s/pdataUse.RDS",inputFolder))
pdata$Grp <- factor(pdata$Grp, levels=c("Grp1","Grp2",'Grp3','Grp4'))

dat_input<-dat_coding[,paste(pdata$ID)];dim(dat_input) #
dat_input<-dat_input[rowSums(dat_input)>0,];dim(dat_input)#
sum(colnames(dat_input)==pdata$ID)
boxplot(log1p(dat_input),las=2)
##############
library(edgeR)
edgeInput1<-DGEList(dat_input,group=as.factor(pdata$Grp))
edgeInput2 <- edgeInput1[rowSums(1e+06*edgeInput1$counts/expandAsMatrix(edgeInput1$samples$lib.size,dim(edgeInput1))>=1)>=ncol(edgeInput1)*0.1, ]
edgeInput2 <- calcNormFactors(edgeInput2)
edgeInput2$samples


pdf(sprintf("%s/MDS_byTreatment.pdf",outputFolder),  width=6,height=6)
plotMDS(cpm(edgeInput2,prior.count=1,log=TRUE), col=pdata$col,pch=19,cex=2,xlab="Dim1",ylab="Dim2")
dev.off()


testCPM<-cpm(edgeInput2,prior.count=1,log=TRUE);head(testCPM)

plt1<-pdata[pdata$Grp %in% c("Grp4",'Grp3'),];plt1
testCPM1<-testCPM[,plt1$ID];head(testCPM1)
plotMDS(testCPM1, col=plt1$col,pch=19,cex=2,xlab="Dim1",ylab="Dim2")

pdf(sprintf("%s/MDS.Parabionet_CD45.1InHost_Vs_CD45.2InHost.pdf",outputFolder),  width=6,height=6)
plotMDS(testCPM1, col=plt1$col,pch=19,cex=2,xlab="Dim1",ylab="Dim2")
dev.off()


#############################################################################################################
############                                                                                       ##########
#############################################################################################################
pTrt<-pdata
design.mat <- model.matrix(~pTrt$Grp)
library(stringr)
(colnames(design.mat) <-str_extract(colnames(design.mat), "Grp\\d+"))
colnames(design.mat)[c(1)]<- c("(Intercept)")
rownames(design.mat) <- colnames(edgeInput2)
cpm_data<-cpm(edgeInput2,prior.count=1,log=TRUE)
write.csv(cpm_data,file=sprintf("%s/cpm.csv",outputFolder))
head(cpm_data)

edgeInput3 <- estimateDisp(edgeInput2, design.mat, robust=TRUE) #
fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE)

Grp2_vs_Grp1<-glmQLFTest(fit, coef=2)
Grp3_vs_Grp1<-glmQLFTest(fit, coef=3)
Grp4_vs_Grp1<-glmQLFTest(fit, coef=4)
Grp3_vs_Grp2 <- glmQLFTest(fit, contrast=c(0,-1,1,0))
Grp4_vs_Grp2<- glmQLFTest(fit, contrast=c(0,-1,0,1))
Grp4_vs_Grp3 <- glmQLFTest(fit, contrast=c(0,0,-1,1))



DE.ALL<-list(Grp2_vs_Grp1=Grp2_vs_Grp1,
             Grp3_vs_Grp1=Grp3_vs_Grp1,
             Grp4_vs_Grp1=Grp4_vs_Grp1,
             Grp3_vs_Grp2=Grp3_vs_Grp2,
             Grp4_vs_Grp2=Grp4_vs_Grp2,
             Grp4_vs_Grp3=Grp4_vs_Grp3)
  
De.table<-lapply(DE.ALL,function(xy) {
    temp<-xy$table
    temp$FDR<-p.adjust(temp$PValue,method="BH")
    temp$stat_gsea<-zscoreT(sign(temp$logFC)*sqrt(temp$F),df=xy$df.total);return(temp)})
lapply(De.table,function(xy) head(xy))
saveRDS(De.table,file=sprintf("%s/TLF_DeTable.rds",outputFolder))
  
FDRCut=0.05;logCut=log2(1.5)
singCount<-do.call("rbind",lapply(De.table,function(xy) {
    tSig<-sum(xy$FDR<FDRCut & abs(xy$logFC)>logCut)
    uSig<-sum(xy$FDR<FDRCut & xy$logFC>logCut)
    dSig<-sum(xy$FDR<FDRCut & xy$logFC<(-logCut))
    out<-c(tSig,uSig,dSig)}))
colnames(singCount)<-c("TotalSig","UpSig","DnSig");singCount
write.table(data.frame(Comparison=row.names(singCount),singCount),file=sprintf("%s/TLR_singCount.xls",outputFolder),col.names=T,row.names=F,sep="\t")
  

###write stat out
lapply(1:length(De.table),function(xy){
  compName<-names(De.table)[xy]
  temp<-De.table[[xy]]
  write.table(data.frame(Gene=row.names(temp),temp),file=sprintf("%s/tlr_StatOut.%s.xls",outputFolder,compName),row.names=F,col.names=T,sep="\t")
})

#############################################################################################################
length(De.table)

sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/TLR_Volcano.tiff",outputFolder), res=600, width=10,height=10,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/TLR_Volcano.pdf",outputFolder),width=10,height=10)}
  par(mfrow=c(2,2))
  lapply(1:length(De.table),function(xy){
    compName<-names(De.table)[xy]
    temp<-De.table[[xy]]
    plot(temp$logFC,-log10(temp$FDR),xlab="Log2 FC",ylab="-log10 FDR",main=compName)
    abline(h=-log10(FDRCut),col="red",lty=2)})
  dev.off()
})

#############################################################################################################
sigGene<-lapply(De.table,function(xy){xy[xy$FDR< FDRCut & abs(xy$logFC)>logCut,]});head(sigGene[[1]])
sigGene<-lapply(sigGene,function(xy) {
  xy$gene<-sapply(strsplit(row.names(xy),split="\\_"),function(y) y[2])
  return(xy)})




De.table<-readRDS(file=sprintf("%s/TLF_DeTable.rds",outputFolder))
deg_g4vsg3<-De.table$Grp4_vs_Grp3;head(deg_g4vsg3)
deg_g4vsg3<-deg_g4vsg3[order(deg_g4vsg3$PValue),]
deg_g4vsg3$gene<-sapply(strsplit(row.names(deg_g4vsg3),split="\\_"),function(x) x[2])
deg_g4vsg3_clean<-deg_g4vsg3[!duplicated(deg_g4vsg3$gene),]
deg_g4vsg3_clean[deg_g4vsg3_clean$gene %in% c("Il1b",'Mki67','Tlr7'),]
sum(deg_g4vsg3_clean$FDR<0.05)

cpm_data<-cpm_data[row.names(cpm_data) %in% row.names(deg_g4vsg3_clean),]
row.names(cpm_data)<-sapply(strsplit(row.names(cpm_data),split="\\_"), function(x) x[2]); head(cpm_data)


sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/Parabionet_MA_45.1CMinHost_vs_45.2TLRCMinHost.tiff",outputFolder), res=600, width=10,height=10,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/Parabionet_MA_45.1CMinHost_vs_45.1TLRCMinHost.pdf",outputFolder),width=10,height=10)}
  xy<-6
  compName<-names(De.table)[xy]
  temp<-De.table[[xy]]
  temp$col<-ifelse(temp$FDR<0.05, 'red','black')
  plot(temp$logCPM, temp$logFC,ylab="Log2 FC",xlab="logCPM", col=temp$col,main=compName)
  dev.off()
  })
  

############
sum(pTrt$ID==colnames(cpm_data));dim(pTrt)
df_p<-data.frame(pTrt,t(cpm_data[c("Il1b",'Mki67','Tlr7'),]));head(df_p)
df_p2<-df_p[df_p$Grp %in% c('Grp3','Grp4','Tlr7'),];df_p2

library(ggplot2)
bp_Plt<-function(GeneID="Il1b",colo=c("blue",'green')){
  df.plot<-data.frame(exp=as.numeric(df_p2[,GeneID]),group=df_p2$trt)
  bp_col<-colo
  
  set.seed(123)
  p <- ggplot(df.plot, aes(x=group, y=exp,col=group)) + 
    geom_boxplot()+
    geom_jitter(shape=1, position=position_jitter(0.3),size=3)+ 
    scale_color_manual(values=bp_col)+
    labs(y = "CPM")+
    labs(x = "")+
    ggtitle(sprintf("%s",GeneID))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.5,size=12,colour = "black"),
          axis.title=element_text(size=12,colour = "black"),
          axis.text.y = element_text(size=12,colour = "black"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position = "none") 
  return(p)
}

print(bp_Plt(GeneID="Il1b",colo=c('green',"blue")))

width=3;height=3
pltG<-c('Il1b','Mki67','Tlr7')
lapply(pltG,function(x){
  pdf(sprintf("%s/Parabio_BoxPlt_%s.pdf", outputFolder,x), width=width,height=height)
  p<-bp_Plt(GeneID=x,colo=c('green',"blue"))
  print(p)
  dev.off()
})
write.csv(cpm_data, file=sprintf("%s/parabio_cpm.csv", outputFolder))

########################################################################################
#             GSEA analysis
########################################################################################
library(msigdbr)
library(fgsea)
library(data.table)

gsea_Ana_edgeR<-function(species = "Homo sapiens",
                         category = "C2",
                         subcategory = "CP:REACTOME",
                         DEinput_df=DEinput_df, ##
                         cut_padj=0.05){
  temp<-DEinput_df
  temp<-temp[order(temp$PValue),]
  temp_df<-data.frame(temp,gene=sapply(strsplit(row.names(temp),split="\\_"),function(xm) xm[2]))
  temp_df<-temp_df[!duplicated(temp_df$gene),]
  row.names(temp_df)<-temp_df$gene
  temp_df<-temp_df[order(-temp_df$stat_gsea),]
  
  if(is.null(subcategory)){
    geneSets <- msigdbr(species = species, category = category)
  }else{
    geneSets <- msigdbr(species = species, category = category, subcategory = subcategory)
  }

  ### filter background to only include genes that we assessed.
  geneSets <- geneSets[geneSets$gene_symbol %in% row.names(temp_df),]
  m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name);length(m_list)

  stats<-temp_df$stat_gsea
  names(stats) <- rownames(temp_df)

  set.seed(1234)
  eaRes <- fgsea(pathways = m_list, stats = stats, nperm = 5e4, minSize = 10)
  eaRes <-eaRes[!is.na(eaRes$padj),]
  eaRes <-eaRes[order(eaRes$padj),]
  sum(eaRes$padj<0.05)

  collapsedPathways <- collapsePathways(eaRes[order(pval)][padj < cut_padj],m_list, stats)
  collDat <- eaRes[pathway %in% collapsedPathways$mainPathways][order(padj),]
  out<-list(eaRes=eaRes,m_list=m_list,stats=stats,collDat=collDat)
  return(out)
}
##################################################################
unique(as.data.frame(msigdbr(species = "Homo sapiens", category ="C5"))$gs_subcat)
unique(as.data.frame(msigdbr(species = "Mus musculus", category ="C5"))$gs_subcat)
##################################################################

hall_list<-lapply(De.table,function(xbp){
  H_gsea<-gsea_Ana_edgeR(species = "Mus musculus",
                          category = "H",
                         subcategory = NULL,
                          DEinput_df=xbp,
                          cut_padj=0.05)
  return(H_gsea)
});names(hall_list)<-names(De.table)

H_grp4vsgrp3<-hall_list$Grp4_vs_Grp3
H_grp4vsgrp3$eaRes
H_grp4vsgrp3$eaRes[H_grp4vsgrp3$eaRes$pathway %in% c("HALLMARK_G2M_CHECKPOINT","HALLMARK_INFLAMMATORY_RESPONSE"),]

fwrite(H_grp4vsgrp3$eaRes, file=sprintf("%s/Hall_Parabio_Cd45.1InHost_vs_Tlr7CD45.2InHost.xls",outputFolder), sep="\t", sep2=c("", " ", ""))

library(ggplot2)
plt<-plotEnrichment(H_grp4vsgrp3$m_list[["HALLMARK_G2M_CHECKPOINT"]], H_grp4vsgrp3$stats) +
  labs(title='G2M_CHECKPOINT')

pdf(sprintf("%s/parabio_GSEA_grp4Vsgrp3.pdf",outputFolder),width=7.87,height=7.87, family="ArialMT")
print(plt)
dev.off()

pdf(sprintf("%s/parabio_GSEA_Cd45.1InHost_Vs_Tlr7CD45.2InHost.pdf",outputFolder),width=7.87,height=7.87, family="ArialMT")
print(plt)
dev.off()

