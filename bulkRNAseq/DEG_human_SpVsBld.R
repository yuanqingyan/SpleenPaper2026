inputFolder<-"./Mar2024/human/RawCount"
outputFolder<-"./DEG"
SaveFolder<-"./Mar2024/human/RDS"
source("~/code/sourceRcode/bulkRNAseq_sourceCode.R")
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
list.files(inputFolder)
dat<-read.delim(sprintf("%s/AllGeneExpCount.txt",inputFolder),header=T,sep="\t",row.names=1,stringsAsFactors = F)
colnames(dat)<-sapply(strsplit(colnames(dat),split="\\."),function(x) sprintf("%s_%s",x[8],x[9]))
#############################################################################################################
############                           protein coding only                                         ##########
#############################################################################################################
anoFile<-read.delim("~/reference/RNAseqRef/RefDownload/mart_export_ensemble110_GRCh38.p14.txt",header=T,sep=",",stringsAsFactors=FALSE);head(anoFile)
anoFile$ID<-paste(anoFile[,1],anoFile[,3],sep="_")
anoFile$Gene.name<-gsub("-",".",anoFile$Gene.name)
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]

dat_coding<-dat[sapply(strsplit(row.names(dat),split="\\_"),function(x) x[1]) %in% codingGene$Gene.stable.ID,];dim(dat_coding)#19797
#############################################################################################################
#############################################################################################################
pdata<-data.frame(ID=colnames(dat_coding),
                  grp=sapply(strsplit(colnames(dat_coding),split="\\_"),function(x) substr(x[1],1,2)));pdata
pdata$trt<-dplyr::case_when(
  pdata$grp==7 ~ 'Sp',
  pdata$grp==8 ~ 'Bl',
  .default="unknown"
);table(pdata$trt)
pdata$trt <- factor(pdata$trt, levels=c("Bl","Sp"))
pdata$col<-c("green","purple")[factor(pdata$trt)]
pdata$patient<-factor(rep(1:6,2));head(pdata)
#############################################################################################################
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
############.                          design matrix                                               ##########
#############################################################################################################
design.mat <- model.matrix(~pdata$trt+pdata$patient)
colnames(design.mat) <- sapply(strsplit(colnames(design.mat),split="\\$"),function(x) x[2])
colnames(design.mat)[c(1,2)]<- c("(Intercept)","Sp")
rownames(design.mat) <- colnames(edgeInput2)
design.mat

#############################################################################################################
############.                      check batch effect                                              ##########
#############################################################################################################
cpm_data<-cpm(edgeInput2,prior.count=1,log=TRUE)
listCPM<-list(cpm=cpm_data,edgeInput2=edgeInput2)
boxplot(cpm_data,col=pdata$col)

tiff(sprintf("%s/MDS_noLabel_human_SpVsBld.tiff",outputFolder), res=600, width=6,height=6,units="in",compression="lzw")
plotMDS(cpm_data, col=c('purple','green4')[factor(pdata$trt)],pch=19,cex=2,xlab="Dim1",ylab="Dim2")
dev.off()

pdf(sprintf("%s/MDS_noLabel_human_SpVsBld.pdf",outputFolder), width=6,height=6)
plotMDS(cpm_data, col=c('purple','green4')[factor(pdata$trt)],pch=19,cex=2,xlab="Dim1",ylab="Dim2")
dev.off()

tiff(sprintf("%s/MDS_Label_human_SpVsBld.tiff",outputFolder), res=600, width=6,height=6,units="in",compression="lzw")
plotMDS(cpm_data,pch=19,cex=1,xlab="Dim1",ylab="Dim2",labels=pdata$ID)
dev.off()

#############################################################################################################
#############################################################################################################
edgeInput3 <- estimateDisp(edgeInput2, design.mat, robust=TRUE) ##
edgeInput3$common.dispersion
plotBCV(edgeInput3)
fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE)
plotQLDisp(fit)

#############################################################################################################
############.                                                                                      ##########
#############################################################################################################
SpvsBl<-glmQLFTest(fit, coef=2)

DE.ALL<-list(SpvsBl=SpvsBl)
names(DE.ALL);lapply(DE.ALL,function(x) summary(decideTests(x)))

De.table<-lapply(DE.ALL,function(x) {
  temp<-x$table
  temp$FDR<-p.adjust(temp$PValue,method="BH")
  temp$stat_gsea<-zscoreT(sign(temp$logFC)*sqrt(temp$F),df=x$df.total);return(temp)})
lapply(De.table,function(x) head(x))
saveRDS(De.table,file=sprintf("%s/SpVsBl_DeTable.rds",SaveFolder))

FDRCut=0.05;logCut=log2(1.5)
singCount<-do.call("rbind",lapply(De.table,function(x) {
  tSig<-sum(x$FDR<FDRCut & abs(x$logFC)>logCut)
  uSig<-sum(x$FDR<FDRCut & x$logFC>logCut)
  dSig<-sum(x$FDR<FDRCut & x$logFC<(-logCut))
  out<-c(tSig,uSig,dSig)}))
colnames(singCount)<-c("TotalSig","UpSig","DnSig");singCount
write.table(data.frame(Comparison=row.names(singCount),singCount),file=sprintf("%s/singCount.xls",outputFolder),col.names=T,row.names=F,sep="\t")

lapply(1:length(De.table),function(x){
  compName<-names(De.table)[x]
  temp<-De.table[[x]]
  write.table(data.frame(Gene=row.names(temp),temp),file=sprintf("%s/StatOut.%s.xls",outputFolder,compName),row.names=F,col.names=T,sep="\t")
})

#############################################################################################################
#############################################################################################################
length(De.table)

sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/MA_human_SpvsBld.tiff",outputFolder), res=600, width=8,height=8,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/MA_human_SpvsBld.pdf",outputFolder),width=8,height=8)}
  lapply(1:length(De.table),function(xy){
    compName<-names(De.table)[xy]
    temp<-De.table[[xy]]
    temp$col<-ifelse(temp$FDR<0.05, 'red','black')
    plot(temp$logCPM, temp$logFC,ylab="Log2 FC",xlab="logCPM", col=temp$col)#,main=compName)
  })
  dev.off()
})

#############################################################################################################
#############################################################################################################
sigGene<-lapply(De.table,function(x){x[x$FDR<0.05 & abs(x$logFC)>log2(1.5),]})[[1]];head(sigGene)
sigGene$Gene<-sapply(strsplit(row.names(sigGene),split="\\_"),function(y) y[2])

stat_spl<-De.table[[1]]
stat_spl$Gene<-sapply(strsplit(row.names(stat_spl),split="\\_"),function(y) y[2])

sum(pdata$ID==colnames(cpm_data))


########################################################################################
#             heat map
########################################################################################
library(gplots)
head(cpm_data)
fdrCut=0.05
logFCcut=log2(1.5)
lapply(De.table,function(x) sum(x$FDR<fdrCut & abs(x$logFC)>logFCcut))
inputGeneList<-lapply(De.table,function(x) row.names(x[(x$FDR<fdrCut & abs(x$logFC)>logFCcut),]))
lapply(inputGeneList,function(x) length(x))

unique(pdata$trt)
sum(colnames(cpm_data)==pdata$ID)
sigPlt<-list(Sp.vs.Bl=cpm_data[inputGeneList[[1]],pdata[,"ID"]])
lapply(sigPlt,function(x) dim(x))


RNAseq_heatmap_update2<-function(plotDat=sigPlt[[1]],outN=outputFolder){
  row.names(plotDat)<-sapply(strsplit(row.names(plotDat),split="\\_"),function(x) x[2])
  color.map<-rep(c("green","magenta"),each=6)
  mini<-min(plotDat);max<-max(plotDat)
  len_bk=100
  bk<-seq(mini,max,by=(max-mini)/len_bk)
  mycols<-c(colorRampPalette(colors = c("darkblue","blue","white"))(length(bk)),
            colorRampPalette(colors = c("white","red","darkred"))(length(bk)))
  dist.my <- function(x) as.dist(1-cor(t(x)))
  hclust.my <- function(x) hclust(x, method="complete")
  for(ip in 1:2){
    if(ip==1){tiff (sprintf("%s.tiff",outN), res=600, width=12,height=12,units="in",compression="lzw")}
    if(ip==2){pdf(sprintf("%s.pdf",outN),width=12,height=12)}
    sidebarcolors <- color.map
    heatmap.2(as.matrix(plotDat),Rowv=TRUE,Colv=TRUE,cexRow=0.2,cexCol=1.2,labRow=row.names(plotDat),
              distfun=dist.my,hclustfun = hclust.my,ColSideColors=sidebarcolors,
              dendrogram=c("both"),col=mycols,key=TRUE,keysize=0.6,symkey=TRUE,scale="row",
              density.info="none",trace="none",margins=c(8,6),main="")
    dev.off()
  }
}

RNAseq_heatmap_update2(plotDat=sigPlt[[1]],outN=sprintf("%s/SpvsBM_v2",outputFolder))

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
                         cut_padj=0.05,
                         iseed=1234){
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
  
  set.seed(iseed)
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

H_list<-lapply(De.table,function(x){
 H_gsea<-gsea_Ana_edgeR(species = "Homo sapiens",
                          category = "H",
                          subcategory = NULL,
                          DEinput_df=x,
                          cut_padj=1,
                        iseed=23)
  return(H_gsea)
});names(H_list)<-paste('HallMark',names(De.table),sep="_")
lapply(H_list,function(x) {sum(x$eaRes$padj<0.05)})
H_SpvsBl<-H_list$HallMark_SpvsBl$eaRes
H_SpvsBl[H_SpvsBl$pathway %in% c("HALLMARK_INFLAMMATORY_RESPONSE"),]


humanList_hall<-list(cpm_data=cpm_data,
                pdata=pdata,
                DE.ALL=DE.ALL,
                De.table=De.table,
                sigGene=sigGene,
                H_list=H_list)
saveRDS(humanList_hall,file=sprintf("%s/humanSPDataList_HallMark.RDS",SaveFolder))



#########
humanList_hall<-readRDS(file=sprintf("%s/humanSPDataList_HallMark.RDS",SaveFolder))
DeTab<-humanList_hall$De.table$SpvsBl
H_list<-humanList_hall$H_list


library(ggplot2)
(plt<-plotEnrichment(H_list$HallMark_SpvsBl$m_list[["HALLMARK_INFLAMMATORY_RESPONSE"]], H_list$HallMark_SpvsBl$stats) +
    labs(title='Inflammatory response'))

pdf(sprintf("%s/Human_SpvBld__Inflam.pdf",outputFolder),width=7.87,height=7.87/1.4, family="ArialMT")
print(plt)
dev.off()


########################################################################################
#             gene plot
########################################################################################
humanList<-readRDS(file=sprintf("%s/humanSPDataList_HallMark.RDS",SaveFolder))


geneToPlot<-c("IL1B","NLRP3")
cpm_data<-humanList$cpm_data;head(cpm_data)
pdata<-humanList$pdata
row.names(cpm_data)[grep(geneToPlot[1],row.names(cpm_data))]
row.names(cpm_data)[grep(geneToPlot[2],row.names(cpm_data))]


PltGeneIndex<-row.names(cpm_data)[sapply(geneToPlot,function(x) {grep(x,row.names(cpm_data))})]
PltGeneData<-t(cpm_data[PltGeneIndex,]);head(PltGeneData)
row.names(PltGeneData)==pdata$ID
df_plt<-data.frame(Organ=pdata$trt,PltGeneData);head(df_plt)

library(ggplot2)
sapply(2:ncol(df_plt),function(x){
  temp<-df_plt[,c("Organ",colnames(df_plt)[x])]
  colnames(temp)[2]<-"Gene"

  pTemp<-ggplot(temp, aes(x = factor(Organ), y = Gene, fill = Organ)) +
    geom_boxplot(alpha = 0.80) +
    geom_point(aes(fill = Organ), size = 5, shape = 21, position = position_jitterdodge()) +
    ggtitle(colnames(df_plt)[x])+
    theme(text = element_text(size = 18),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  pdf(sprintf("%s/Human_%s.pdf",outputFolder,colnames(df_plt)[x]), width=6,height=4)
  print(pTemp)
  dev.off()
})




