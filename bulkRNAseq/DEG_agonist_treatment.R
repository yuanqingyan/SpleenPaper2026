
# Pam3CSK4 – TLR1/2
# HKLM – TLR2
# Poly(I:C) HMW – TLR3
# LPS-EK – TLR4
# FLA-ST – TLR5
# FSL-1 – TLR2/6
# R837, R848, ssRNA40 – TLR7 agonists
# ODN 1826 – TLR9
# 
# For TLR7, we used three different agonists (R837, R848, and ssRNA40). 


inputFolder<-"./Nov2025/mouse/RawCount"
outputFolder<-"./DEG"
SaveFolder<-"./Nov2025/mouse/RDS"
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
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]

dat_coding<-dat[sapply(strsplit(row.names(dat),split="\\_"),function(x) x[1]) %in% codingGene$Gene.stable.ID,];dim(dat_coding)#
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
pdata<-data.frame(ID=colnames(dat_coding),
                  organ=sapply(strsplit(colnames(dat_coding),split="\\_"),function(x) x[1]),
                  trt=sapply(strsplit(colnames(dat_coding),split="\\_"),function(x) x[2]));pdata
unique(pdata$trt)
unique(pdata$organ)
table(pdata$trt,pdata$organ)
pdata$trt<-factor(pdata$trt,c("NT",unique(pdata$trt[!pdata$trt %in% 'NT'])))#
pdata$Grp<-pdata$trt

pdata$col1<-c("green","purple")[factor(pdata$organ)]

library(RColorBrewer)
pdata$col2<-brewer.pal(n=length(unique(as.character(pdata$trt))),name="Paired")[pdata$trt]

#############################################################################################################
############                                                                                       ##########
#############################################################################################################
organ<-unique(pdata$organ)
trt<-unique(pdata$trt)

library(edgeR)
library(stringr)
library(msigdbr)
library(fgsea)
library(data.table)
library(ggplot2)
########################################
degFun<-function(dat_coding=dat_coding, pIn=pIn, name='BM', outputFolder=outputFolder){
  dat_input<-dat_coding[,paste(pIn$ID)] 
  dat_input<-dat_input[rowSums(dat_input)>0,]#
  ##############
  
  edgeInput1<-DGEList(dat_input,group=as.factor(pIn$Grp))
  edgeInput2 <- edgeInput1[rowSums(1e+06*edgeInput1$counts/expandAsMatrix(edgeInput1$samples$lib.size,dim(edgeInput1))>=1)>=ncol(edgeInput1)*0.1, ]
  edgeInput2 <- calcNormFactors(edgeInput2)

  pdf(sprintf("%s/MDS_byTreatment_%s.pdf",outputFolder,name),  width=6,height=6)
  plotMDS(cpm(edgeInput2,prior.count=1,log=TRUE), col=pIn$col2,pch=19,cex=2,xlab="Dim1",ylab="Dim2")
  dev.off()
  
  pTrt<-pIn
  design.mat <- model.matrix(~pTrt$Grp)
  (colnames(design.mat) <-sapply(strsplit(colnames(design.mat),split='Grp'),function(xxx) xxx[2]))
  colnames(design.mat)[c(1)]<- c("(Intercept)")
  rownames(design.mat) <- colnames(edgeInput2)
  cpm_data<-cpm(edgeInput2,prior.count=1,log=TRUE)
  write.csv(cpm_data,file=sprintf("%s/cpm_%s.csv",outputFolder, name))
  
  
  edgeInput3 <- estimateDisp(edgeInput2, design.mat, robust=TRUE) ##
  fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE)
  
  outList<-list(edgeInput3=edgeInput3,
                fit=fit,
                cpm_data=cpm_data,
                design.mat=design.mat)
  return(outList)
}


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
  
  geneSets <- msigdbr(species = species, category = category, subcategory = subcategory)
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
############

lapply(1:length(organ),function(iOrg){
  name<-organ[iOrg];print(name)
  pIn<-pdata[pdata$organ %in% name,]
  
  out1<-degFun(dat_coding=dat_coding, 
               pIn=pIn, 
               name=name, 
               outputFolder=outputFolder)
  
  fit=out1$fit
  out1$design.mat
  cpm_data=out1$cpm_data
  
  FLA.ST_vs_NT<-glmQLFTest(fit, coef=2)
  FSL_vs_NT<-glmQLFTest(fit, coef=3)
  HKML_vs_NT<-glmQLFTest(fit, coef=4)
  LPS_vs_NT<-glmQLFTest(fit, coef=5)
  ODN_vs_NT<-glmQLFTest(fit, coef=6) ###
  Pam_vs_NT<-glmQLFTest(fit, coef=7) 
  P_vs_NT<-glmQLFTest(fit, coef=8)
  R837_vs_NT<-glmQLFTest(fit, coef=9)
  R848_vs_NT<-glmQLFTest(fit, coef=10)
  ssRNA_vs_NT<-glmQLFTest(fit, coef=11)
  
  
  DE.ALL<-list(FLA.ST_vs_NT=FLA.ST_vs_NT,
               FSL_vs_NT=FSL_vs_NT,
               HKML_vs_NT=HKML_vs_NT,
               LPS_vs_NT=LPS_vs_NT,
               ODN_vs_NT=ODN_vs_NT,
               Pam_vs_NT=Pam_vs_NT,
               P_vs_NT=P_vs_NT,
               R837_vs_NT=R837_vs_NT,
               R848_vs_NT=R848_vs_NT,
               ssRNA_vs_NT=ssRNA_vs_NT)
  
  De.table<-lapply(DE.ALL,function(xy) {
    temp<-xy$table
    temp$FDR<-p.adjust(temp$PValue,method="BH")
    temp$stat_gsea<-zscoreT(sign(temp$logFC)*sqrt(temp$F),df=xy$df.total);return(temp)})
  lapply(De.table,function(xy) head(xy))
  saveRDS(De.table,file=sprintf("%s/TLF_DeTable_%s.rds",outputFolder,name))
  
  ###write stat out
  lapply(1:length(De.table),function(xy){
    compName<-names(De.table)[xy]
    temp<-De.table[[xy]]
    write.table(data.frame(Gene=row.names(temp),temp),file=sprintf("%s/%s_StatOut.%s.xls",outputFolder,name,compName),row.names=F,col.names=T,sep="\t")
  })
  
  
  ##################################################################
  unique(as.data.frame(msigdbr(species = "Homo sapiens", category ="C5"))$gs_subcat)
  unique(as.data.frame(msigdbr(species = "Mus musculus", category ="C5"))$gs_subcat)
  ##################################################################
  bp_list<-lapply(De.table,function(xbp){
    bp_gsea<-gsea_Ana_edgeR(species = "Mus musculus",
                            category = "C5",
                            subcategory = "GO:BP",
                            DEinput_df=xbp,
                            cut_padj=0.05)
    return(bp_gsea)
  });names(bp_list)<-names(De.table)
  
  lapply(bp_list,function(xt) {sum(xt$eaRes$padj<0.05)})
  
 
  lapply(names(bp_list),function(xy){
    print(xy)
    temp<-bp_list[[xy]]
    fwrite(temp$eaRes, file=sprintf("%s/%s_BP.%s.xls",outputFolder,name,xy), sep="\t", sep2=c("", " ", ""))
  })
  
  out<-list(bp_list=bp_list,
            De.table=De.table,
            cpm_data=cpm_data,
            pdata=pIn)
  saveRDS(out,file=sprintf("%s/%s_Mouse_Nov2025.RDS",SaveFolder,name))
})


########################################################
list.files(SaveFolder)

BM<-readRDS(sprintf("%s/BM_Mouse_Nov2025.RDS", SaveFolder))$De.table
SP<-readRDS(sprintf("%s/SP_Mouse_Nov2025.RDS", SaveFolder))$De.table

trtName<-c("R837",'R848','ssRNA')
BM_sel<-lapply(BM[c("R837_vs_NT", "R848_vs_NT", "ssRNA_vs_NT")], function(x) {x$gene<-row.names(x); x<-x[,c('gene','logFC','FDR')]; return(x)})
BM_sel<-lapply(1:length(BM_sel), function(x) {colnames(BM_sel[[x]])[2:3]<-sprintf("%s.%s", trtName[x], colnames(BM_sel[[x]])[2:3]); return(BM_sel[[x]])})

SP_sel<-lapply(SP[c("R837_vs_NT", "R848_vs_NT", "ssRNA_vs_NT")], function(x) {x$gene<-row.names(x); x<-x[,c('gene','logFC','FDR')]; return(x)})
SP_sel<-lapply(1:length(SP_sel), function(x) {colnames(SP_sel[[x]])[2:3]<-sprintf("%s.%s", trtName[x], colnames(SP_sel[[x]])[2:3]); return(SP_sel[[x]])})


# Merge the list by row names using Reduce and merge
merged_BM <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), BM_sel);head(merged_BM)
merged_SP <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), SP_sel);head(merged_SP)

rownames(merged_BM) <- merged_BM$gene; merged_BM$gene <- NULL;head(merged_BM)
rownames(merged_SP) <- merged_SP$gene; merged_SP$gene <- NULL;head(merged_SP)

selGene<-c("Cxcl10",'Il12a','Il6','Il1b','Csf3','Il1a','Tnf','Il15','Stat1','Nfkbia','Il12b','Ccl2','Irf1','Il18')
lapply(selGene, function(x) rownames(merged_BM)[grepl(x, rownames(merged_BM))])

pltGene<-c("ENSMUSG00000034855_Cxcl10",
           "ENSMUSG00000027776_Il12a",
           "ENSMUSG00000025746_Il6",
           "ENSMUSG00000027398_Il1b",
           "ENSMUSG00000038067_Csf3",
           "ENSMUSG00000027399_Il1a",
           "ENSMUSG00000024401_Tnf",
           "ENSMUSG00000031712_Il15",
           "ENSMUSG00000026104_Stat1",
           "ENSMUSG00000021025_Nfkbia",
           "ENSMUSG00000004296_Il12b",
           "ENSMUSG00000035385_Ccl2",
           "ENSMUSG00000018899_Irf1",
           "ENSMUSG00000039217_Il18"
           )

pltBM<-merged_BM[pltGene,];pltBM
pltSP<-merged_SP[pltGene,];pltSP

library(tidyverse)

pltAll<-bind_rows(
  pltBM %>% rownames_to_column('Gene') %>% mutate(Organ='BM'),
  pltSP %>% rownames_to_column('Gene') %>% mutate (Organ='SP')
) %>% pivot_longer(
  cols=-c(Gene,Organ),
  names_to=c("Group",'.value'),
  names_pattern="(.*)\\.(.*)"
) %>% mutate(Status=case_when(
  FDR<0.05 & logFC>0 ~'Up',
  FDR<0.05 & logFC<0 ~'Dn',
  TRUE~'NS'
));head(pltAll)
range(pltAll$logFC)

pltAll$Gene<-factor(pltAll$Gene, levels=rev(pltGene))

gp<-ggplot(pltAll, aes(x=Organ, y=Gene)) + 
  geom_point(aes(size=abs(logFC), fill=Status), shape=21, color='black')+
  facet_wrap(~Group)+
  scale_fill_manual(values=c(
  "Up"='red',
  'Dn'='blue',
  "NS"='grey'))+
  scale_size_continuous(range=c(0, 10))+
  theme_minimal()+
  labs(
    title='Different Treatments',
    size='|logFC|',
    fill='FDR<0.05'
  )

pdf(sprintf("%s/BM_SP_DiffTrt_DotPlot.pdf",outputFolder),width=10,height=10)
print(gp)
dev.off()

###########
names(BM)
trtNameAll<-c('Pam','HKML','P','LPS','FLA.ST','FSL',"R837",'R848','ssRNA','ODN')
trtN_order<-c("Pam_vs_NT",
              "HKML_vs_NT",
              "P_vs_NT",
              "LPS_vs_NT",
              "FLA.ST_vs_NT",
              "FSL_vs_NT",
              "R837_vs_NT", 
              "R848_vs_NT", 
              "ssRNA_vs_NT",
              "ODN_vs_NT")
              
BM_All<-lapply(BM[trtN_order], function(x) {x$gene<-row.names(x); x<-x[,c('gene','logFC','FDR')]; return(x)})
BM_All<-lapply(1:length(BM_All), function(x) {colnames(BM_All[[x]])[2:3]<-sprintf("%s.%s", trtNameAll[x], colnames(BM_All[[x]])[2:3]); return(BM_All[[x]])})

SP_All<-lapply(SP[trtN_order], function(x) {x$gene<-row.names(x); x<-x[,c('gene','logFC','FDR')]; return(x)})
SP_All<-lapply(1:length(SP_All), function(x) {colnames(SP_All[[x]])[2:3]<-sprintf("%s.%s", trtNameAll[x], colnames(SP_All[[x]])[2:3]); return(SP_All[[x]])})


# Merge the list by row names using Reduce and merge
merged_BM_All <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), BM_All);head(merged_BM_All)
merged_SP_All <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), SP_All);head(merged_SP_All)

rownames(merged_BM_All) <- merged_BM_All$gene; merged_BM_All$gene <- NULL;head(merged_BM_All)
rownames(merged_SP_All) <- merged_SP_All$gene; merged_SP_All$gene <- NULL;head(merged_SP_All)

selGene<-c("Cxcl10",'Il12a','Il6','Il1b','Csf3','Il1a','Tnf','Il15','Stat1','Nfkbia','Il12b','Ccl2','Irf1','Il18')
lapply(selGene, function(x) rownames(merged_BM_All)[grepl(x, rownames(merged_BM_All))])

pltGene<-c("ENSMUSG00000034855_Cxcl10",
           "ENSMUSG00000027776_Il12a",
           "ENSMUSG00000025746_Il6",
           "ENSMUSG00000027398_Il1b",
           "ENSMUSG00000038067_Csf3",
           "ENSMUSG00000027399_Il1a",
           "ENSMUSG00000024401_Tnf",
           "ENSMUSG00000031712_Il15",
           "ENSMUSG00000026104_Stat1",
           "ENSMUSG00000021025_Nfkbia",
           "ENSMUSG00000004296_Il12b",
           "ENSMUSG00000035385_Ccl2",
           "ENSMUSG00000018899_Irf1",
           "ENSMUSG00000039217_Il18"
)

pltBM_All<-merged_BM_All[pltGene,];pltBM_All
pltSP_All<-merged_SP_All[pltGene,];pltSP_All

library(tidyverse)

pltAll_All<-bind_rows(
  pltBM_All %>% rownames_to_column('Gene') %>% mutate(Organ='BM'),
  pltSP_All %>% rownames_to_column('Gene') %>% mutate (Organ='SP')
) %>% pivot_longer(
  cols=-c(Gene,Organ),
  names_to=c("Group",'.value'),
  names_pattern="(.*)\\.(.*)"
) %>% mutate(Status=case_when(
  FDR<0.05 & logFC>0 ~'Up',
  FDR<0.05 & logFC<0 ~'Dn',
  TRUE~'NS'
));head(pltAll_All)
range(pltAll_All$logFC)

pltAll_All$Gene<-factor(pltAll_All$Gene, levels=rev(pltGene))
pltAll_All$ShortGene<-sapply(strsplit(as.character(pltAll_All$Gene), split="\\_"),function(x) x[2])
pltAll_All$ShortGene<-factor(pltAll_All$ShortGene, 
                             levels=rev(sapply(strsplit(pltGene, split="\\_"), function(x) x[2])))


gp_All<-ggplot(pltAll_All, aes(x=Organ, y=ShortGene)) + 
  geom_point(aes(size=abs(logFC), fill=Status), shape=21, color='black')+
  facet_wrap(~Group, ncol = 5)+
  scale_fill_manual(values=c(
    "Up"='red',
    'Dn'='blue',
    "NS"='grey'))+
  scale_size_continuous(range=c(0, 10))+
  theme_minimal()+
  labs(
    title='Different Treatments',
    size='|logFC|',
    fill='FDR<0.05'
  )

pdf(sprintf("%s/BM_SP_DiffTrt_DotPlot_all.pdf",outputFolder),width=12,height=10)
print(gp_All)
dev.off()


######################
dim(merged_BM_All)
head(merged_BM_All[,1:5])
sum(row.names(merged_BM_All)==row.names(merged_SP_All))

interGene<-intersect(row.names(merged_SP_All), row.names(merged_BM_All));length(interGene)

uniqueTrt<-unique(pltAll_All$Group);uniqueTrt
allScattPlot<-lapply(1:length(uniqueTrt), function(x){
  tempTrt <- uniqueTrt[x]
  temp_BM <- data.frame(
    BM_logFC = merged_BM_All[interGene, sprintf("%s.logFC", tempTrt)],
    SP_logFC = merged_SP_All[interGene, sprintf("%s.logFC", tempTrt)],
    BM_FDR   = merged_BM_All[interGene, sprintf("%s.FDR", tempTrt)],
    SP_FDR   = merged_SP_All[interGene, sprintf("%s.FDR", tempTrt)]
  )
  
  temp_BM <- temp_BM %>%
    mutate(Status = case_when(
      BM_FDR < 0.05 & SP_FDR < 0.05 ~ "Both",
      BM_FDR < 0.05 | SP_FDR < 0.05 ~ "Either",
      TRUE ~ "NS"
    ))
  
  plt <- ggplot(temp_BM, aes(x = BM_logFC, y = SP_logFC, color = Status)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c(
      "Both" = "purple",
      "Either" = "cyan",
      "NS" = "grey80"
    )) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      title = tempTrt,
      x = "BM logFC",
      y = "SP logFC",
      color = "FDR < 0.05"
    ) +
    theme_minimal()
  
  return(plt)
})

require(gridExtra)

pdf(sprintf("%s/ScatterPlot_logFC_supF5.pdf",outputFolder),width=7.87*4/1.2,height=7.87*3/1.2, family="ArialMT")
do.call(grid.arrange, 
        c(allScattPlot, ncol=4,nrow=3))
dev.off()


