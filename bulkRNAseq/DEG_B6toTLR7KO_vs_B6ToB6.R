
inputFolder<-"./Dec2024/mouse/RawCount"
outputFolder<-"./DEG"
SaveFolder<-"./Dec2024/mouse/RDS"
if (!dir.exists(SaveFolder)) {dir.create(SaveFolder, recursive = TRUE)}
source("~/code/sourceRcode/bulkRNAseq_sourceCode.R")


#############################################################################################################
############                                                                                       ##########
############                           read the data                                               ##########
############                                                                                       ##########
#############################################################################################################
list.files(inputFolder)
dat<-read.delim(sprintf("%s/AllGeneExpCount.txt",inputFolder),header=T,sep="\t",row.names=1,stringsAsFactors = F)
colnames(dat)<-sapply(strsplit(colnames(dat),split="\\."),function(x) sprintf("%s_%s",x[8],x[9]))
#############################################################################################################
############                                                                                       ##########
############                           protein coding only                                         ##########
############                                                                                       ##########
#############################################################################################################
anoFile<-read.delim("~/code/sourceRcode/sourceData/mart_export_grcm38.p6.txt",header=T,sep="\t",stringsAsFactors=FALSE);head(anoFile)
anoFile$ID<-paste(anoFile[,1],anoFile[,3],sep="_")
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]#;#twoGene<-c("ENSG00000017427_IGF1","ENSG00000100985_MMP9")

dat_coding<-dat[sapply(strsplit(row.names(dat),split="\\_"),function(x) x[1]) %in% codingGene$Gene.stable.ID,];dim(dat_coding)#21832
#############################################################################################################
############                                                                                       ##########
#############################################################################################################
pdata<-data.frame(ID=colnames(dat_coding),
                  trt=rep(c("Tlr7KOToB6",'B6ToTlr7KO','B6ToB6','Tlr7KOToTlr7KO'),each=3));pdata
pdata$col<-c("magenta","cyan","tomato",'green')[factor(pdata$trt)]
pdata


#############################################################################################################
############                                                                                       ##########
#############################################################################################################
pdata<-pdata[pdata$trt %in% c('B6ToTlr7KO','B6ToB6'),]; pdata
pdata$trt<-factor(pdata$trt, levels=c('B6ToB6','B6ToTlr7KO'));pdata
dat_input<-dat_coding[,paste(pdata$ID)];dim(dat_input) #21832*10
dat_input<-dat_input[rowSums(dat_input)>0,];dim(dat_input)#15801
sum(colnames(dat_input)==pdata$ID)
boxplot(log1p(dat_input),las=2)
##############
library(edgeR)
edgeInput1<-DGEList(dat_input,group=as.factor(pdata$trt))
edgeInput2 <- edgeInput1[rowSums(1e+06*edgeInput1$counts/expandAsMatrix(edgeInput1$samples$lib.size,dim(edgeInput1))>=1)>=ncol(edgeInput1)*0.1, ]
edgeInput2 <- calcNormFactors(edgeInput2)
edgeInput2$samples


pdf(sprintf("%s/MDS_B6ToTlr7Ko_vs_B6ToB6.pdf",outputFolder),  width=6,height=6)
plotMDS(cpm(edgeInput2,prior.count=1,log=TRUE), col=pdata$col,pch=19,cex=2,xlab="Dim1",ylab="Dim2")
dev.off()


testCPM<-cpm(edgeInput2,prior.count=1,log=TRUE);head(testCPM)
plotMDS(testCPM, col=pdata$col,pch=19,cex=2,xlab="Dim1",ylab="Dim2")

plotMDS(cpm(edgeInput2,prior.count=1,log=TRUE), col=pdata$col,pch=19,cex=2,xlab="Dim1",ylab="Dim2")

#############################################################################################################
############                                                                                       ##########
#############################################################################################################
pTrt<-pdata
design.mat <- model.matrix(~pTrt$trt)
library(stringr)
(colnames(design.mat) <-c("(Intercept)","B6ToTlr7KO"))#
rownames(design.mat) <- colnames(edgeInput2)
cpm_data<-cpm(edgeInput2,prior.count=1,log=TRUE)
write.csv(cpm_data,file=sprintf("%s/cpm_B6ToTlr7KO_vs_B6ToB6__Apr2026.csv",outputFolder))


edgeInput3 <- estimateDisp(edgeInput2, design.mat, robust=TRUE) #
fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE)

B6ToTlr7KO_vs_B6ToB6<-glmQLFTest(fit, coef=2)


DE.ALL<-list(B6ToTlr7KO_vs_B6ToB6=B6ToTlr7KO_vs_B6ToB6)

saveRDS(DE.ALL, file=sprintf("%s/MouseTLR.DE.B6ToTlr7KO_vs_B6ToB6.RDS",SaveFolder))


De.table<-lapply(DE.ALL,function(xy) {
    temp<-xy$table
    temp$FDR<-p.adjust(temp$PValue,method="BH")
    temp$stat_gsea<-zscoreT(sign(temp$logFC)*sqrt(temp$F),df=xy$df.total);return(temp)})
lapply(De.table,function(xy) head(xy))
saveRDS(De.table,file=sprintf("%s/TLF_DeTable.B6ToTlr7KO_vs_B6ToB6.rds",outputFolder))
  
De.table<-readRDS(file=sprintf("%s/TLF_DeTable.B6ToTlr7KO_vs_B6ToB6.rds",outputFolder))
FDRCut=0.05;logCut=log2(1.5)
singCount<-do.call("rbind",lapply(De.table,function(xy) {
    tSig<-sum(xy$FDR<FDRCut & abs(xy$logFC)>logCut)
    uSig<-sum(xy$FDR<FDRCut & xy$logFC>logCut)
    dSig<-sum(xy$FDR<FDRCut & xy$logFC<(-logCut))
    out<-c(tSig,uSig,dSig)}))
colnames(singCount)<-c("TotalSig","UpSig","DnSig");singCount
write.table(data.frame(Comparison=row.names(singCount),singCount),file=sprintf("%s/TLR_singCount.B6ToTlr7KO_vs_B6ToB6.xls",outputFolder),col.names=T,row.names=F,sep="\t")
  

###write stat out
lapply(1:length(De.table),function(xy){
  compName<-names(De.table)[xy]
  temp<-De.table[[xy]]
  write.table(data.frame(Gene=row.names(temp),temp),file=sprintf("%s/tlr_StatOut.%s.xls",outputFolder,compName),row.names=F,col.names=T,sep="\t")
})

#############################################################################################################
length(De.table)

sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/TLR_MA_B6ToTlr7KO_vs_B6ToB6.tiff",outputFolder), res=600, width=10,height=10,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/TLR_MA_B6ToTlr7KO_vs_B6ToB6.pdf",outputFolder),width=10,height=10)}
  par(mfrow=c(2,2))
  lapply(1:length(De.table),function(xy){
    compName<-names(De.table)[xy]
    temp<-De.table[[xy]]
    temp$col<-ifelse(temp$FDR<0.05, 'red','black')
    plot(temp$logCPM, temp$logFC,ylab="Log2 FC",xlab="logCPM", col=temp$col)#,main=compName)
    })
  dev.off()
})

#############################################################################################################
sigGene<-lapply(De.table,function(xy){xy[xy$FDR< FDRCut & abs(xy$logFC)>logCut,]});head(sigGene[[1]])
sigGene<-lapply(sigGene,function(xy) {
  xy$gene<-sapply(strsplit(row.names(xy),split="\\_"),function(y) y[2])
  return(xy)})



#####################################################

library(msigdbr)
library(fgsea)
library(data.table)

gsea_Ana_edgeR<-function(species = "Homo sapiens",
                         category = "C2",
                         subcategory = "CP:REACTOME",
                         DEinput_df=DEinput_df, ##data.frame; should have stat_gsea colum
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


DE.ALL2<-DE.ALL

De.table2<-lapply(DE.ALL2,function(xy) {
  temp<-xy$table
  temp$FDR<-p.adjust(temp$PValue,method="BH")
  temp$stat_gsea<-zscoreT(sign(temp$logFC)*sqrt(temp$F),df=xy$df.total);return(temp)})


hall_list<-lapply(De.table2,function(xbp){
  H_gsea<-gsea_Ana_edgeR(species = "Mus musculus",
                         category = "H",
                         subcategory = NULL,
                         DEinput_df=xbp,
                         cut_padj=2)
  return(H_gsea)
});names(hall_list)<-names(De.table2)

H_B6ToTLRKO<-hall_list$B6ToTlr7KO_vs_B6ToB6


H_B6ToTLRKO$eaRes
fwrite(H_B6ToTLRKO$eaRes, file=sprintf("%s/Hall_B6ToTlr7KO_vs_B6ToB6.xls",outputFolder), sep="\t", sep2=c("", " ", ""))

# pathway         pval         padj         ES        NES nMoreExtreme  size
# <char>        <num>        <num>      <num>      <num>        <num> <int>
#   1:                        HALLMARK_DNA_REPAIR 3.886212e-05 0.0002069793  0.4721943  2.6411713            0   142
# 2:                       HALLMARK_E2F_TARGETS 3.863092e-05 0.0002069793  0.5641253  3.3196330            0   195
# 3:             HALLMARK_FATTY_ACID_METABOLISM 3.881536e-05 0.0002069793  0.4711265  2.5555384            0   119
# 4:                    HALLMARK_G2M_CHECKPOINT 3.867724e-05 0.0002069793  0.3340913  1.9497691            0   185
# 5:                   HALLMARK_MITOTIC_SPINDLE 4.139587e-05 0.0002069793 -0.3392894 -2.0003788            0   187
# 6:                  HALLMARK_MTORC1_SIGNALING 3.865930e-05 0.0002069793  0.4446464  2.5931267            0   184
# 7:                    HALLMARK_MYC_TARGETS_V1 3.862346e-05 0.0002069793  0.6600773  3.9011393            0   200
# 8:                    HALLMARK_MYC_TARGETS_V2 3.983746e-05 0.0002069793  0.4813550  2.2498731            0    58
# 9:         HALLMARK_OXIDATIVE_PHOSPHORYLATION 3.863092e-05 0.0002069793  0.6483690  3.8153703            0   195
# 10:         HALLMARK_UNFOLDED_PROTEIN_RESPONSE 3.906403e-05 0.0002069793  0.3926115  2.0798582            0   105
# 11:                      HALLMARK_ADIPOGENESIS 7.724692e-05 0.0003423579  0.3230009  1.8607976            1   170
# 12:                    HALLMARK_UV_RESPONSE_DN 8.216589e-05 0.0003423579 -0.3707664 -1.9520454            1    98
# 13:             HALLMARK_XENOBIOTIC_METABOLISM 1.632399e-03 0.0062784586  0.2993469  1.6481681           41   130
# 14:           HALLMARK_ESTROGEN_RESPONSE_EARLY 1.937186e-03 0.0069185204 -0.3085607 -1.6652806           46   112
# 15:   HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY 2.866470e-03 0.0095549009  0.3917758  1.7621630           71    49
# 16:                HALLMARK_HEDGEHOG_SIGNALING 3.493947e-03 0.0102763133 -0.5320732 -1.8842270           85    20
# 17:                        HALLMARK_PEROXISOME 3.480641e-03 0.0102763133  0.3314500  1.6621898           88    79
# 18:        HALLMARK_WNT_BETA_CATENIN_SIGNALING 5.441624e-03 0.0151156232 -0.4729641 -1.7893485          133    25
# 19:                HALLMARK_TGF_BETA_SIGNALING 7.653575e-03 0.0201409869 -0.3795079 -1.6926376          189    46
# 20:                        HALLMARK_GLYCOLYSIS 8.190358e-03 0.0204758947  0.2651816  1.4954689          210   149
# 21:                 HALLMARK_PROTEIN_SECRETION 1.820388e-02 0.0433425798  0.2918802  1.4912898          464    87
# 22:                    HALLMARK_APICAL_SURFACE 2.281400e-02 0.0498856013 -0.4768198 -1.6402678          561    18
# 23:                    HALLMARK_UV_RESPONSE_UP 2.294738e-02 0.0498856013  0.2604555  1.4231506          589   124
# 24:                   HALLMARK_NOTCH_SIGNALING 3.025381e-02 0.0630287648 -0.4140688 -1.5665319          744    25
# 25:                 HALLMARK_KRAS_SIGNALING_DN 4.870600e-02 0.0974119916 -0.3142394 -1.4240848         1211    49
# 26:                         HALLMARK_APOPTOSIS 7.207382e-02 0.1371938734  0.2343498  1.2874960         1858   128
# 27:         HALLMARK_INTERFERON_ALPHA_RESPONSE 7.408469e-02 0.1371938734 -0.2535226 -1.3105801         1808    89
# 28:                   HALLMARK_APICAL_JUNCTION 1.142281e-01 0.2039787985 -0.2304428 -1.2337563         2776   107
# 29:               HALLMARK_PANCREAS_BETA_CELLS 1.422958e-01 0.2453375598  0.4296440  1.3294260         3584    13
# 30:                   HALLMARK_HEME_METABOLISM 1.504754e-01 0.2507923384 -0.2076090 -1.1760506         3639   147
# 31:           HALLMARK_TNFA_SIGNALING_VIA_NFKB 3.276490e-01 0.5284661397 -0.1835472 -1.0562581         7915   161
# 32:           HALLMARK_PI3K_AKT_MTOR_SIGNALING 3.931145e-01 0.5956279679 -0.1981832 -1.0289977         9602    91
# 33:                   HALLMARK_SPERMATOGENESIS 3.906754e-01 0.5956279679 -0.2165377 -1.0342231         9694    61
# 34:            HALLMARK_ESTROGEN_RESPONSE_LATE 4.569632e-01 0.6720046407  0.1837373  0.9982631        11769   120
# 35:               HALLMARK_ALLOGRAFT_REJECTION 5.172227e-01 0.7388895046  0.1716654  0.9724136        13333   153
# 36:                 HALLMARK_ANDROGEN_RESPONSE 5.622309e-01 0.7808762217 -0.1857063 -0.9458551        13709    83
# 37:                      HALLMARK_ANGIOGENESIS 6.449840e-01 0.8646053386  0.2426425  0.8702370        16310    21
# 38:              HALLMARK_BILE_ACID_METABOLISM 7.435606e-01 0.8646053386  0.1727819  0.8493395        18965    72
# 39:           HALLMARK_CHOLESTEROL_HOMEOSTASIS 7.205590e-01 0.8646053386  0.1809450  0.8558450        18147    61
# 40:                        HALLMARK_COMPLEMENT 6.720353e-01 0.8646053386 -0.1610101 -0.9100276        16314   145
# 41: HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION 6.894245e-01 0.8646053386  0.1740722  0.8834298        17620    84
# 42:                        HALLMARK_MYOGENESIS 7.357511e-01 0.8646053386 -0.1681318 -0.8648652        17994    87
# 43:                       HALLMARK_P53_PATHWAY 7.360993e-01 0.8646053386  0.1549941  0.8830124        19036   158
# 44:                       HALLMARK_COAGULATION 9.543211e-01 0.9830730600  0.1408237  0.6963002        24359    74
# 45:                           HALLMARK_HYPOXIA 9.521200e-01 0.9830730600 -0.1343316 -0.7503868        23106   135
# 46:               HALLMARK_IL2_STAT5_SIGNALING 9.053343e-01 0.9830730600 -0.1397181 -0.8018164        21842   159
# 47:           HALLMARK_IL6_JAK_STAT3_SIGNALING 8.951593e-01 0.9830730600 -0.1548304 -0.7542204        22079    67
# 48:         HALLMARK_INTERFERON_GAMMA_RESPONSE 9.634116e-01 0.9830730600 -0.1296556 -0.7595238        23223   180
# 49:                 HALLMARK_KRAS_SIGNALING_UP 9.622766e-01 0.9830730600 -0.1335000 -0.7264242        23314   117
# 50:             HALLMARK_INFLAMMATORY_RESPONSE 9.926246e-01 0.9926246395 -0.1197727 -0.6745076        24090   142


library(ggplot2)
(plt<-plotEnrichment(H_B6ToTLRKO$m_list[["HALLMARK_G2M_CHECKPOINT"]], H_B6ToTLRKO$stats) +
  labs(title='G2M_CHECKPOINT'))

pdf(sprintf("%s/B6ToTlr7KO_vs_B6ToB6__G2M.pdf",outputFolder),width=7.87,height=7.87, family="ArialMT")
print(plt)
dev.off()


(plt2<-plotEnrichment(H_B6ToTLRKO$m_list[["HALLMARK_INFLAMMATORY_RESPONSE"]], H_B6ToTLRKO$stats) +
  labs(title='INFLAMMATORY_RESPONSE'))




testTab<-De.table2$B6ToTlr7KO_vs_B6ToB6;head(testTab)
testTab[grep("Il1b",row.names(testTab)),]
testTab[grep("Top2a",row.names(testTab)),]
testTab[grep("Mki67",row.names(testTab)),]


dim(cpm_data)
cpm_data[grep("Il1b",row.names(cpm_data)),]


