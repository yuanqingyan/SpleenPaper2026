FolderDay1to5<-"./BulkRNAseq/TO/Oct2022/RawCount"

outputFolder<-"./RNAseq/Com2022/"
SaveFolder<-"./BulkRNAseq/TO/Oct2022/RDS"
source("~/code/sourceRcode/bulkRNAseq_sourceCode.R")
#############################################################################################################
############                                                                                       ##########
############                           read the data                                               ##########
############                                                                                       ##########
#############################################################################################################
Day15<-read.delim(sprintf("%s/AllGeneExpCount.txt",FolderDay1to5),header=T,sep="\t",row.names=1,stringsAsFactors = F)
colnames(Day15)<-sapply(strsplit(colnames(Day15),split="\\."),function(x) x[8])
colnames(Day15)<-gsub("_",".",paste("Oct2022",sapply(strsplit(colnames(Day15),split="_S"),function(x) x[1]),sep="."))

dat_com<-Day15;dim(dat_com)
#############################################################################################################
############                                                                                       ##########
############                           protein coding only                                         ##########
############                                                                                       ##########
#############################################################################################################
anoFile<-read.delim("~/code/sourceRcode/sourceData/mart_export_grcm38.p6.txt",header=T,sep="\t",stringsAsFactors=FALSE);head(anoFile)
anoFile$ID<-paste(anoFile[,1],anoFile[,3],sep="_")
codingGene<-anoFile[anoFile$Gene.type=="protein_coding",]

dat_coding<-dat_com[sapply(strsplit(row.names(dat_com),split="\\_"),function(x) x[1]) %in% codingGene$Gene.stable.ID,];dim(dat_coding)
#############################################################################################################
############                                                                                       ##########
############.                          read pheno data                                             ##########
############                                                                                       ##########
#############################################################################################################
pdata<-data.frame(ID=colnames(dat_coding),
                  batch=sapply(strsplit(colnames(dat_coding),split="\\."),function(x) x[1]),
                  day=sapply(strsplit(colnames(dat_coding),split="\\."),function(x) x[2]));pdata
pdata$day <- factor(pdata$day, levels=c("WT","1D","2D","3D","5D"))
pdata$col<-c("black","green","purple","blue","red")[factor(pdata$day)]

#############################################################################################################
############                                                                                       ##########
############.                         prepare input dat                                            ##########
############                                                                                       ##########
#############################################################################################################
dat_input<-dat_coding[,paste(pdata$ID)];dim(dat_input) #
dat_input<-dat_input[rowSums(dat_input)>0,];dim(dat_input)#
sum(colnames(dat_input)==pdata$ID)
boxplot(log1p(dat_input),las=2)
##############
library(edgeR)
edgeInput1<-DGEList(dat_input,group=as.factor(pdata$day))
edgeInput2 <- edgeInput1[rowSums(1e+06*edgeInput1$counts/expandAsMatrix(edgeInput1$samples$lib.size,dim(edgeInput1))>=1)>=ncol(dat_input)*0.1, ]
edgeInput2 <- calcNormFactors(edgeInput2)
edgeInput2$samples

#############################################################################################################
#############################################################################################################
design.mat <- model.matrix(~pdata$day)
colnames(design.mat) <- sapply(strsplit(colnames(design.mat),split="day"),function(x) x[2])
colnames(design.mat)[c(1)]<- c("(Intercept)")
rownames(design.mat) <- colnames(edgeInput2)
design.mat

#############################################################################################################
cpm_data<-cpm(edgeInput2,prior.count=1,log=TRUE)
listCPM<-list(cpm=cpm_data,edgeInput2=edgeInput2)
saveRDS(listCPM,file=sprintf("%s/listCPM_5Days.rds",SaveFolder))










