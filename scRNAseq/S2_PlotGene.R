library(dplyr)
library(Seurat)
library(Matrix)
library(singleGEO)
library(SoupX)
library(reticulate)
library(anndata)
library(scDblFinder)
library(harmony)
library(cowplot)

minFeature <- 200  
maxFeature <- 7500
minCount <- 400    
maxCount <- 40000
maxMT <- 10       


###################################
CountFolder="./BmSp_Apr2024/count"
outFolder="./BmSp_Apr2024/RDS"
ifelse(!dir.exists(outFolder), dir.create(outFolder), FALSE)
tmpOutput="./temp"

(AllSample<-list.files(CountFolder,pattern="-"))
pData<-read.delim(sprintf("%s/pCount.txt",CountFolder),header=F,sep="\t")
sum(pData$V1 %in% AllSample)

objName="BmSp_Apr2024"
numRef=10

###################################
saveObj<-readRDS(file=sprintf("%s/BM_SP.saveObj.RDS",outFolder))
seu<-saveObj$seu
harF<-saveObj$harF

#########################################################
harF$grp<-paste(harF$Organ,harF$Stage,sep="_")
dp<-DimPlot(harF,group.by="RNA_snn_res.0.2",label=T,split.by="grp",ncol=2,raster=F)+NoLegend()

pdf(sprintf("%s/DimPlot_may22.pdf",tmpOutput),width=8,height=8, family="ArialMT")
print(dp)
dev.off()

dp<-DimPlot(harF,group.by="RNA_snn_res.0.2",label=F,split.by="grp",ncol=2,raster=F,cols=c('grey50','blue'))+NoLegend()

pdf(sprintf("%s/DimPlot_may22_color.pdf",tmpOutput),width=8,height=8, family="ArialMT")
print(dp)
dev.off()


dp0<-DimPlot(harF,group.by="RNA_snn_res.0.2",label=F,raster=F)+NoLegend()

pdf(sprintf("%s/DimPlot_Sep2024.pdf",tmpOutput),width=6,height=6, family="ArialMT")
print(dp0)
dev.off()

(dp01<-DimPlot(harF,group.by="RNA_snn_res.0.2",label=F,raster=F,cols=c('grey50','blue'))+NoLegend())

pdf(sprintf("%s/DimPlot_Sep2024_changeColor.pdf",tmpOutput),width=6,height=6, family="ArialMT")
print(dp01)
dev.off()



(FP_Paper<-FeaturePlot(harF,
                       features=c('Mki67','Top2a','Il1b'),
                       ncol=2,
                       order=T,
                       col=c('grey','red')))

pdf(sprintf("%s/CM_IL1B_mki67.pdf",tmpOutput),width=10,height=8)
print(FP_Paper)
dev.off()


#####################################
library(msigdbr)
library(fgsea)
library(data.table)
library(UCell)
unique(as.data.frame(msigdbr(species = "Mus musculus", category ="C5"))$gs_subcat)
Hall<-geneSets <- msigdbr(species = 'Mus musculus', category = 'H');head(Hall)
unique(Hall$gs_name)
Hall_Inf<-Hall[Hall$gs_name=="HALLMARK_INFLAMMATORY_RESPONSE",];head(Hall_Inf)
Hall_G2M<-Hall[Hall$gs_name=="HALLMARK_G2M_CHECKPOINT",];head(Hall_G2M)

##################################################################
Inflammaotry <- list(Inflammatory = Hall_Inf$gene_symbol, G2M=Hall_G2M$gene_symbol)
harF<-AddModuleScore_UCell(harF, features =Inflammaotry)

unique(harF@meta.data$Organ)
wilcox.test(Inflammatory_UCell~Organ,harF@meta.data)
wilcox.test(G2M_UCell~Organ,harF@meta.data)

harFMeta<-harF@meta.data
unique(harFMeta$grp)
wilcox.test(Inflammatory_UCell~RNA_snn_res.0.2,harFMeta[harFMeta$grp=="SP_Resting",])
wilcox.test(Inflammatory_UCell~RNA_snn_res.0.2,harFMeta[harFMeta$grp=="SP_Day7",])

fpI<-FeaturePlot(harF,features=c("Inflammatory_UCell"),col=c('grey','red'))+ggtitle('INFLAMMATORY_RESPONSE_SCORE')
(vpI<-VlnPlot(harF,features=c("Inflammatory_UCell"), col=c('grey','blue')))+ggtitle('INFLAMMATORY_RESPONSE_SCORE')
(vpG<-VlnPlot(harF,features=c("G2M_UCell"), col=c('grey','blue')))+ggtitle('G2M_CHECKPOINT')


pdf(sprintf("%s/Inflam_G2M_Violin.pdf",tmpOutput),width=8,height=6,family="ArialMT")
print(vpI+ggtitle('INFLAMMATORY_RESPONSE_SCORE'))
print(vpG+ggtitle('G2M_CHECKPOINT'))
dev.off()

table(harF$RNA_snn_res.0.2,harF$ID)
(vpI__2<-VlnPlot(harF,features=c("Inflammatory_UCell"), split.by="ID"))+ggtitle('INFLAMMATORY_RESPONSE_SCORE')
(vpG__2<-VlnPlot(harF,features=c("G2M_UCell"), col=c('grey','blue')))+ggtitle('G2M_CHECKPOINT')


(vpI2<-VlnPlot(harF,features=c("Inflammatory_UCell"),split.by='grp'))+ggtitle('INFLAMMATORY_RESPONSE_SCORE')

pdf(sprintf("%s/InflamScore_plot1.pdf",tmpOutput),width=8,height=6,family="ArialMT")
print(fpI+ggtitle('INFLAMMATORY_RESPONSE_SCORE'))
dev.off()

pdf(sprintf("%s/InflamScore_Violin.pdf",tmpOutput),width=8,height=6,family="ArialMT")
print(vpI+ggtitle('INFLAMMATORY_RESPONSE_SCORE'))
print(vpI2+ggtitle('INFLAMMATORY_RESPONSE_SCORE'))
dev.off()

harF2<-harF
Idents(harF2)<-"grp"
(vpI3<-VlnPlot(harF2,features=c("Inflammatory_UCell"),split.by='RNA_snn_res.0.2'))+ggtitle('INFLAMMATORY_RESPONSE_SCORE')

Idents(harF2)<-"RNA_snn_res.0.2"
(vpI4<-VlnPlot(harF2,features=c("Inflammatory_UCell"),split.by='grp'))+ggtitle('INFLAMMATORY_RESPONSE_SCORE')

pdf(sprintf("%s/InflamScore_Violin_split.pdf",tmpOutput),width=8,height=6,family="ArialMT")
print(vpI3+ggtitle('INFLAMMATORY_RESPONSE_SCORE'))
print(vpI4+ggtitle('INFLAMMATORY_RESPONSE_SCORE'))
dev.off()

##############################################################
DimPlot(harF,group.by="CellType",label=T)
colnames(harF@meta.data)[13]<-"InflammatoryScore"
colnames(harF@meta.data)[14]<-"G2MScore"

saveRDS(harF,file=sprintf("%s/Seu_MouseCM_SP.RDS",tmpOutput))
harF<-readRDS(file=sprintf("%s/Seu_MouseCM_SP.RDS",tmpOutput))


seu@meta.data<-seu@meta.data[,!colnames(seu@meta.data) %in% 'CellType']
colnames(seu@meta.data)[11]<-"FinalCellType"
DimPlot(seu,group.by="FinalCellType",label=T)
saveRDS(seu,file=sprintf("%s/allMouseSpCells.RDS",tmpOutput))


ridgeP<-RidgePlot(harF, 
          features = c("InflammatoryScore","G2MScore"), 
          ncol = 2, 
          col=c('grey','blue'))


pdf(sprintf("%s/DistScore.pdf",tmpOutput),width=12,height=6,family="ArialMT")
print(ridgeP)
dev.off()



(vpI<-VlnPlot(harF,features=c("InflammatoryScore"), col=c('grey','blue'), pt.size=0))+ggtitle('INFLAMMATORY_RESPONSE_SCORE')
(vpG<-VlnPlot(harF,features=c("G2MScore"), col=c('grey','blue'), pt.size=0))+ggtitle('G2M_CHECKPOINT')


pdf(sprintf("%s/Paper__Inflam_G2M_Violin_noDot.pdf",tmpOutput),width=8,height=6,family="ArialMT")
print(vpI+ggtitle('INFLAMMATORY_RESPONSE_SCORE'))
print(vpG+ggtitle('G2M_CHECKPOINT'))
dev.off()


######### pseudobluk ##
harF$NewCT<-ifelse(harF$RNA_snn_res.0.2=="0",'BM','SP')
table(harF$ID, harF$Organ)
table(harF$ID, harF$Stage)

(cellProp<-as.data.frame(prop.table(table(harF$ID, harF$NewCT),margin=1)))
cellProp$Stage<-dplyr::case_when(
  cellProp$Var1 %in% c("AB-24-11","AB-24-12","AB-24-16",'AB-24-17')~"D7",
  .default='Resting'
)

cellProp$Var2<-factor(cellProp$Var2, levels=c('BM','SP'))
cellProp$Stage<-factor(cellProp$Stage, levels=c('Resting','D7'))
cellProp<-cellProp[order(cellProp$Stage),]
cellProp$Var1<-factor(cellProp$Var1, levels=c("AB-24-8","AB-24-9","AB-24-10",
                                              "AB-24-14","AB-24-13","AB-24-15",
                                              "AB-24-11","AB-24-12",
                                              "AB-24-16","AB-24-17"))
p<-ggplot(cellProp, 
       aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("BM" = "cyan", "SP" = "magenta"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(title = "",
       x = "ID",
       y = "Fraction")

pdf(sprintf("%s/Paper__barplot_fraction_BMSP.pdf",tmpOutput),width=6,height=6,family="ArialMT")
print(p)
dev.off()



head(harF);table(harF$ID, harF$Organ)
pseudo_counts <- AggregateExpression(harF,
                                     group.by =c("ID"), #c("ID", "NewCT"),
                                     assays = "RNA",
                                     return.seurat=FALSE)$RNA
head(pseudo_counts)

table(harF$ID, harF$NewCT)

pseudo<-pseudo_counts
### gsva score were calcuated in python; can not install gsva pacakge in R
write.csv(pseudo, file=sprintf("%s/pseudo_BMSP.csv",tmpOutput))





harF_twoDis<-subset(harF, ((Organ %in% 'BM') & (NewCT %in% 'BM')) | ((Organ %in% 'SP') & (NewCT %in% 'SP')))
DimPlot(harF_twoDis, group.by='NewCT', split.by='Organ')
DimPlot(harF, group.by='NewCT', split.by='Organ')

pseudo_counts_TwoDis <- AggregateExpression(harF_twoDis,
                                     group.by =c("ID"), #c("ID", "NewCT"),
                                     assays = "RNA",
                                     return.seurat=FALSE)$RNA
head(pseudo_counts_TwoDis)

table(harF_twoDis$ID, harF_twoDis$NewCT)

pseudo_twoDis<-pseudo_counts_TwoDis
### gsva score were calcuated in python; can not install gsva pacakge in R
write.csv(pseudo_twoDis, file=sprintf("%s/pseudo_BMSP_twoDis.csv",tmpOutput))



