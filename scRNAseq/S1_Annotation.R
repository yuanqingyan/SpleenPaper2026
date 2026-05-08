
library(dplyr)
library(Seurat)
library(Matrix)
library(singleGEO)
library(SoupX)
library(reticulate)
library(anndata)
library(sceasy)
library(scDblFinder)


minFeature <- 200  
maxFeature <- 7500
minCount <- 400    #
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
seu<-readRDS(file=sprintf("%s/%s_Soupx_beforeAnnotation.rds",outFolder,objName))
unique(seu$orig.ident)
seu$Organ<-case_when(
  seu$orig.ident %in% c("AB-24-8","AB-24-9","AB-24-10","AB-24-11","AB-24-12")~"BM",
  seu$orig.ident %in% c("AB-24-13","AB-24-14","AB-24-15","AB-24-16","AB-24-17")~"SP")

seu$Stage<-case_when(
  seu$orig.ident %in% c("AB-24-8","AB-24-9","AB-24-10","AB-24-13","AB-24-14","AB-24-15")~"Resting",
  seu$orig.ident %in% c("AB-24-11","AB-24-12","AB-24-16","AB-24-17")~"Day7")

DimPlot(seu,group.by="Organ",label=T)

###################################
DimPlot(seu,group.by="ID")
DefaultAssay(seu)<-"RNA";Idents(seu)<-"integrated_snn_res.1.5"
DimPlot(seu,label=T)+NoLegend()
FeaturePlot(seu,features=c("Ptprc",'Ms4a1','Ccr7','Igkc'),col=c('grey','red'))
FeaturePlot(seu,features=c("Mrc1"),col=c('grey','red'))

Marker<-FindAllMarkers(seu, slot='data', test.use='wilcox', only.pos=T)
Marker %>% group_by(cluster) %>% top_n(5, avg_log2FC) %>% as.data.frame()



################
DimPlot(seu,group.by="integrated_snn_res.1.5",label=T)+NoLegend()

Marker2<-Marker[order(-Marker$avg_log2FC),]
head(Marker2[Marker2$cluster=="18",],20)

FeaturePlot(seu,features=c("Ptprc"),col=c("grey","red"),raster=F) #For immune cells
FeaturePlot(seu,features=c("percent.mt","nCount_RNA"),col=c("grey","red"),raster=F) #For immune cells
FeaturePlot(seu,features=c("Ly6c2","F13a1","S100a4"),col=c("grey","red"),raster=F) #For cm

FeaturePlot(seu,features=c("Alas2","Hba.a1","Car1"),col=c("grey","red"),raster=F) #RBC-----cluster 17; 10;15
FeaturePlot(seu,features=c("Myb","Cdk6"),col=c("grey","red"),raster=F) #HSC-----7
FeaturePlot(seu,features=c("Nkg7","Ctsg","Mpo"),col=c("grey","red"),raster=F) #HPC
FeaturePlot(seu,features=c("Mki67","Top2a"),col=c("grey","red"),raster=F) #
FeaturePlot(seu,features=c("Igkc","Cd79a"),col=c("grey","red"),raster=F) #Plasma/B-----cluster 21
FeaturePlot(seu,features=c("Ms4a2","Mcpt8"),col=c("grey","red"),raster=F) #Mast cell-----cluster 20
FeaturePlot(seu,features=c("Mmp9","S100a9","Ly6g"),col=c("grey","red"),raster=F) #Neutrophils---cluster 2;19
FeaturePlot(seu,features=c("Cd74","H2.Aa","H2.Ab1","Cd209a"),col=c("grey","red"),raster=F) #Dc---14


#####cluster 18
FeaturePlot(seu,features=c("Crip1","Gng2"),col=c('grey','red'))

#### cluster 3
FeaturePlot(seu,features=c("Ifi27l2a","Csf1r"),col=c('grey','red'))

#### cluster 12
FeaturePlot(seu,features=c("Ccl6","Ccr2"),col=c('grey','red'))

#### cluster 4
FeaturePlot(seu,features=c("Ccl6","Apoe","Ly6c2"),col=c('grey','red'))



########
seu$CellType<-dplyr::case_when(
  seu$integrated_snn_res.1.5 %in% c("0","1","2","5","13","19")~"Neutrophils",
  seu$integrated_snn_res.1.5 %in% c("21")~"Plasma",
  seu$integrated_snn_res.1.5 %in% c("20")~"Mast",
  seu$integrated_snn_res.1.5 %in% c("10","15","17")~"RBC",
  seu$integrated_snn_res.1.5 %in% c("14")~"DC",
  seu$integrated_snn_res.1.5 %in% c("7","11","16")~"HSCs",
  seu$integrated_snn_res.1.5 %in% c("3","4","12","6","8","9","18")~"CM",
  .default='Unknown'
);table(seu$CellType)

saveRDS(seu,file=sprintf("%s/BM_SP.Annotated.RDS",outFolder))
########################################
p1<-DimPlot(seu,group.by="integrated_snn_res.1.5",label=T,repel=T)+NoLegend()
DefaultAssay(seu)<-"RNA"

tiff(sprintf("%s/BM_SP_Dim.tiff",tmpOutput),res=300,width=7,height=7,compression = "lzw",unit="in")
print(p1)
dev.off()

tiff(sprintf("%s/FP_immune.tiff",tmpOutput),res=300,width=7,height=7,compression = "lzw",unit="in")
print(FeaturePlot(seu,features=c("Ptprc"),col=c("grey","red"),raster=F)) #For immune cells
dev.off()

tiff(sprintf("%s/FP_RBCs.tiff",tmpOutput),res=300,width=12,height=10,compression = "lzw",unit="in")
fp<-FeaturePlot(seu,features=c("Alas2","Hba.a1","Car1"),col=c("grey","red"),raster=F,ncol=2) #RBC-----cluster 17; 10;15
print(fp)
dev.off()

tiff(sprintf("%s/FP_HSCs.tiff",tmpOutput),res=300,width=12,height=10,compression = "lzw",unit="in")
fp<-FeaturePlot(seu,features=c("Mpo","Cdk6","Nkg7","Elane"),col=c("grey","red"),raster=F,ncol=2) #HSC-----7
print(fp)
dev.off()

tiff(sprintf("%s/FP_proliferative.tiff",tmpOutput),res=300,width=12,height=5,compression = "lzw",unit="in")
fp<-FeaturePlot(seu,features=c("Mki67","Top2a"),col=c("grey","red"),raster=F) #
print(fp)
dev.off()

tiff(sprintf("%s/FP_plasma.tiff",tmpOutput),res=300,width=12,height=5,compression = "lzw",unit="in")
fp<-FeaturePlot(seu,features=c("Igkc","Cd79a"),col=c("grey","red"),raster=F) #Plasma/B-----cluster 21
print(fp)
dev.off()

tiff(sprintf("%s/FP_Mast.tiff",tmpOutput),res=300,width=12,height=5,compression = "lzw",unit="in")
fp<-FeaturePlot(seu,features=c("Ms4a2","Mcpt8"),col=c("grey","red"),raster=F) #Mast cell-----cluster 20
print(fp)
dev.off()

tiff(sprintf("%s/FP_Neutrophils.tiff",tmpOutput),res=300,width=12,height=10,compression = "lzw",unit="in")
fp<-FeaturePlot(seu,features=c("Mmp9","S100a9","Ly6g"),col=c("grey","red"),raster=F) #Neutrophils---cluster 2;19
print(fp)
dev.off()


tiff(sprintf("%s/FP_Dcs.tiff",tmpOutput),res=300,width=12,height=10,compression = "lzw",unit="in")
fp<-FeaturePlot(seu,features=c("Cd74","H2.Aa","H2.Ab1","Cd209a"),col=c("grey","red"),raster=F) #Dc---14
print(fp)
dev.off()

tiff(sprintf("%s/FP_CM.tiff",tmpOutput),res=300,width=12,height=10,compression = "lzw",unit="in")
fp<-FeaturePlot(seu,features=c("Ly6c2","F13a1","S100a4"),col=c("grey","red"),raster=F) #For cm
print(fp)
dev.off()


tiff(sprintf("%s/BM_SP_CellType.tiff",tmpOutput),res=300,width=7,height=7,compression = "lzw",unit="in")
print(DimPlot(seu,group.by="CellType",label=T,repel=T)+NoLegend())
dev.off()


############################################################################################################
seu<-readRDS(file=sprintf("%s/BM_SP.Annotated.RDS",outFolder))
Idents(seu)<-'CellType';DefaultAssay(seu)<-"RNA"
DimPlot(seu,label=T,repel=T)+NoLegend()
FeaturePlot(seu,features=c("Cd74","Ccr7","Nr4a1"),col=c('grey','red'),ncol=2)



######## function ###########
generateSubCl<-function(obj=seuST,ct=c("AdvF"),proN="AdvF",nFeature=3000,npcs=35,res=0.9,kparam=50){
  temp_CT<-subset(obj,CellType %in% ct)
  rawCount_temp<-temp_CT@assays$RNA@counts
  har_temp<-CreateSeuratObject(counts = rawCount_temp, project =proN, min.cells = 3) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nFeature) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE)
  har_temp<-AddMetaData(har_temp,metadata=temp_CT@meta.data[,c(5:ncol(temp_CT@meta.data))])
  har_temp<-har_temp %>% RunHarmony("ID", plot_convergence = FALSE)
  har_temp <- har_temp %>%
    RunUMAP(reduction = "harmony", dims = 1:npcs) %>%
    FindNeighbors(reduction = "harmony", dims = 1:npcs, k.param = kparam) %>%
    FindClusters(resolution = res) %>%
    identity()
  return(har_temp)
}
######## function ###########

#############################################
library(harmony)
library(cowplot)

har_SpBm<-generateSubCl(obj=seu, ct=c("CM"), proN="CM", nFeature=3500, npcs=50, res=0.2, kparam=20)
DimPlot(har_SpBm,label=T)
#har_SpBm<-FindClusters(har_SpBm,resolution = 0.2)
DimPlot(har_SpBm,label=T)
DimPlot(har_SpBm,group.by="Organ")
FeaturePlot(har_SpBm,features=c("Mki67"),col=c('grey','red'))


######cluster 2 is Neutrophils; Cluster 3 is DC
neID<-row.names(subset(har_SpBm,RNA_snn_res.0.2=="2")@meta.data);length(neID)
DCID<-row.names(subset(har_SpBm,RNA_snn_res.0.2=="3")@meta.data);length(DCID)

har2<-subset(har_SpBm, RNA_snn_res.0.2 %in% c("0","1"))

generateSubCl2<-function(obj=seuST,ct=c("AdvF"),proN="AdvF",nFeature=3000,npcs=35,res=0.9,kparam=50){
  temp_CT<-subset(obj,CellType %in% ct)
  rawCount_temp<-temp_CT@assays$RNA@counts
  har_temp<-CreateSeuratObject(counts = rawCount_temp, project =proN, min.cells = 3) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nFeature) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE)
  har_temp<-AddMetaData(har_temp,metadata=temp_CT@meta.data[,c(4:ncol(temp_CT@meta.data))])
  har_temp<-har_temp %>% RunHarmony("ID", plot_convergence = FALSE)
  har_temp <- har_temp %>%
    RunUMAP(reduction = "harmony", dims = 1:npcs) %>%
    FindNeighbors(reduction = "harmony", dims = 1:npcs, k.param = kparam) %>%
    FindClusters(resolution = res) %>%
    identity()
  return(har_temp)
}

harF<-generateSubCl2(obj=har2, ct=c("CM"), proN="CM", nFeature=1500, npcs=50, res=0.2, kparam=20)
DimPlot(harF,label=T)

tiff(sprintf("%s/BMSP_CM_Dim.tiff",tmpOutput),res=300,width=10,height=5,compression = "lzw",unit="in")
print(DimPlot(harF,label=T))
dev.off()

tiff(sprintf("%s/BMSP_CM_Dim_2.tiff",tmpOutput),res=300,width=10,height=5,compression = "lzw",unit="in")
print(DimPlot(harF,label=T,group.by="Organ"))
dev.off()

harF$Stage<-case_when(
  harF$ID %in% c("AB-24-8","AB-24-9","AB-24-10","AB-24-13","AB-24-14","AB-24-15")~"Resting",
  harF$ID %in% c("AB-24-11","AB-24-12","AB-24-16","AB-24-17")~"Day7")

tiff(sprintf("%s/BMSP_CM_Dim_Stage.tiff",tmpOutput),res=300,width=10,height=5,compression = "lzw",unit="in")
print(DimPlot(harF,label=F,group.by="Stage"))
dev.off()


tiff(sprintf("%s/BMSP_CM_Dim_Stage_plot2.tiff",tmpOutput),res=300,width=10,height=5,compression = "lzw",unit="in")
print(DimPlot(harF,group.by="Organ",split.by="Stage"))
dev.off()


tiff(sprintf("%s/BMSP_CM_FP_MKI67.tiff",tmpOutput),res=300,width=10,height=5,compression = "lzw",unit="in")
print(FeaturePlot(harF,features=c("Mki67","Top2a"),col=c("grey",'red')))
dev.off()
#################################################
(uniCl<-unique(harF$RNA_snn_res.0.2))
har_cl<-lapply(uniCl,function(x){
  temp<-subset(harF, RNA_snn_res.0.2 %in% x)
  tempSub<-generateSubCl2(obj=temp, ct=c("CM"), proN="CM", nFeature=1500, npcs=50, res=0.2, kparam=20)
  return(tempSub)
})

c0<-DimPlot(har_cl[[1]],label=T,group.by="Stage")
c1<-DimPlot(har_cl[[2]],label=T,group.by="Stage")

tiff(sprintf("%s/BMSP_CM_Dim_SubclusterForStage.tiff",tmpOutput),res=300,width=10,height=5,compression = "lzw",unit="in")
print(c0+c1)
dev.off()

#################################################
#################################################
seu$CellType2<-ifelse(row.names(seu@meta.data) %in% neID,"Neutrophils",seu$CellType)
seu$CellType2<-ifelse(row.names(seu@meta.data) %in% DCID,"DC",seu$CellType)
seu$Stage<-case_when(
  seu$ID %in% c("AB-24-8","AB-24-9","AB-24-10","AB-24-13","AB-24-14","AB-24-15")~"Resting",
  seu$ID %in% c("AB-24-11","AB-24-12","AB-24-16","AB-24-17")~"Day7")
head(seu)


Idents(seu)<-'CellType2';DefaultAssay(seu)<-"RNA"
DimPlot(seu,label=T,repel=T)+NoLegend()
FeaturePlot(seu,features=c("Cd74","Ccr7","Nr4a1"),col=c('grey','red'),ncol=2)
VlnPlot(seu,features=c("Cd74","Ccr7","Nr4a1"))

####CD45+, MHCII-, CD11c-, CD11b+, Ly6C
FeaturePlot(seu,features=c("Ptprc","H2.Ab1","Itgam","Itgax","Ly6c1"),col=c('grey','red'),ncol=2)



saveObj<-list(seu=seu,
              harF=harF)
saveRDS(saveObj,file=sprintf("%s/BM_SP.saveObj.RDS",outFolder))









###########################################################################################################
###########################################################################################################
saveObj<-readRDS(file=sprintf("%s/BM_SP.saveObj.RDS",outFolder))
seu<-saveObj$seu
harF<-saveObj$harF

generateSubCl2<-function(obj=seuST,ct=c("AdvF"),proN="AdvF",nFeature=3000,npcs=35,res=0.9,kparam=50){
  temp_CT<-subset(obj,CellType %in% ct)
  rawCount_temp<-temp_CT@assays$RNA@counts
  har_temp<-CreateSeuratObject(counts = rawCount_temp, project =proN, min.cells = 3) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nFeature) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE)
  har_temp<-AddMetaData(har_temp,metadata=temp_CT@meta.data[,c(4:ncol(temp_CT@meta.data))])
  har_temp<-har_temp %>% RunHarmony("ID", plot_convergence = FALSE)
  har_temp <- har_temp %>%
    RunUMAP(reduction = "harmony", dims = 1:npcs) %>%
    FindNeighbors(reduction = "harmony", dims = 1:npcs, k.param = kparam) %>%
    FindClusters(resolution = res) %>%
    identity()
  return(har_temp)
}


library(harmony)
library(cowplot)
head(harF)
DimPlot(harF,group.by="CellType",label=T)+NoLegend()
DimPlot(seu,group.by="CellType2",label=T)+NoLegend()
DimPlot(harF,group.by="Organ",label=T)+NoLegend()
DimPlot(harF,group.by="Stage",label=T)+NoLegend()

sp_day7<-subset(harF,Organ %in% "SP" & Stage=="Day7")
DimPlot(sp_day7,group.by="RNA_snn_res.0.2")
table(sp_day7$RNA_snn_res.0.2,sp_day7$ID) ##170 vs 672

sp_rest<-subset(harF,Organ %in% "SP" & Stage=="Resting")
DimPlot(sp_rest,group.by="RNA_snn_res.0.2")
table(sp_rest$RNA_snn_res.0.2,sp_rest$ID) ##99 vs 346


dp1<-DimPlot(sp_rest,group.by="RNA_snn_res.0.2",label=T);dp1
dp2<-DimPlot(sp_day7,group.by="RNA_snn_res.0.2",label=T);dp2

spData<-subset(harF,Organ %in% "SP")
dp1<-DimPlot(spData,group.by="RNA_snn_res.0.2",label=T,split.by="Stage",raster=F);dp1

tiff(sprintf("%s/spDataOnly_Dimplot.tiff",tmpOutput),res=300,width=14,height=7,compression = "lzw",unit="in")
print(dp1)
dev.off()


dfProp<-as.data.frame(prop.table(table(spData$RNA_snn_res.0.2,spData$ID),margin=2) )
colnames(dfProp)<-c("Cluster","ID","Prop")
dfProp$group<-ifelse(dfProp$ID %in% c("AB-24-16", "AB-24-17"),"Day7","Resting")



har_spD7<-generateSubCl2(obj=sp_day7, ct=c("CM"), proN="CM", nFeature=300, npcs=50, res=0.2, kparam=20)
DimPlot(har_spD7,label=T)+NoLegend()

###
DefaultAssay(har_spD7)<-"RNA";Idents(har_spD7)<-"RNA_snn_res.0.2"
MarkerD7<-FindAllMarkers(har_spD7, slot='data', test.use='wilcox', only.pos=T)
MarkerD7 %>% group_by(cluster) %>% top_n(30, avg_log2FC) %>% as.data.frame()
saveRDS(MarkerD7,file=sprintf("%s/MarkerD7.RDS",tmpOutput))

tiff(sprintf("%s/sp_day7_Dimplot.tiff",tmpOutput),res=300,width=7,height=7,compression = "lzw",unit="in")
print(DimPlot(har_spD7,label=T)+NoLegend())
dev.off()

tiff(sprintf("%s/sp_day7_mki67.tiff",tmpOutput),res=300,width=10,height=5,compression = "lzw",unit="in")
print(FeaturePlot(har_spD7,features=c("Mki67","Top2a"),col=c('grey','red')))
dev.off()


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

library(gprofiler2)
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

har_spD7 <- CellCycleScoring(har_spD7, s.features = mmus_s, g2m.features = mmus_g2m, set.ident = TRUE)
head(har_spD7)

tiff(sprintf("%s/sp_day7_cellcycl.tiff",tmpOutput),res=300,width=10,height=5,compression = "lzw",unit="in")
print(FeaturePlot(har_spD7,features=c("S.Score", "G2M.Score"),col=c('grey','red')))
dev.off()


########################
Rest<-subset(harF,Stage=="Resting")
DefaultAssay(Rest)<-"RNA";Idents(Rest)<-"Organ"
M.BM<-FindAllMarkers(Rest, slot='data', test.use='wilcox', only.pos=T)
M.BM %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% as.data.frame()
head(M.BM2[bmGene,])
c0Gene<-M.BM2[M.BM2$cluster=="SP","gene"];length(c0Gene)
head(M.BM2[c0Gene,])


############
M.BM<-M.BM[order(-M.BM$avg_log2FC),]
M.BM2<-M.BM[M.BM$p_val_adj<0.05 & M.BM$avg_log2FC>log2(1.5),];table(M.BM2$cluster)
bmGene<-M.BM2[M.BM2$cluster=="BM","gene"];length(bmGene)

har_spD7_score<-AddModuleScore(
  object=har_spD7,
  features=list(bmgene=bmGene),
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name ="BMgene",
  seed = 1,
  search = FALSE,
  slot = "data")
head(har_spD7_score)


tiff(sprintf("%s/sp_day7_FP_bmGene.tiff",tmpOutput),res=300,width=7,height=7,compression = "lzw",unit="in")
print(FeaturePlot(har_spD7_score,features=c("BMgene1"),col=c('grey','red')))
dev.off()

sp_day7_sub<-subset(har_spD7_score,RNA_snn_res.0.2==0);DimPlot(sp_day7_sub,label=T)
spD7Met<-sp_day7_sub@meta.data
spD7Met<-spD7Met[order(spD7Met$BMgene1),]
spD7Met$indes<-1:nrow(spD7Met)
spD7Met$newCat<-ifelse(spD7Met$BMgene1>0.25,"High","Mid")
spD7Met$newCat<-ifelse(spD7Met$BMgene1<0.05,"low",spD7Met$newCat)
table(spD7Met$newCat)


tiff(sprintf("%s/sp_day7_ScatterPlot.tiff",tmpOutput),res=300,width=10,height=7,compression = "lzw",unit="in")
plot(BMgene1~indes,data=spD7Met,xlab="index",ylab="Score for BM gene")
abline(h=0.05,col='blue',lty=2)
abline(h=0.25,col='red',lty=2)
dev.off()


spD7Met_back<-spD7Met[row.names(sp_day7_sub@meta.data),]
sum(row.names(spD7Met_back)==row.names(sp_day7_sub@meta.data))

sp_day7_sub$NewCat<-spD7Met_back$newCat
Idents(sp_day7_sub)<-"NewCat"

dp<-DimPlot(sp_day7_sub,group.by="NewCat")
vn<-VlnPlot(sp_day7_sub,"BMgene1",split.by='NewCat')

tiff(sprintf("%s/sp_day7_Dp_NewCat.tiff",tmpOutput),res=300,width=10,height=7,compression = "lzw",unit="in")
print(dp)
dev.off()

tiff(sprintf("%s/sp_day7_vn_NewCat.tiff",tmpOutput),res=300,width=10,height=7,compression = "lzw",unit="in")
print(vn)
dev.off()


DefaultAssay(sp_day7_sub)<-"RNA";Idents(sp_day7_sub)<-"NewCat"
M.sp<-FindAllMarkers(sp_day7_sub, slot='data', test.use='wilcox', only.pos=T)
M.sp %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% as.data.frame()

tiff(sprintf("%s/sp_day7_vn_3sigGene.tiff",tmpOutput),res=300,width=10,height=7,compression = "lzw",unit="in")
print(VlnPlot(sp_day7_sub,features=c("Lgals1","Ly6c2","Tubb5"),split.by="NewCat"))
dev.off()

head(sp_day7_sub)
M.hvsl<-FindMarkers(sp_day7_sub, 
                    ident.1="High",
                    ident.2='low',
                    slot='data', test.use='wilcox', only.pos=F)
M.hvsl<-M.hvsl[order(M.hvsl$p_val_adj),]
M.hvsl_sig<-M.hvsl[M.hvsl$p_val_adj<0.05,]
write.csv(M.hvsl_sig,file=sprintf("%s/M.hvsl_sig.csv",tmpOutput))



#########################################################
geneToPlot<-c("Il1b","Nlrp3","Tgfbr1","Tgfbr2","Itgam")
fp<-FeaturePlot(seu,features=geneToPlot,col=c('grey','red'),ncol=3)

tiff(sprintf("%s/additionGene_May9.tiff",tmpOutput),res=300,width=18,height=12,compression = "lzw",unit="in")
print(fp)
dev.off()

