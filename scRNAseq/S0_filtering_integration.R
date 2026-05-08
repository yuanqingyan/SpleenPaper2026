library(dplyr)
library(Seurat)
library(Matrix)
library(singleGEO)
library(SoupX)
library(reticulate)
library(anndata)
library(sceasy)
library(scDblFinder)


minFeature <- 200  #
maxFeature <- 7500
minCount <- 400    #
maxCount <- 40000
maxMT <- 10        #


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
RawDataList_soupx<-sapply(pData$V1,function(x) {
  sample_i<-x;print(sample_i)
  
  SCReadFolder<-sprintf("%s/%s/outs",CountFolder,sample_i)
  sc = load10X(SCReadFolder)
  sc = autoEstCont(sc)
  out_forSeu = adjustCounts(sc)

  raw_dat <- CreateSeuratObject(counts = out_forSeu,
                                project =sprintf("%s",sample_i),
                                min.cells = 1,
                                min.features = 1)
  
  set.seed(12345)
  raw_sce <- scDblFinder(GetAssayData(raw_dat, slot="counts"))
  raw_dat$scDblFinder.class=raw_sce$scDblFinder.class
  raw_dat$scDblFinder.score=raw_sce$scDblFinder.score
  
  raw_dat[['percent.mt']] <- PercentageFeatureSet(raw_dat, pattern = '^mt-') #
  raw_dat[['percent.rps']] <- PercentageFeatureSet(raw_dat, pattern = '^rps') 
  raw_dat[['percent.rpl']] <- PercentageFeatureSet(raw_dat, pattern = '^rpl') 
  raw_dat$percent.rp <- raw_dat$percent.rps + raw_dat$percent.rpl 

  raw_dat[['patient']] <- sample_i
  raw_dat[['ID']] <- sample_i

  tiff(sprintf("%s/%s_FeaturePlot_rawData_soupx.tiff",tmpOutput,sample_i),res=300,width=10,height=7,compression = "lzw",unit="in")
  VlnPlot<-VlnPlot(raw_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(VlnPlot)
  dev.off()

  return(raw_dat)
})
if(is.null(names(RawDataList_soupx))){names(RawDataList_soupx)<-pData$V1}
saveRDS(RawDataList_soupx, file =sprintf("%s/RawDataList_%s.Soupx_beforeFiltering.rds",tmpOutput,objName))

#########################################################################################################
############################   Filter data                           ####################################
#########################################################################################################
list.files(outFolder)
rawIs<-readRDS(sprintf("%s/RawDataList_%s.Soupx_beforeFiltering.rds",tmpOutput,objName))
lapply(rawIs,function(x) dim(x))
filter.dat<-lapply(rawIs,function(x){
  temp<-subset(x,
               subset = nFeature_RNA > minFeature &
                 nFeature_RNA < maxFeature &
                 nCount_RNA > minCount &
                 nCount_RNA < maxCount &
                 percent.mt < maxMT &
                 scDblFinder.class %in% "singlet")
  filterDat<-GetAssayData(temp,slot="counts")
  return(filterDat)
});lapply(filter.dat,function(x) dim(x))

filter.Seu<-lapply(rawIs,function(x){
  temp<-subset(x,
               subset = nFeature_RNA > minFeature &
                 nFeature_RNA < maxFeature &
                 nCount_RNA > minCount &
                 nCount_RNA < maxCount &
                 percent.mt < maxMT &
                 scDblFinder.class %in% "singlet")
  return(temp)
});lapply(filter.Seu,function(x) dim(x))
saveRDS(filter.Seu, file =sprintf("%s/filter_%s.Seu.rds",tmpOutput,objName))

#####################################################################################################################
rawSeu<-MakeSeuObj_FromRawRNAData(RawList=filter.dat,
                                  GSE.ID=sprintf("%s",objName),
                                  MtPattern='^mt-',
                                  MinFeature=1,
                                  MaxFeature=750000,
                                  MinCount=1,
                                  MaxCount=4000000,
                                  MaxMT=10,
                                  Norm.method="lognorm",
                                  Feature.selection.method = "vst",
                                  Resolution=1.2,
                                  Nfeatures = 3000)
saveRDS(rawSeu, file=sprintf("%s/%s_Soupx_rawSeu.rds",outFolder,objName))

seu <- SeuObj_integration(Object.list =rawSeu,
                          Object.list2 = NULL,
                          Frow.which.Norm="lognorm",
                          SampleNameAsReference=NULL,
                          NumberOfSampleForReference=numRef,
                          Nfeatures = 3000,
                          Do.scale = TRUE,
                          Do.center = TRUE,
                          Anchor.reduction = 'rpca',
                          Dims.anchor=1:30,
                          Dims.umap=1:30,
                          Resolution=1.5,
                          Future.globals.maxSize = 96000*1024^2)

#######################################################################
saveRDS(seu,file=sprintf("%s/%s_Soupx_beforeAnnotation.rds",outFolder,objName))


