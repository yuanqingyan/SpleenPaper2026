
library(edgeR)
FolderDay1to5<-"./BulkRNAseq/TO/Oct2022/RawCount"

outputFolder<-"./RNAseq/Com2022/"
SaveFolder<-"./BulkRNAseq/TO/Oct2022/RDS"
source("~/code/sourceRcode/bulkRNAseq_sourceCode.R")
listCPM<-readRDS(file=sprintf("%s/listCPM_5Days.rds",SaveFolder))
cpm<-listCPM$cpm
write.csv(cpm,file=sprintf("%s/CPM_DiffDay_timecourse.csv",outputFolder))
cpm_noLog<-cpm(listCPM[[2]],log=FALSE,normalized.lib.size=TRUE);sum(cpm_noLog[,1]<0)

#################################################################################################################################
################################ post analysis--plot some genes                ##################################################
#################################################################################################################################
targetGene<-"Csf1r"
tarGenedf<-data.frame(ID=colnames(cpm),data=as.numeric(cpm[grepl("Csf1r",row.names(cpm)),]));head(tarGenedf)
tarGenedf$cat<-sapply(strsplit(tarGenedf$ID,split="\\."),function(x) sprintf("%s.%s",x[1],x[2]));tarGenedf

typeN="Csf1r"
tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=6,height=6,units="in",compression="lzw")
par(mar=c(8,4,1,1))
boxplot(data~cat,tarGenedf,las=2,ylab="CPM",xlab="")
dev.off()


targetGene<-"P2rx7"
tarGenedf<-data.frame(ID=colnames(cpm),data=as.numeric(cpm[grepl("Csf1r",row.names(cpm)),]));head(tarGenedf)
tarGenedf$cat<-sapply(strsplit(tarGenedf$ID,split="\\."),function(x) sprintf("%s.%s",x[1],x[2]));tarGenedf

typeN="P2rx7"
tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=6,height=6,units="in",compression="lzw")
par(mar=c(8,4,1,1))
boxplot(data~cat,tarGenedf,las=2,ylab="CPM",xlab="")
dev.off()


#################################################################################################################################
################################ Filter data by cv                ###############################################################
#################################################################################################################################
library(dplyr);library(tidyverse);library(gplots);library("RColorBrewer")
cpm_noLog<-cpm_noLog[!apply(cpm_noLog,1,function(x){sum(x==0)==ncol(cpm_noLog)}),]
library(matrixStats)
cv<-rowSds(cpm_noLog)/rowMeans(cpm_noLog)

typeN="Histogram_CV"
sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=6,height=6,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=6,height=6)}
  hist(cv,breaks=30);abline(v=0.45,col="red",lty=2);sum(cv>=0.8)
  dev.off()
})

cpm_sel<-cpm_noLog[cv>=0.45,];dim(cpm_sel) ##2412 gene
scaledata <- t(scale(t(cpm_sel),center=TRUE,scale=TRUE)) #
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # 
TreeC = as.dendrogram(hc, method="complete")
plot(TreeC, main = "Sample Clustering",  ylab = "Height")

hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") #
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") #
heatmap.2(scaledata,Rowv=as.dendrogram(hr),Colv=NULL,col="bluered",scale="row",margins = c(7, 7),cexCol = 0.7,labRow = F,main = "",trace = "none")
heatmap.2(scaledata,Rowv=as.dendrogram(hr),Colv=NULL,col=rev(brewer.pal(11,"RdBu")),scale="row",margins = c(7, 7),cexCol = 0.7,labRow = F,main = "Heatmap.2",trace = "none")

#################################################################################################################################
################################ K-Means: How many clusters?      ###############################################################
#################################################################################################################################
##
wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata, centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

typeN="SSE"
sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=6,height=6,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=6,height=6)}
  plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
  dev.off()
})


library(cluster)
sil <- rep(0, 20)
for(i in 2:20){
  k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(scaledata))
  sil[i] <- mean(ss[, 3])}#

typeN="silhouette"
sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=6,height=6,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=6,height=6)}
  plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
  abline(v = which.max(sil), lty = 2)#
  dev.off()
})
cat("Average silhouette width optimal number of clusters:", which.max(sil), "\n")


library(vegan)
fit <- cascadeKM(scaledata, 1, 20, iter = 100)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
typeN="Calinsky"
sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=8,height=6,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=8,height=6)}
  plot(fit, sortg = TRUE, grpmts.plot = TRUE)
  dev.off()
})
calinski.best <- as.numeric(which.max(fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

#Gap statistic
library(cluster)
set.seed(1234);gap <- clusGap(scaledata, kmeans, 20, B = 100, verbose = interactive())
plot(gap, main = "Gap statistic");abline(v=which.max(gap$Tab[,3]), lty = 2)

# Find the optimal k using the 1-standard-error rule
best_k <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")
print(best_k)

plot(gap, main = "Gap statistic")
abline(v = best_k, col = "red", lty = 2)

#################################################################################################################################
################################         Clustering the data      ###############################################################
#################################################################################################################################
library(ggplot2);library(reshape)
set.seed(12345)
kClust <- kmeans(scaledata, centers=2, nstart = 300, iter.max = 20)
kClusters <- kClust$cluster;table(kClusters)
dfCluster<-as.data.frame(kClusters)
write.csv(dfCluster,file=sprintf("%s/ClusterGenes_March2026.csv",outputFolder))###
dfCluster[grepl("Il1b",row.names(dfCluster)),]
dfCluster[grepl("Nrlp3",row.names(dfCluster)),]
dim(cpm_sel)

#cluster ‘cores’ aka centroids
clust.centroid = function(i, dat, clusters) {ind = (clusters == i);colMeans(dat[ind,])}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata, kClusters)
Kmolten <- melt(kClustcentroids);colnames(Kmolten) <- c('sample','cluster','value')
Kmolten$Day<-sapply(strsplit(as.character(Kmolten$sample),split="\\."),function(x) x[2]);table(Kmolten$Day)


df_km <- Kmolten %>%
  group_by(Day, cluster) %>%
  summarise(mean_val = mean(value), .groups = 'drop'); head(df_km)


set.seed(1235)
(p0<-ggplot(Kmolten, aes(x = factor(Day), y = value, color = factor(cluster), fill = factor(cluster))) +
  geom_boxplot(alpha = 0.2, position = position_dodge(width = 0.0), outlier.shape = NA) +
  geom_jitter(position = position_dodge(width = 0.5), alpha = 1) +
  stat_summary(fun = mean, geom = "point", 
    size = 3,position = position_dodge(width = 0)) +
  stat_summary(
    fun = mean, 
    geom = "line", 
    aes(group = factor(cluster)), 
    linewidth = 1,
    position = position_dodge(width = 0)
  ) +
  xlab("Days") + 
  ylab("Expression") + 
  labs(title= "Gene Clustering",color = "Cluster")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw())

typeN="ClusterCentroids_kmeans"
sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=6,height=4,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=6,height=4)}
  print(p0)
  dev.off()
})



(p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + geom_point() +  geom_line() + xlab("Different Samples/Days") + ylab("Expression") + labs(title= "Gene Clustering",color = "Cluster")+theme(axis.text.x = element_text(angle = 45, hjust = 1)))
typeN="ClusterCentroids_kmeans_2"
sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=8,height=6,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=8,height=6)}
  print(p1)
  dev.off()
})
cor(kClustcentroids)


#################################################################################################################################
################################   Using a cluster score to identify core genes    #############################################
#################################################################################################################################
sapply(unique(kClusters),function(ic){
  core <- Kmolten[Kmolten$cluster==ic,]
  K <- (scaledata[kClusters==ic,])
  corscore <- function(x){cor(x,core$value)}#
  score <- apply(K, 1, corscore)
  Kmolten <- melt(K)#
  colnames(Kmolten) <- c('gene','sample','value')
  Kmolten <- merge(Kmolten,score, by.x='gene',by.y='row.names', all.x=T)#
  colnames(Kmolten) <- c('gene','sample','value','score')
  Kmolten$order_factor <- 1:length(Kmolten$gene)#
  Kmolten <- Kmolten[order(Kmolten$score),]#
  Kmolten$order_factor <- factor(Kmolten$order_factor , levels = Kmolten$order_factor)#s
  (p2 <- ggplot(Kmolten, aes(x=sample,y=value)) +  geom_line(aes(colour=score, group=gene)) + scale_colour_gradientn(colours=c('blue1','red2')) +  geom_line(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +xlab("Samples/Days") + ylab("Expression") + labs(title= sprintf("Cluster %s Expression by Samples/Days",ic),color = "Score")+theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  typeN=sprintf("Cluster_%s_expression_3",ic)
  sapply(1:2,function(ip){
    if(ip==1){tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=8,height=6,units="in",compression="lzw")}
    if(ip==2){pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=8,height=6)}
    print(p2)
    dev.off()
  })
})

#################################################################################################################################
################################      heatmap                                      #############################################
#################################################################################################################################
library(pheatmap)
pSca<-as.data.frame(cbind(scaledata,kClusters));head(pSca)
pSca<-pSca[order(pSca$kClusters),];head(pSca)

pheat<-pheatmap(pSca[,1:(which(colnames(pSca)=="kClusters")-1)], 
         cluster_rows = F, 
         cluster_cols = F, 
         scale="row",
         col=colorRampPalette(c("blue", "white", "red"))(200), 
         legend=TRUE, 
         show_rownames=FALSE, 
         show_colnames=TRUE)

typeN="pheat_timecourse_V1"
sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=8,height=6,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=8,height=6)}
  print(pheat)
  dev.off()
})

pheat2<-pheatmap(pSca[,1:(which(colnames(pSca)=="kClusters")-1)], 
         cluster_rows = F, 
         cluster_cols = F, 
         scale="row",
         gaps_row = cumsum(as.numeric(table(pSca$kClusters)))[1:2],
         col=colorRampPalette(c("blue", "white", "red"))(200), 
         legend=TRUE, 
         show_rownames=FALSE, 
         show_colnames=TRUE)

typeN="pheat_timecourse_ver2"
sapply(1:2,function(ip){
  if(ip==1){tiff(sprintf("%s/%s.tiff",outputFolder,typeN), res=600, width=8,height=6,units="in",compression="lzw")}
  if(ip==2){pdf(sprintf("%s/%s.pdf",outputFolder,typeN),width=8,height=6)}
  print(pheat2)
  dev.off()
})

pheatmap(pSca[,1:(which(colnames(pSca)=="kClusters")-1)], 
         cluster_rows = F, 
         cluster_cols = F, 
         scale="row",
         kmeans_k=3,
         col=colorRampPalette(c("blue", "white", "red"))(200), 
         #breaks=pSca$kClusters,
         legend=TRUE, 
         show_rownames=FALSE, 
         show_colnames=TRUE)

#################################################################################################################################
################################      pcluster gene and pathway                     #############################################
#################################################################################################################################
uniqMod<-unique(pSca$kClusters)
moduleGene<-lapply(uniqMod,function(x){
  ID.gene<-row.names(pSca[pSca$kClusters==x,])
  gene<-sapply(strsplit(ID.gene,split="\\_"),function(y) y[2])
  out<-list(IDGene=ID.gene,gene=gene);return(out)});names(moduleGene)<-uniqMod
InGene<-lapply(moduleGene,function(x) x$gene);InGene
names(InGene)<-sprintf("Cluster_%s",uniqMod)

library(clusterProfiler)
library(msigdbr);msigdbr_species();
M.df<- as.data.frame(msigdbr(species = "Mus musculus"));unique(M.df$gs_cat)
set_H <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, gene_symbol)

enrich_H<-lapply(InGene,function(x){
  H_temp<-enricher(gene=x,
                   TERM2GENE = set_H,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   minGSSize = 15,
                   maxGSSize = 500,
                   qvalueCutoff = 1,
                   TERM2NAME = NA)
  return(H_temp)})
lapply(enrich_H,function(x) sum(x$p.adjust<0.05))
enrich_H.sig<-lapply(enrich_H,function(x) x[x$p.adjust<0.05,])

lapply(1:length(enrich_H),function(x) {
  nam<-names(enrich_H)[x]
  temp<-enrich_H[[x]]
  write.table(temp, file=sprintf("%s/%s_H.xls",outputFolder,nam), sep="\t", row.names=F,col.names=T)
})

lapply(1:length(enrich_H.sig),function(x) {
  nam<-names(enrich_H.sig)[x]
  temp<-enrich_H.sig[[x]]
  write.table(temp, file=sprintf("%s/%s_H.sig.xls",outputFolder,nam), sep="\t", row.names=F,col.names=T)
})

lapply(enrich_H.sig,function(x) {nrow(x)})

library(ggplot2)
lapply(1:length(enrich_H.sig),function(x) {
  nam<-names(enrich_H.sig)[x]
  temp<-enrich_H.sig[[x]]
  temp$Negtive_log_adj.p<-(-log10(temp$p.adjust))
  temp$shortP<-gsub("HALLMARK_","",temp$ID)
  temp$shortP <- factor(temp$shortP, levels = temp$shortP[order(temp$Negtive_log_adj.p)])
  temp$col<-"grey"
  temp$height<-ifelse(nrow(temp)>50,24,8)
  
  sapply(1:2,function(iloop){
    if(iloop==1){tiff(sprintf("%s/H_%s_forpaper.tiff",outputFolder,nam),res=300,width=10,height=temp$height[1],compression = "lzw",unit="in")}
    if(iloop==2){pdf(sprintf("%s/H_%s_forpaper.pdf",outputFolder,nam),width=10,height=temp$height[1])}
    
    
    p1 <- ggplot(data=temp, aes(x=shortP, y=Negtive_log_adj.p)) +
      geom_bar(stat="identity", aes(fill=col)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
      geom_text(aes(label=shortP, y= 0.05), vjust=0.5, color="black", size=3.5, hjust=0) +
      coord_flip() +
      theme_classic() + 
      theme(
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.y = element_blank()      
      )
    
    print(p1)
    dev.off()
  })
})






