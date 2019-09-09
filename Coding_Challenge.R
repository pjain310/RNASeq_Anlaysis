#Load required libraries
library(limma)
library(edgeR)
library(Glimma)
library(RColorBrewer)
library(fpc)
library(cluster)

#Set wd as datset folder
setwd("Downloads/Tempus_RNA_coding_challenge_internship/datasets")

#Load file into dataset, get counts matrix
data=read.csv("data.csv")
countdata=data[,-1]

#Filter genes with lower cpm than threshold (cpm > 1)
data.cpm <- cpm(countdata)
threshold <- data.cpm > 1

#Filter long non-coding RNA seqs (represented as: c<chromosome number>orf<orf number>)
noncoding <- grepl("orf",data[,1])
keep <- (rowSums(threshold)) > 30

#Filtering rows based on criteria
counts <- countdata[keep,]
genes <- data[keep,1]

#Find correlation between samples using Pearson Correlation and cluster to identify outliers
sample_cor <- cor(counts,use="p")
cluster1<-hclust(as.dist(1-sample_cor),method="complete")
plot(cluster1,cex=0.7,labels=names(countdata),main="Cluster Dendrogram for Outlier Identification",xlab="Samples",ylab="Distance")

#Remove outliers based on clustering (find a way to automate)
counts$cancer8=NULL
counts$cancer58=NULL
counts$normal87=NULL
counts$cancer98=NULL
counts$normal27=NULL
counts$normal34=NULL
counts$normal=NULL
counts$cancer100=NULL

#Store counts matrix as DGEList
y<-DGEList(counts=counts,lib.size = colSums(counts), genes=genes, group=ifelse(grepl("normal",names(counts)),"normal","cancer"))
design<-model.matrix(~0+y$samples$group)

#Estimate Dispersion of gene counts using negative binomial distribution
y <- estimateDisp(y,design)
sqrt(y$common.dispersion) # biological coefficient of variation
plotBCV(y) #From BCV, we see that variation does not follow a strict trend, therefore we do not develop an LMM

#Use Fischer's Exact test to estimate dispersion between two sample groups
et <- exactTest(y,pair = c("normal","cancer"))

#Use FDR correction to find significantly expressed genes 
decide_r<-decideTestsDGE(et, adjust.method = "BH")

#Find top upregulated and downregulated genes
results_edgeR <- topTags(et,n=nrow(et))
sig_genes <- results_edgeR[as.logical(decide_r),] 
upreg_genes <- sig_genes$table$logFC > 0
dnreg_genes <- sig_genes$table$logFC < 0

#Print top 10 upregulated and downregulated genes
sig_genes[upreg_genes,][1:10,]
sig_genes[dnreg_genes,][1:10,]

#Show trends in data before/after normalisation
col<-brewer.pal(12,"Paired") #Colour box plots (12 is max)

lcpm <-cpm(y,log=TRUE)
boxplot(lcpm, las=2, col=col, main="Box Plot: Before Normalisation")

#Normalise counts using RLE
y <- calcNormFactors(y,method="RLE")

lcpm <-cpm(y,log=TRUE)
boxplot(lcpm, las=2, col=col, main="Box Plot: After Normalisation")

#Plot MDS plot for only tumor samples to observe subtypes
glMDSPlot(lcpm[,49:144], labels = colnames(lcpm[,49:144]), main="MDS Plot: Tumor Subtyping")

#Find correlation between samples using Pearson Correlation and cluster to identify outliers
tumor_cor <- cor(lcpm[,49:144],use="p")

#Use hclust to give cluster dendrogram
cluster_t<-hclust(as.dist(1-tumor_cor),method="complete")
plot(cluster_t,cex=0.7,labels=names(counts)[49:144],main="Cluster Dendrogram for Tumor Subtype Identification",xlab="Samples",ylab="Distance")

#Use PAM (partitioning around medoids) to croscheck maximum number of clusters
pamk.tumor<-pamk((1-tumor_cor),krange = 1:6, criterion = "asw")
cat("number of clusters estimated by optimum average silhouette width:", pamk.tumor$nc, "\n")
clusplot(pam((1-tumor_cor),pamk.tumor$nc),main = "Classification of Tumor into Subtypes by Partitioning Around Medoids")
 
#Bootstrap samples to increase confidence 
pamk.boot<-clusterboot((1-tumor_cor),10,distances=TRUE,bootmethod = "subset",clustermethod=pamkCBI)

#Calculate standard deviation for twos group
cat("The standard deviation for subtype 1: ",sd(tumor_cor[t_1,t_1]))
cat("The standard deviation for subtype 2: ",sd(tumor_cor[t_2,t_2]))
cat("The standard deviation for tumor group: ",sd(tumor_cor))
par(mfrow=c(2,2))

#Create subsets of samples in group 1 (larger) and draw cluster plots (12 iterations)
t_1<-pamk.tumor$pamobject$clustering==1
t_2<-pamk.tumor$pamobject$clustering==2
y<-c(8:70)
tnames<-names(pamk.tumor$pamobject$clustering)
for (j in c(1:12)){
subtypes<-c()
for (i in c(8:70)){
  t1<-logical(length=sum(t_1))
  t2<-logical(length=sum(t_2))
  names(t1)<-tnames[t_1]
  names(t2)<-tnames[t_2]
  
  length=sum(t_1)
  
  t<-logical(length=length(tnames))
  names(t)<-tnames
  
  #Sample randomly from subset 1
  names1<-sample(names(t1),i)
  t[names1]=TRUE
  
  #Sample all from subset 2
  t[names(t2)]=TRUE
  
  mat<-1-tumor_cor[t,t]
  
  pk.tum<-pamk(mat,krange = 2:7, criterion = "asw")
  subtypes<-c(subtypes,pk.tum$nc)
}

plot(y,subtypes,type="l",xlab = "Number of samples taken from group 1")
abline(v=20,col="blue")
  }