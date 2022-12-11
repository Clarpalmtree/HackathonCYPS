#!/usr/bin/Rscript --slave

# rm(list=objects())
# graphics.off()
set.seed(123)

library("DESeq2")

args = commandArgs(T)
counts_file = args[1]
dir = args[2]

meta<-read.table(paste(dir, "/metaData.txt", sep=""), header = TRUE, sep = ",",row.names = 1)
Counts = read.table(paste(dir, "/count/output.counts", sep=""), header=T)    # Loading count matrix file
# ## DATA ------------------------------------------------------------------------
# meta<-read.table("metaData.txt", header = TRUE, sep = ",",row.names = 1)
# Counts<-read.table("output.counts", header = TRUE, sep = "",row.names = 1)

# Data formatting and Data cleaning --------------------------------------------

## Dataframe
# Removing the first 5 columns (unecessary for analysis) 
new_Counts= Counts[,-c(1:5)]
# Removing ".bam" extension in column names
colnames(new_Counts)[c(1:dim(new_Counts)[2])] = gsub('.{0,4}$', '', colnames(new_Counts)[c(1:dim(new_Counts)[2])])
# Sorting in order : 
new_Counts<-new_Counts[,rownames(meta)]

print("Is data correctly sorted as in metadata ?")
print(all(rownames(meta) == colnames(new_Counts)))

## Labels
# Getting samples ID
labels=c()
for (i in colnames(new_Counts)){
  labels <- c(labels,i)
}
# Getting mutation status
mutations=c()
for (label in labels){
  mutations <-c(mutations,meta[label,"Type"])
}
# Sample ID and tissue type correspondance
labels2 = data.frame(labels, mutations)
labels2$mutations<-factor(labels2$mutations,levels=c("WT","M"))


new_Counts = new_Counts[rowSums(new_Counts)>0,]   
new_Counts = as.matrix(new_Counts) 
new_Counts <- data.frame(new_Counts)

## Removal of all unexpressed genes
dds <- DESeqDataSetFromMatrix(countData=new_Counts, colData=labels2, design= ~ mutations)
dds$mutations <- relevel(dds$mutations, ref ="WT") #set WT as reference for log2FC
dds <- DESeq(dds) # Levels: WT M
res <- results(dds)
## Removing NA values
res <- na.omit(res)

## PCA -------------------------------------------------------------------------

tran_data<-t(new_Counts)
vsdata <- vst(dds, blind=FALSE)
#export
png("PCA_DE.png",width=600)
plotPCA(vsdata, intgroup="mutations")
dev.off()

## Heatmap ---------------------------------------------------------------------

rld<-rlog(dds,blind=TRUE)
# Extract the rlog matrix from the object
rld_mat <- assay(rld)
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)
png("heatmap_de.png",width=700, height=600)
heatmap(rld_cor)
dev.off()

## Differential analysis -------------------------------------------------------
## Calculating genes parameters
res_mat = as.matrix(res)
write.csv(res_mat, "DESeq_results.csv")
## DE genes
DEgenes = which(res$padj<0.05) # DE genes with pvalue < 5%
significative_DEgenes = res[DEgenes,]   # matrix containing names, log2FC and adjusted p-value of significantly differentially expressed genes
significative_DEgenes = significative_DEgenes[, c(2,6)]
significative_DEgenes <- significative_DEgenes[order(significative_DEgenes$padj),]
write.csv(significative_DEgenes, "significative_DEgenes.csv") 
top_10=significative_DEgenes[1:10, ]
write.csv(top_10,"top_10_de_genes.csv")

## Summary table -----------------------------------------------------------------

summary.df = as.data.frame(matrix("", ncol= 2, nrow = 2))  # DataFrame whose summarise previous matrix
rownames(summary.df) = c("overexpressed genes in Wild Types", "overexpressed genes in Mutants")
colnames(summary.df) = c("Amount", "Ratio (significantly differentially expressed)")
tot = dim(significative_DEgenes)[1]
summary.df[1,1] = length(which(significative_DEgenes$log2FoldChange < -1)) #overexpressed in WT
summary.df[2,1] = length(which(significative_DEgenes$log2FoldChange > 1)) #overexpressed in M

summary.df[1,2] = round(as.numeric(summary.df[1,1]) / tot, digits = 2)
summary.df[2,2] = round(as.numeric(summary.df[2,1]) / tot, digits = 2)
summary.df
write.csv(summary.df, "summary.csv")
