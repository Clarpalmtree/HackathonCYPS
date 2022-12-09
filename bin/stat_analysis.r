## HACKATON ANALYSIS

#### Data Loading
# Counts Matrix
args = commandArgs(T)
counts_file = args[1]
dir = args[2]
path1 = paste(dir, "/count/", sep="")
path2 = paste(dir, sep="")

meta<-read.table(paste(path2,"metaData.txt"), header = TRUE, sep = ",",row.names = 1)
Counts = read.table(paste(path2, counts_file, sep=""), header=T)    # Loading count matrix file


##### Set Seed #####
set.seed(123)

#### PACKAGES
# install.packages("pheatmap")
# install.packages("factoextra")
# if (!requireNamespace('BiocManager', quietly = TRUE)){
#   install.packages('BiocManager')
# }
# BiocManager::install('EnhancedVolcano')
# BiocManager::install('apeglm')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
install.packages("plotly")
library("DESeq2")
library("pheatmap")
library("FactoMineR")
library("factoextra")
library("tidyverse")
library("FactoMineR")
library("apeglm")
library ("RColorBrewer") 
library("ggrepel")
library("EnhancedVolcano")
library("magrittr") # needs to be run every time you start R and want to use %>%
library("dplyr")    # alternatively, this also loads %>%
library("ggplot2")
library("plotly")


# Data formatting and Data cleaning --------------------------------------------

## Dataframe
# Removing the first 5 columns (unecessary for analysis) 
new_Counts= Counts[,-c(1:5)]
# Removing ".bam" extension in column names
colnames(new_Counts)[c(1:8)] = gsub('.{0,4}$', '', colnames(new_Counts)[c(1:8)])
colnames(new_Counts)
# Sorting in order : 
new_Counts<-new_Counts[,rownames(meta)]
colnames(new_Counts)
all(rownames(meta) == colnames(new_Counts))

## Labels
labels = c("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
mutations = c("WT", "M", "WT", "WT", "WT", "WT", "M", "M")
labels2 = data.frame(labels, mutations)

## Ordering the tissue types so that it is sensible and make sure the control sample is first: WT -> M:mutÃ©
labels2$mutations<-factor(labels2$mutations,levels=c("WT","M"))
## Removal of all unexpressed genes
dds <- DESeqDataSetFromMatrix(countData=new_Counts, colData=labels2, design= ~ mutations)
dds <- DESeq(dds) # Levels: M WT
dds$sizeFactor
res <- results(dds)
summary(res)
new_Counts = new_Counts[rowSums(new_Counts)>0,]   
new_Counts=as.matrix(new_Counts) 
new_Counts<- data.frame(new_Counts)
#head(new_Counts)

# I- RNA-seq count distribution ------------------------------------------------

hist=ggplot(new_Counts) +geom_histogram(aes(x =SRR628582), stat = "bin", bins = 200) + xlim(-5, 500)  +xlab("Raw expression counts") + ylab("Number of genes")
png("rawcount_vs_nbgene.png",width=600)
plot(hist)
dev.off()

## Modeling count data for Mutated samples -------------------------------------

mean_counts <- apply(new_Counts[c(2,7,8)], 1, mean)
variance_counts <- apply(new_Counts[c(2,7,8)], 1, var)
df <- data.frame(mean_counts, variance_counts)
Disp=ggplot(df) + geom_point(aes(x=mean_counts, y=variance_counts)) + geom_line(aes(x=mean_counts, y=mean_counts, color="red")) + scale_y_log10() + scale_x_log10()
# exporting figure
png("MeanVariance_Mutatedsamples.png")
plot(Disp)
dev.off()

## Modeling count data for WT samples ------------------------------------------

mean_counts <- apply(new_Counts[c(1,3,4,5,6)], 1, mean)
variance_counts <- apply(new_Counts[c(1,3,4,5,6)], 1, var)
df <- data.frame(mean_counts, variance_counts)
Disp=ggplot(df) + geom_point(aes(x=mean_counts, y=variance_counts)) + geom_line(aes(x=mean_counts, y=mean_counts, color="red")) + scale_y_log10() + scale_x_log10()
# exporting figure
png("MeanVariance_WTsamples.png")
plot(Disp)
dev.off()

# II - ANALYSIS ----------------------------------------------------------------

## ACP -------------------------------------------------------------------------

tran_data<-t(new_Counts)
res_pca<-PCA(tran_data,graph = FALSE)
#res_pca$eig
plot_pca=fviz_pca_ind (res_pca,
                       repel = TRUE,
                       col.ind = mutations, # coloring by tissue type
                       palette = c("#009999", "#0000FF"),
                       legend.title = "Sample Type",
                       addEllipses = TRUE,
                       ellipse.level=0.95,
                       mean.point = FALSE)
png("PCA.png",width=600)
plot(plot_pca)
dev.off()

## HEATMAP ---------------------------------------------------------------------

## Data
#Transform normalized counts using the rlog transformation
rld<-rlog(dds,blind=TRUE)
# Extract the rlog matrix from the object
rld_mat <- assay(rld)
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)

## Plot heatmap
heat.colors <- brewer.pal(6, "Purples")
# exporting
# png("heatmap.png",width=600)
# heatmap=pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
#                  fontsize_row = 10, height=20)
# dev.off()

##
sizeFactors(dds)
res <- results(dds)
res<- results(dds,contrast = c("mutations","M","WT"))
class(res)
summary(res)
mcols(res, use.names=T)
#res %>% data.frame() %>% View()

# Remoing 'NA' values
res_NO_NA<-res[complete.cases(res[,6]),]
#res_NO_NA %>% data.frame() %>% View()


## Parameters
p_adj<-0.05
LFC<-0.5
res_sig<-res_NO_NA%>%data.frame() %>%rownames_to_column(var="gene") %>% as_tibble()
sig_DE_gene <- res_sig %>%filter(padj < p_adj & abs(log2FoldChange) >LFC)
write.csv(sig_DE_gene, file="gÃ¨nes_df_sig.csv", row.names = FALSE )
summary(sig_DE_gene)

## Pvalue histogram
#exporting
png("histogram_pvalue.png")
hist(res_sig$pvalue)
dev.off()

sig_DE_gene 

top10_sigDE_genes <- sig_DE_gene %>% data.frame() %>%arrange(padj) %>% pull(gene) %>% head(n=10) 
info_top_10_gene<- res_sig %>%filter(padj < p_adj & abs(log2FoldChange) >LFC)  %>%arrange(padj) %>% head(n=10)
DF_to_10_gene<-as.data.frame(info_top_10_gene)
write.csv(DF_to_10_gene, file="DF_TOP_10.csv", row.names = FALSE )

lab_italics <- paste0("italic('", rownames(res_NO_NA), "')")
selectLab_italics = paste0(
  "italic('",
  c(top10_sigDE_genes),
  "')")
png("volcano.png",width= 700, height = 500)
EnhancedVolcano(res_NO_NA,
                lab = lab_italics,
                x = 'log2FoldChange',
                y = 'pvalue',
                title='WT vs Mutated with p_adj',
                selectLab = selectLab_italics,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = p_adj,
                FCcutoff = LFC,
                pointSize = 3.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('purple', 'red3', 'black', 'pink'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()

dev.off()
