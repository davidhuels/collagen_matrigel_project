### Analysis of small intestine during development ----


#1 mustata fetal signature...when highly expressed
#2 extracellular matrix changes over time


library(tidyverse)
library(here)
library(DESeq2)
library(viridis)
#library(apeglm)
cts<-read.delim(here("GSE115541_RNAseq_RawRead_Counts_All_Samples.txt"))
head(cts)
colnames(cts)
cts_wt<-cts%>%
  select(-contains("Cdx2"))%>%
  select(Genes, matches("SmallIntestine|Ileum"))


# remove lines with no name
cts_wt<-cts_wt%>%
  filter(Genes !="")
dim(cts_wt)
# 2 rows contain NA, just remove these
dim(cts_wt[complete.cases(cts_wt),])
cts_wt<-cts_wt[complete.cases(cts_wt),]
head(cts_wt)
table(duplicated(cts_wt$Genes))
rownames(cts_wt)<-cts_wt$Genes
cts_m<-as.matrix(cts_wt[,-1])
head(cts_m)
# how many genes are expressed in each sample, more than 0 or 10 reads...
colSums(cts_m>0)
colSums(cts_m>10)

# seems these are really the raw read counts...need to feed it to DESeq2 process

# make sample data.frame
samples<-data.frame(group=str_split(colnames(cts_m), pattern="_", simplify = TRUE)[,1])
rownames(samples)<-colnames(cts_m)
samples
#check the order of sample_ids and columns in count matrix
all(rownames(samples) %in% colnames(cts_m))
all(rownames(samples) == colnames(cts_m))


### create DESeqDataSet (dds)

dds<-DESeqDataSetFromMatrix(countData = cts_m,
                            colData = samples,
                            design = ~ group)
dds

# Pre-filtering, removing genes with low reads
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]
levels(dds$group)[1]

# for DE analysis set factor levels in correct order
# set the factor levels straight for contrast build

dds$group<-factor(dds$group, levels = c(levels(dds$group)[2],
                                        levels(dds$group)[3],
                                        levels(dds$group)[4],
                                        levels(dds$group)[5],
                                        levels(dds$group)[1]))
dds$group
# Extracting transformed values
# vst = variance stabilizing transformation
# The transformed data is on the log2 scale for large counts
vsd<-vst(dds, blind = FALSE)
vsd_embryo<-assay(vsd)

# Heatmap of the sample-to-sample distance
# dist creates 12x12 matrix with distances between samples
sampleDist<-dist(t(assay(vsd)))


library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDist)
sampleDistMatrix
rownames(sampleDistMatrix) <- vsd$group

colnames(sampleDistMatrix) <- paste(colnames(sampleDistMatrix), vsd$group, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDist,
         clustering_distance_cols=sampleDist,
         col=colors,
         main="Distance Matrix")

# PCA plots

pcaData<-plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar<-round(100*attr(pcaData, "percentVar"))
cols<-brewer.pal(5, "Dark2")
cols
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=5) +
  geom_point(shape=1, size=5, colour="black")+
  scale_color_manual(values=cols)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()+
  theme(panel.grid = element_blank())


# Heatmap of Mustata fetal signature
# Mustata fetal signature



library(org.Mm.eg.db)
fetal<-read.csv(here("mustata_fetal_genes.csv"), stringsAsFactors = FALSE)
head(fetal)
colnames(fetal)<-c("log2FC", "symbol", "genename")
fetal$entrez<-mapIds(org.Mm.eg.db, fetal$symbol,
                     keytype = "SYMBOL",
                     column = "ENTREZID")
dim(fetal)
fetal<-fetal[!is.na(fetal$entrez),]
head(fetal)
# first clustering based on the subset of genes
#install.packages("dendextend")
library(dendextend)
head(cts_m)
head(vsd_embryo)
vsd_embryo_df<-as.data.frame(vsd_embryo)
head(vsd_embryo_df)
vsd_embryo_df$symbol<-rownames(vsd_embryo_df)
head(fetal)
tail(fetal)

mustata_fetal_up<-fetal[fetal$log2FC>1, "symbol"]
mustata_fetal_down<-fetal[fetal$log2FC<1, "symbol"]

fetal_matrix_up<-vsd_embryo_df%>%
  filter(symbol %in% mustata_fetal_up)%>%
  dplyr::select(-symbol)%>%
  as.matrix()
rownames(fetal_matrix_up)<-vsd_embryo_df%>%
  filter(symbol %in% mustata_fetal_up)%>%
  pull(symbol)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}


fetal_norm_up <- t(apply(fetal_matrix_up, 1, cal_z_score))

pheatmap(fetal_norm_up, clustering_method = "ward.D2",
         cutree_rows = 3, show_rownames = FALSE)
pheatmap(fetal_matrix_up, clustering_method = "ward.D2")


fetal_matrix_down<-vsd_embryo_df%>%
  filter(symbol %in% mustata_fetal_down)%>%
  dplyr::select(-symbol)%>%
  as.matrix()
rownames(fetal_matrix_down)<-vsd_embryo_df%>%
  filter(symbol %in% mustata_fetal_down)%>%
  pull(symbol)

fetal_norm_down <- t(apply(fetal_matrix_down, 1, cal_z_score))
pheatmap(fetal_norm_down, clustering_method = "ward.D2")

library(dendextend)
d<-dist(fetal_norm_up)
hc<-hclust(d, method = "ward.D2")
plot(as.dendrogram(hc))


# check expression for selected genes
selected.genes<-c("Lgr5", "Olfm4", "Lyz1", "Lyz2", 
                  "Clu", "Msln", "Ly6a", "Col4a1", "Col4a2", 
                  "Itga6", "Itgb4", "Itgb1")
selected.genes
selected.genes.df<-data.frame(selected.genes, stringsAsFactors = FALSE)
head(selected.genes.df)
class(selected.genes.df[1,1])
selected.genes.df$entrez<-mapIds(org.Mm.eg.db, selected.genes.df$selected.genes,
                     keytype = "SYMBOL",
                     column = "ENTREZID")
selected.genes.df


head(vsd_embryo)
vsd_embryo_df<-as.data.frame(vsd_embryo)
vsd_embryo_df$symbol<-rownames(vsd_embryo_df)


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

vsd_embryo_df%>%
  filter(symbol %in% selected.genes.df$selected.genes)%>%
  dplyr::select(-symbol)%>%
  as.matrix()%>%
  apply(1, cal_z_score)%>%
  t()%>%
  pheatmap(clustering_method = "ward.D2", 
           cluster_cols = FALSE,
           main = "Selected Genes during development\n(rows scaled)")

vsd_embryo_df%>%
  filter(symbol %in% selected.genes.df$selected.genes)%>%
  dplyr::select(-symbol)%>%
  as.matrix()%>%
  # apply(1, cal_z_score)%>%
  # t()%>%
  pheatmap(clustering_method = "ward.D2", 
           cluster_cols = FALSE)

# all integrin receptors and Dystroglycan (other Laminin receptor) also include main non-integrin receptor Rpsa
vsd_embryo_df%>%
  filter(str_detect(symbol, "Itga|Itgb|Dag|Rpsa"))%>%
  dplyr::select(-symbol)%>%
  pheatmap(clustering_method = "ward.D2",
           cluster_cols = FALSE,
           main = "Integrin and Dystroglycan receptors during development\n(rows not scaled)")

cluster1<-c("Dag1", "Itgb1", "Itgb4", "Itga3", "Itga6")
vsd_embryo_df%>%
  filter(symbol %in% cluster1)%>%
  dplyr::select(-symbol)%>%
  pheatmap(clustering_method = "ward.D2",
           cluster_cols = FALSE,
           color = plasma(10))

