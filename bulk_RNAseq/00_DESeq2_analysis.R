# RNAseq R1.fastq.files were aligned to mm10primay assembly using Rsubread and featureCounts
# counts were directly extracted from 
# fc<-featureCounts(paste0("bam/",all.bam), annot.inbuilt="mm10")



# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# BiocManager::install("DESeq2")
# BiocManager::install("org.Mm.eg.db")
# 
# install.packages("tidyverse")
# install.packages("here")
# install.packages("pheatmap")
# install.packages("RColorBrewer")
# install.packages("ggrepel")

library(DESeq2)
library(tidyverse)
library(here)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(org.Mm.eg.db)
library(fgsea)

# the rownames(or here 1st column) from featurecounts is Entrez gene identifiers
cts<-read.csv(here("output_files/wt_collagen_matrigel_RNAseq_counts.csv"), stringsAsFactors = FALSE)
cts[1:5, 1:5]
dim(cts)
rownames(cts)<-cts$X
cts<-cts[,-1]
dim(cts)
columns<-sapply(colnames(cts), function(x) {
  unlist(strsplit(x, "[.]"))[1]
})
colnames(cts)<-as.character(columns)
cts<-as.matrix(cts)


# how many genes are expressed in each sample, more than 0 or 10 reads...
colSums(cts>0)
colSums(cts>10)

# get sample information data
samples<-read.csv(here("sample_information.csv"), stringsAsFactors = FALSE)
head(samples)
samples<-samples[-24,]
seq(1:13)
tail(samples)
samples$RNAseq_id<-paste0("R", seq(1:23))
rownames(samples)<-samples$RNAseq_id
samples<-samples[,-1]

head(samples)
# important, get the count matrix in the order of the sample rownames!!!
head(cts)
cts <- cts[, rownames(samples)]


#check the order of sample_ids and columns in count matrix
all(rownames(samples) %in% colnames(cts))
all(rownames(samples) == colnames(cts))

samples$group<-paste(samples[,3], samples[,4], samples[,5], sep = "_")
unique(samples$group)
### create DESeqDataSet (dds)

dds<-DESeqDataSetFromMatrix(countData = cts,
                            colData = samples,
                            design = ~ group)
dds

# Pre-filtering, removing genes with low reads
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]
dds
#saveRDS(dds, here("output_files/dds_LgrApc_invitro_structure.rds"))
#dds<-readRDS(here("output_files/dds_LgrApc_invitro_structure.rds"))
class(dds$group)

# set the factor levels straight for contrast build
# for strucure already 2D vs 3D
# dds$structure<-factor(dds$structure, levels = c("2D", "3D"))

### Data quality assessment by sample clustering and visualization --------------------------

# Extracting transformed values
# vst = variance stabilizing transformation
# The transformed data is on the log2 scale for large counts
vsd<-vst(dds, blind = FALSE)
vsd_wt<-assay(vsd)
vsd_si<-vsd[,vsd$tissue == "Si"]
#write.csv(as.data.frame(vsd_invitro), here("output_files", "vst_invitro.csv"))

# Heatmap of the sample-to-sample distance
# dist creates 12x12 matrix with distances between samples


sampleDist<-dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)
sampleDistMatrix
rownames(sampleDistMatrix) <- vsd$group

colnames(sampleDistMatrix) <- paste(colnames(sampleDistMatrix), vsd$treatment, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDist,
         clustering_distance_cols=sampleDist,
         col=colors,
         main="Distance Matrix")

# PCA plots
pcaData<-plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar<-round(100*attr(pcaData, "percentVar"))
cols<-brewer.pal(8, "Dark2")
#cols<-c("#999999", "#E69F00", "#56B4E9", "#009E73")
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=5) +
  geom_point(shape=1, size=5, colour="black")+
  #geom_text_repel(aes(label=group), nudge_y = 5)+
  scale_color_manual(values=cols)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()+
  theme(panel.grid = element_blank())

# splitting the deseqdatafile for Si and colon
vsd_si<-vsd[,vsd$tissue == "Si"]

pcaData_si<-plotPCA(vsd_si, intgroup="group", returnData=TRUE)
percentVar_si<-round(100*attr(pcaData_si, "percentVar"))
cols<-brewer.pal(8, "Dark2")
#cols<-c("#999999", "#E69F00", "#56B4E9", "#009E73")
ggplot(pcaData_si, aes(PC1, PC2, color=group)) +
  geom_point(size=5) +
  geom_point(shape=1, size=5, colour="black")+
  #geom_text_repel(aes(label=group), nudge_y = 5)+
  scale_color_manual(values=cols)+
  xlab(paste0("PC1: ",percentVar_si[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_si[2],"% variance")) + 
  coord_fixed()+
  theme_bw()+
  theme(panel.grid = element_blank())

# save the file for additional analyses...so no need to do the 1st QCs again
dds<-DESeq(dds)
#saveRDS(dds, file = "RNAseq_all_samples_groups.rds")

### Differential expression anaysis -------------------------------------------------------
# here with FDR <=0.05, could also use FDR <=0.1






resultsNames(dds)
res<- results(dds, alpha = 0.05, contrast = c("group", "Si_ENR_onColl", "Si_ENR_inMG"))

# if a row has a low mean normalized count (baseMean) then p-adj will be NA...remove
res<-res[!is.na(res$padj),]

rownames(res)[1:5]
res$symbol<-mapIds(org.Mm.eg.db, rownames(res),
                        keytype = "ENTREZID",
                        column = "SYMBOL")
# remove all non-matched entrez ids...lncRNA etc.
res<-res[!is.na(res$symbol),]
dim(res)
table(res$padj<0.05)
summary(res)
head(res)

head(vsd_wt)

colnames(vsd_wt)<-samples$group
all(rownames(res) %in% rownames(vsd_wt))

res_3D2D<-merge(res, counts_normal, by="row.names")
head(res_3D2D)
dim(res_3D2D)

#write.csv(res_3D2D, here("output_files/res_3D_2D_normCounts.csv"))

### fgsea with MSigDB and and yap genes--------------------------------------------------
#load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"))
# This will load Mm.c2, which is a list of gene sets, each a vector of Entrez Ids. 

#saveRDS(Mm.c2, file = "MSigDB_mouse_C2_v5p2.rdata")
Mm.c2<-readRDS(file="MSigDB_mouse_C2_v5p2.rdata")
length(Mm.c2)

# make fgsea with Entrez Ids

res$entrez<-rownames(res)

res_ranks<-res %>% 
  as.data.frame()%>%
  dplyr::select(entrez, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(entrez) %>% 
  summarize(stat=mean(stat))
dim(res_ranks)
# creating a named, numeric vector with gene list
ranks <- res_ranks[["stat"]]
names(ranks)<-res_ranks[["entrez"]]
head(ranks)

# Mustata fetal signature
fetal<-read.csv(here("mustata_fetal_genes.csv"), stringsAsFactors = FALSE)
head(fetal)
colnames(fetal)<-c("log2FC", "symbol", "genename")
fetal$entrez<-mapIds(org.Mm.eg.db, fetal$symbol,
                   keytype = "SYMBOL",
                   column = "ENTREZID")
dim(fetal)
fetal<-fetal[!is.na(fetal$entrez),]
head(fetal)
mustata_fetal_up<-fetal[fetal$log2FC>1, "entrez"]
mustata_fetal_down<-fetal[fetal$log2FC<1, "entrez"]
mustata_fetal<-list("mustata_fetal_up"=mustata_fetal_up,
                    "mustata_fetal_down"=mustata_fetal_down)

#Yap signature from Gregorieff
yapgs<-read.csv(here("gregorieff_nature2015_yap.csv"), stringsAsFactors = FALSE)
head(yapgs)
yapgs$entrez <- mapIds(org.Mm.eg.db, yapgs$GeneID, keytype="ENSEMBL", column="ENTREZID")
head(yapgs)
yapgs<-yapgs[!is.na(yapgs$entrez),]

yap_suppressed<-yapgs[yapgs$Combined.fold.change>4,"entrez"]
yap_activated<-yapgs[yapgs$Combined.fold.change<0.2,"entrez"]
length(yap_activated)

gregorieff_yap<-list("gregorieff_yap_suppressed"=yap_suppressed,
                     "gregorieff_yap_activated"=yap_activated)

length(gregorieff_yap)

# fetal gene signatures
si_dev<-read.csv(here("Si_development_top100_per_group.csv"), stringsAsFactors = FALSE)
head(si_dev)
si_dev[si_dev$group=="E12_14", "entrez"]
si_dev$entrez<-mapIds(org.Mm.eg.db, si_dev$symbol,
                     keytype = "SYMBOL",
                     column = "ENTREZID")
si_dev_list<-list("E12_E14_top100_genes"=si_dev[si_dev$group=="E12_14", "entrez"],
                  "E16_top100_genes"=si_dev[si_dev$group=="E16", "entrez"],
                  "E18_top100_genes"=si_dev[si_dev$group=="E18", "entrez"],
                  "adult_top100_genes"=si_dev[si_dev$group=="adult", "entrez"])

str(si_dev_list)

MyPathways<-c(Mm.c2, gregorieff_yap, mustata_fetal, si_dev_list)
length(MyPathways)


fgseaRes <- fgsea(pathways = MyPathways, 
                  stats = ranks,
                  minSize = 15,
                  maxSize = 850,
                  nperm=100000)
dim(fgseaRes)
fgseaResTidy <- fgseaRes %>%
  filter(padj < 0.05)%>%
  as_tibble() %>%
  arrange(desc(abs(NES)))
fgseaResTidy

fgseaRes%>%
  filter(str_detect(pathway, pattern = "mustata|yap"))%>%
  as_tibble()%>%
  print(n=Inf)

fgseaRes %>%
  filter(padj < 0.05 & NES > abs(2.8))%>%
  as_tibble()%>%
  arrange(desc(NES))
  
# extract leadingEdge from inflammatory pathways to check which genes are shared
cell_junction<-fgseaRes%>%
  filter(str_detect(pathway, "REACTOME_CELL_JUNCTION_ORGANIZATION"))%>%
  pull(leadingEdge)

head(cell_junction)
norm_counts<-as.data.frame(assay(vsd_si))
head(norm_counts)
vsd_si$group

colnames(norm_counts)<-vsd_si$group
head(unlist(cell_junction))

z_score <- function(x){
  (x - mean(x)) / sd(x)
}

norm_counts%>%
  rownames_to_column('entrez')%>%
  filter(entrez %in% unlist(cell_junction))%>%
  mutate(symbol = mapIds(org.Mm.eg.db, entrez, keytype="ENTREZID", column="SYMBOL"))%>%
  column_to_rownames('symbol')%>%
  dplyr::select(-entrez)%>%
  as.matrix()%>%
  apply(1, z_score)%>%
  t()%>%
  pheatmap(clustering_method = "ward.D2")

norm_counts%>%
  rownames_to_column('entrez')%>%
  filter(entrez %in% unlist(cell_junction))%>%
  mutate(symbol = mapIds(org.Mm.eg.db, entrez, keytype="ENTREZID", column="SYMBOL"))%>%
  column_to_rownames('symbol')%>%
  dplyr::select(-entrez)%>%
  as.matrix()%>%
  pheatmap(clustering_method = "ward.D2")





# the top100 genes per timepoint didn't work very well, only few genes in that list are differentially expressed 
# between matrigel and collagen cultures


fgseaResTidy%>%
  filter(abs(NES)>2.8)%>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  scale_fill_gradient2()+
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="significant MSigDB pathways \nfrom fGSEA (padj<0.05 & |NES|> 2.8)") + 
  theme_minimal()+
  theme(axis.text =element_text(size=10, colour = "black"))


#----  pheatmap of fetal signature across all samples -----
# Mustata fetal signature
fetal<-read.csv(here("mustata_fetal_genes.csv"), stringsAsFactors = FALSE)
head(fetal)
fetal$Gene.Symbol[2]<-"Gja1"
fetal$Gene.Symbol[5]<-"Tacstd2"
colnames(fetal)<-c("log2FC", "symbol", "genename")
fetal$entrez<-mapIds(org.Mm.eg.db, fetal$symbol,
                     keytype = "SYMBOL",
                     column = "ENTREZID")
dim(fetal)
fetal<-fetal[!is.na(fetal$entrez),]
head(fetal)
mustata_fetal_up<-fetal[fetal$log2FC>1, "entrez"]
mustata_fetal_down<-fetal[fetal$log2FC<1, "entrez"]


z_score <- function(x){
  (x - mean(x)) / sd(x)
}

norm_counts<-as.data.frame(assay(vsd))
head(norm_counts)
vsd$group
all(row.names(colData(dds))==colnames(norm_counts))
colnames(norm_counts)<-vsd$group
head(fetal)
norm_counts%>%
  rownames_to_column('entrez')%>%
  filter(entrez %in% mustata_fetal_up)%>%
  column_to_rownames('entrez')%>%
  as.matrix()%>%
  apply(1, z_score)%>%
  t()%>%
  pheatmap(clustering_method = "ward.D2")
dim(fetal)


# get the most variable genes, calculate sd for each gene in the matrix
fetal_matrix<-norm_counts%>%
  rownames_to_column('entrez')%>%
  filter(entrez %in% fetal$entrez)%>%
  column_to_rownames('entrez')%>%
  as.matrix()

sd_fetalgenes<-apply(fetal_matrix, 1, sd)
head(sd_fetalgenes)
# heatmap of most varible genes (sd>1.5) across samples

# for labels of fetal genes
tail(fetal)
fetal_gene_label <- data.frame(cluster = ifelse(fetal$log2FC > 1, "up fetal organoids",  "up adult organoids"))
rownames(fetal_gene_label)<-fetal$symbol
head(fetal_gene_label)

colnames(norm_counts)

sample_col_label <- data.frame(sample = ifelse(grepl(pattern="MG", colnames(norm_counts)), "matrigel", "collagen"))
sample_col_label
row.names(sample_col_label) <-make.unique(colnames(norm_counts))

norm_counts%>%
  rownames_to_column('entrez')%>%
  filter(entrez %in% names(sd_fetalgenes[sd_fetalgenes>1.5]))%>%
  left_join(fetal[,c(2,4)], by="entrez")%>%
  dplyr::select(-entrez)%>%
  column_to_rownames('symbol')%>%
  as.matrix()%>%
  apply(1, z_score)%>%
  t()%>%
  pheatmap(clustering_method = "ward.D2",
           annotation_row = fetal_gene_label,
           annotation_col = sample_col_label,
           main = "Fetal Genes (Mustata) in Small Intestinal samples\n - highly variable genes (sd > 1.5)")

sd_fetal_up<-sd_fetalgenes[names(sd_fetalgenes) %in% mustata_fetal_up]
length(sd_fetalgenes)
length(sd_fetal_up)
head(sd_fetal_up)
class(mustata_fetal_up[1])
class(names(sd_fetalgenes)[1])
norm_counts%>%
  rownames_to_column('entrez')%>%
  filter(entrez %in% names(sd_fetal_up[sd_fetal_up>1.5]))%>%
  left_join(fetal[,c(2,4)], by="entrez")%>%
  dplyr::select(-entrez)%>%
  column_to_rownames('symbol')%>%
  as.matrix()%>%
  apply(1, z_score)%>%
  t()%>%
  pheatmap(clustering_method = "ward.D2",
           annotation_row = fetal_gene_label,
           annotation_col = sample_col_label,
           main = "Fetal Genes (Mustata) in Small Intestinal samples\n - highly variable genes (sd > 1.5)")
