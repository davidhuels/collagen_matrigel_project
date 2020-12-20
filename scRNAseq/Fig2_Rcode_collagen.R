
# This is the final R code for publication figures
# Collagen cells separately
# several chunks of code from Dave Tang's blog

library(here)
library(Seurat)
library(tidyverse)
#install.packages("viridis")
library("viridis") 
library(RColorBrewer)

# --------------- collagen analysis ------------------------------------

# files need to be original gz files, matrix, barcodes..., project id will be later used for identification of set
collagen.data<-Read10X(data.dir = here("src/AMC_SC_20_2_collagen/"))
collagen <- CreateSeuratObject(counts = collagen.data, project = "collagen")


# QC, normalisation and scaling
# mitochondria genes conveniently start with MT, or for mouse data mt
head(rownames(x=collagen@assays$RNA@data))
mito.genes <- grep(pattern = "^mt-", x = rownames(x = collagen@assays$RNA@data), value = TRUE, ignore.case = TRUE)
length(mito.genes)
mito.genes

collagen[["percent.mt"]] <- PercentageFeatureSet(collagen, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(collagen, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(collagen, features = c("nFeature_RNA", "percent.mt"), ncol = 2)

FeatureScatter(collagen, feature1 = "nFeature_RNA", feature2 = "percent.mt")+
  geom_hline(yintercept = 15)+
  geom_vline(xintercept = 1000)
FeatureScatter(collagen, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# manual check for filtering parameters, I use Empirical Cumulative Density Function to quickly 
# see effect of percent.mt on % of cells included
head(collagen@meta.data)
collagen@meta.data%>%
  ggplot(aes(percent.mt, colour=orig.ident))+
  stat_ecdf(geom="point")+
  theme_bw()


collagen@meta.data%>%
  group_by(orig.ident)%>%
  summarize(n_perc(percent.mt<0.25))
collagen@meta.data%>%
  group_by(orig.ident)%>%
  summarise(n=n(),
            median.mt=median(percent.mt),
            cells_mt25=sum(percent.mt<25),
            cells_mt20=sum(percent.mt<20),
            cells_mt15=sum(percent.mt<15),
            cells_mt10=sum(percent.mt<10),
            cells_mt05=sum(percent.mt<5))



# subsetting the data, 1st run as permissive as possible
collagen<-subset(collagen, subset = percent.mt < 15 & nFeature_RNA > 1000)


# ------- Normalizing the data -------
before_norm<-hist(colSums(as.matrix(collagen@assays$RNA@data)),
                  breaks = 100,
                  main = "Total expression before normalisation",
                  
                  xlab = "Sum of expression")
collagen <- NormalizeData(collagen, normalization.method = "LogNormalize", scale.factor = 10000)


after_norm<-hist(colSums(as.matrix(collagen@assays$RNA@data)),
                 breaks = 100,
                 main = "Total expression after normalisation",
                 xlab = "Sum of expression")

# ----- Identification of highly variable features (feature selection)

collagen <- FindVariableFeatures(collagen, selection.method = "vst", nfeatures = 2000)
collagen

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(collagen), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(collagen)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


# ----- Scaling and PCA ------

all.genes <- rownames(collagen)
# also try with removing unwanted sources of heterogeneity, e.g. mt-contamination, cell-cycle etc
# culture.combined <- ScaleData(culture.combined, vars.to.regress = "percent.mt")
collagen <- ScaleData(collagen, features = all.genes)

# PCA
collagen <- RunPCA(collagen, features = VariableFeatures(object = collagen))
DimPlot(collagen, reduction = "pca")

# Determine dimensionality of the dataset
ElbowPlot(collagen)

DimHeatmap(collagen, dims = 1:7, cells = 500, balanced = TRUE)

# -------- Cluster the cells --------
# try with cluster 1:6 for PCA !
# also try resolution 0.3 - 0.9
collagen <- FindNeighbors(collagen, dims = 1:6)
collagen <- FindClusters(collagen, resolution = 0.5)


# UMAP/tSNE
#collagen<-RunUMAP(collagen, dims = 1:6)
#DimPlot(collagen, reduction = "umap",pt.size = 1.5, label = TRUE)

collagen<-RunTSNE(collagen, dims = 1:6)
DimPlot(collagen, reduction = "tsne", pt.size = 1.5, label = TRUE)

collagen.markers <- FindAllMarkers(collagen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
collagen.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_logFC) %>%
  print (n=Inf)



cell_types<-c("Lgr5", "Olfm4", "Mki67", "Ccnd2", "Apoa1", "Muc2", "Lyz1", "Dclk1", "Chga", "Ly6a", "Clu")
# for supplementary figure
FeaturePlot(object = collagen,
            features = cell_types,
            reduction= "tsne")&
  scale_color_gradientn(colors = viridis(n=10))
enterocyte<-c("Alpi", "Apoa1", "Apoa4", "Fabp1", "Aldob")
FeaturePlot(object = collagen,
            features = c("Ly6a", "Clu", "Msln", "Col4a1"),
            reduction= "tsne")&
  scale_color_gradientn(colors = viridis(n=10))



# ------------ Figures for publication --------------------------



DimPlot(collagen, reduction = "tsne", pt.size = 1.5, label = TRUE)
levels(x = collagen)
collagen <- RenameIdents(object = collagen, 
                         '7' = '0',
                         '2' = '1',
                         '5' = '2',
                         '3' = '3',
                         '0' = '4',
                         '1' = '5',
                         '6' = '6',
                         '4' = '7')
levels(x = collagen)
DimPlot(collagen, reduction = "tsne", pt.size = 1.5, label = TRUE)


colours_coltSNE<-c("#d82e2e", "#884664", "#4876a4", "#498f8a", "#53a361",
                   "#d97141", "#f79a4b", "#bb3978")

dpc<-DimPlot(collagen, reduction = "tsne", pt.size = 1.5,
             label = FALSE, cols = colours_coltSNE)+
  # scale_color_discrete(name="Cluster Id",
  #                      type = colours_tSNE,
  #                      breaks=as.character(seq(0,7)),
  #                      labels=c("0 - Lgr5 CBCs",
  #                               "1 - early TA cells",
  #                               "2 - late TA cells",
  #                               "3 - enterocytes, early progenitors",
  #                               "4 - enterocytes 1",
  #                               "5 - enterocytes 2",
  #                               "6 - secretory cells",
  #                               "7 - Ly6a+/Clu+ cells"))+
xlab("tSNE 1")+ylab("tSNE 2")

LabelClusters(dpc, id="ident", size=5, box = TRUE, fill="white")

### Figure of cluster characteristics
VlnPlot(object = collagen, features = cell_types, stack = TRUE, flip = TRUE,
        cols = colours_coltSNE, fill.by = 'ident')+
  NoLegend()

## check for cell type composition
# find percentage of clusters per group


cluster_freq<-data.frame(table(collagen@active.ident))
cluster_freq$relfreq<-cluster_freq$Freq/sum(cluster_freq$Freq)
cluster_freq
cluster_freq$Var1<-factor(cluster_freq$Var1, levels = seq(7,0))
cluster_freq%>%
  ggplot(aes(x=Var1, y=relfreq, fill=Var1))+
  geom_col(position = "dodge", fill = colours_coltSNE)+ 
  xlab("cluster")+ylab("Frequency")+
  scale_x_discrete(breaks=seq(0,7),
                   labels=c("0 - Lgr5 CBCs",
                            "1 - early TA cells",
                            "2 - late TA cells",
                            "3 - enterocytes, early progenitors",
                            "4 - enterocytes 1",
                            "5 - enterocytes 2",
                            "6 - secretory cells",
                            "7 - Ly6a+/Clu+ cells"))+
  coord_flip()+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text = element_text(face = 'bold',
                                 colour = 'black',
                                 size = 10))
