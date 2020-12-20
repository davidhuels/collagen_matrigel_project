

# This is the final R code for publication figures
# Matrigel cells separately
# several chunks of code from Dave Tang's blog

library(here)
library(Seurat)
library(tidyverse)
#install.packages("viridis")
library("viridis") 
library(RColorBrewer)

# files need to be original gz files, matrix, barcodes..., project id will be later used for identification of set
matrigel.data<-Read10X(data.dir = here("src/AMC_SC_20_1_matrigel/"))
matrigel <- CreateSeuratObject(counts = matrigel.data, project = "matrigel")


# summary of total expression per single cell
summary(colSums(as.matrix(matrigel.data)))


h1<-hist(colSums(as.matrix(matrigel.data)),
         breaks = 100, main = "Expression sum per cell",
         xlab = "Sum expression")


# QC, normalisation and scaling
# mitochondria genes conveniently start with MT, or for mouse data mt
head(rownames(x=matrigel@assays$RNA@data))
mito.genes <- grep(pattern = "^mt-", x = rownames(x = matrigel@assays$RNA@data), value = TRUE, ignore.case = TRUE)
length(mito.genes)
mito.genes

matrigel[["percent.mt"]] <- PercentageFeatureSet(matrigel, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(matrigel, features = c("nFeature_RNA", "percent.mt"), ncol = 2)


FeatureScatter(matrigel, feature1 = "nFeature_RNA", feature2 = "percent.mt")+
  geom_hline(yintercept = 15)+
  geom_vline(xintercept = 1000)
FeatureScatter(matrigel, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# manual check for filtering parameters, using Empirical Cumulative Density Function to quickly 
# see effect of percent.mt on % of cells included
head(matrigel@meta.data)
matrigel@meta.data%>%
  ggplot(aes(percent.mt, colour=orig.ident))+
  stat_ecdf(geom="point")+
  theme_bw()

matrigel@meta.data%>%
  group_by(orig.ident)%>%
  summarise(n=n(),
            median.mt=median(percent.mt),
            cells_mt25=sum(percent.mt<25),
            cells_mt20=sum(percent.mt<20),
            cells_mt15=sum(percent.mt<15),
            cells_mt10=sum(percent.mt<10),
            cells_mt05=sum(percent.mt<5))


# subsetting the data, with mt-RNA cut-off
matrigel<-subset(matrigel, subset = percent.mt < 15 & nFeature_RNA > 1000)


# ------- Normalizing the data -------
before_norm<-hist(colSums(as.matrix(matrigel@assays$RNA@data)),
                  breaks = 100,
                  main = "Total expression before normalisation",
                  xlab = "Sum of expression")
matrigel <- NormalizeData(matrigel, normalization.method = "LogNormalize", scale.factor = 10000)

after_norm<-hist(colSums(as.matrix(matrigel@assays$RNA@data)),
                 breaks = 100,
                 main = "Total expression after normalisation",
                 xlab = "Sum of expression")

# how many genes are expressed per cell
hist(colSums(as.matrix(matrigel@assays$RNA@data)!=0),
     main = "How many genes are expressed per cell",
     xlab = "expressed genes")

# ----- Identification of highly variable features (feature selection)

matrigel <- FindVariableFeatures(matrigel, selection.method = "vst", nfeatures = 2000)
matrigel

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(matrigel), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(matrigel)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


# ----- Scaling and PCA ------

all.genes <- rownames(matrigel)
# also try with removing unwanted sources of heterogeneity, e.g. mt-contamination, cell-cycle etc
# culture.combined <- ScaleData(culture.combined, vars.to.regress = "percent.mt")
matrigel <- ScaleData(matrigel, features = all.genes)

# PCA
matrigel <- RunPCA(matrigel, features = VariableFeatures(object = matrigel))
DimPlot(matrigel, reduction = "pca")

# Determine dimensionality of the dataset
ElbowPlot(matrigel)

DimHeatmap(matrigel, dims = 1:8, cells = 500, balanced = TRUE)

# -------- Cluster the cells --------
# according to Elbow plot try with dims 1:5 or more!
# also try resolution 0.3 - 0.9
matrigel <- FindNeighbors(matrigel, dims = 1:5)
matrigel <- FindClusters(matrigel, resolution = 0.5)

# UMAP/tSNE
#matrigel<-RunUMAP(matrigel, dims = 1:5)
#DimPlot(matrigel, reduction = "umap",pt.size = 1.5, label = TRUE)

matrigel<-RunTSNE(matrigel, dims = 1:5)
DimPlot(matrigel, reduction = "tsne", pt.size = 1.5, label = TRUE)



# find markers for every cluster compared to all remaining cells, report only the positive ones

matrigel.markers <- FindAllMarkers(matrigel, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(matrigel.markers)
matrigel.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_logFC)
matrigel_cluster_top100<-matrigel.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_logFC)

#write.csv(matrigel_cluster_top100, file="matrigel_cluster_top100.csv")



# gene signatures for cell types
stem<-c("Lgr5", "Ascl2", "Slc12a2", "Axin2", "Olfm4", "Gkn3")
cell_cycle<-c("Mki67", "Cdk4", "Mcm5", "Mcm6", "Pcna")
enterocyte<-c("Alpi", "Apoa1", "Apoa4", "Fabp1", "Aldob")
goblet<-c("Muc2", "Clca3", "Tff3", "Agr2")
paneth<-c("Lyz1", "Defa17", "Defa22", "Defa24", "Ang4")
tuft<-c("Dclk1", "Trpm5", "Gfi1b", "Il25")
enteroendo<-c("Chga", "Chgb", "Tac1", "Tph1", "Neurog3")

cell_types<-c("Lgr5", "Olfm4", "Mki67", "Ccnd2", "Apoa1", "Muc2", "Lyz1", "Dclk1", "Chga")

colours_tSNE<-c("#d82e2e", "#884664", "#4876a4", "#498f8a", "#53a361",
                "#d97141", "#f79a4b", "#f5c643", "#bf9041", "#a25b3c", "#e684b5")

### Figure of cluster characteristics
VlnPlot(object = matrigel, features = cell_types, stack = TRUE, flip = TRUE,
        cols = colours_tSNE, fill.by = 'ident')+
  NoLegend()



DimPlot(matrigel, reduction = "tsne", pt.size = 1.5, label = TRUE)
# change the order of numbering, also for later violin plot and bar plot
# rename clusters with assigned cell-type, we have 9 clusters
#new.cluster.ids<-c(0, 3, 1, 5, 2, 4, 8, 6, 7)

Idents(matrigel)
levels(x = matrigel)
matrigel <- RenameIdents(object = matrigel, 
                           '0' = '0',
                           '3' = '1',
                           '1' = '2',
                           '5' = '3',
                           '2' = '4',
                           '4' = '5',
                           '8' = '6',
                           '6' = '7',
                           '7' = '8')
levels(x = matrigel)
levels(matrigel@active.ident)


DimPlot(matrigel, reduction = "tsne", pt.size = 1.5, label = TRUE)
#old.cluster.ids<-as.character(seq(0,8))
levels(matrigel)

FeaturePlot(object = matrigel,
            features = enterocyte,
            reduction= "tsne")&
  scale_color_gradientn(colors = viridis(n=10))
# --------------------- Figures for publication ---------------------------------

dp<-DimPlot(matrigel, reduction = "tsne", pt.size = 1.5,
        label = FALSE, cols = colours_tSNE)+
  scale_color_discrete(name="Cluster Id",
                       type = colours_tSNE,
                       breaks=as.character(seq(0,8)),
                       labels=c("0 - Lgr5 CBCs",
                                         "1 - early TA cells",
                                         "2 - late TA cells",
                                         "3 - enterocytes, early progenitors",
                                         "4 - enterocytes, late progenitors",
                                         "5 - enterocytes",
                                         "6 - secretory cell progenitors",
                                         "7 - Paneth and goblet cells",
                                         "8 - enteroendocrine cells"))+
  xlab("tSNE 1")+ylab("tSNE 2")
  
LabelClusters(dp, id="ident", size=5, box = TRUE, fill="white")

### Figure of cluster characteristics
VlnPlot(object = matrigel, features = cell_types, stack = TRUE, flip = TRUE,
        cols = colours_tSNE, fill.by = 'ident')+
  NoLegend()

## check for cell type composition
# find percentage of clusters per group


cluster_freq<-data.frame(table(matrigel@active.ident))
cluster_freq$relfreq<-cluster_freq$Freq/sum(cluster_freq$Freq)

cluster_freq$Var1<-factor(cluster_freq$Var1, levels = seq(8,0))
cluster_freq%>%
ggplot(aes(x=Var1, y=relfreq, fill=Var1))+
  geom_col(position = "dodge", fill = colours_tSNE[1:9])+ 
  xlab("cluster")+ylab("Frequency")+
  scale_x_discrete(breaks=seq(0,8),
                   labels=c("0 - Lgr5 CBCs",
                            "1 - early TA cells",
                            "2 - late TA cells",
                            "3 - enterocytes, early progenitors",
                            "4 - enterocytes, late progenitors",
                            "5 - enterocytes",
                            "6 - secretory cell progenitors",
                            "7 - Paneth and goblet cells",
                            "8 - enteroendocrine cells"))+
  coord_flip()+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text = element_text(face = 'bold',
                                 colour = 'black',
                                 size = 10))



