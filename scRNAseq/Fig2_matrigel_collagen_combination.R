
# This is for the combination of both datasets (matrigel and collagen)

# needs conda or similar installed for python, miniconda made problems
library(here)
library(Seurat)
library(tidyverse)
library("viridis") 
library(RColorBrewer)
here()



# files need to be original gz files, matrix, barcodes..., project id will be later used for identification of set
matrigel.data<-Read10X(data.dir = here("src/AMC_SC_20_1_matrigel/"))
matrigel <- CreateSeuratObject(counts = matrigel.data, project = "matrigel")
collagen.data<-Read10X(data.dir = here("src/AMC_SC_20_2_collagen/"))
collagen <- CreateSeuratObject(counts = collagen.data, project = "collagen")

class(matrigel.data)
?colSums
# summary of total expression per single cell
summary(colSums(as.matrix(matrigel.data)))
summary(colSums(as.matrix(collagen.data)))

h1<-hist(colSums(as.matrix(matrigel.data)),
         breaks = 100, main = "Expression sum per cell",
         xlab = "Sum expression")
h2<-hist(colSums(as.matrix(collagen.data)),
         breaks = 100)#,add=T)


# Merge more than two Seurat Objects
culture.combined <- merge(matrigel, y = collagen, 
                          add.cell.ids = c("matrigel", "collagen"), project = "intestine_culture")
culture.combined
table(culture.combined$orig.ident)

# from here normal processing of Seurat object as described in tutorial

# QC, normalisation and scaling
# mitochondria genes conveniently start with MT, or for mouse data mt
head(rownames(x=culture.combined@assays$RNA@data))
mito.genes <- grep(pattern = "^mt-", x = rownames(x = culture.combined@assays$RNA@data), value = TRUE, ignore.case = TRUE)
length(mito.genes)
mito.genes

culture.combined[["percent.mt"]] <- PercentageFeatureSet(culture.combined, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(culture.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


FeatureScatter(culture.combined, feature1 = "nFeature_RNA", feature2 = "percent.mt")+
  geom_vline(xintercept = 800)+
  geom_hline(yintercept = 25)
FeatureScatter(culture.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# manual check for filtering parameters, I use Empirical Cumulative Density Function to quickly 
# see effect of percent.mt on % of cells included
head(culture.combined@meta.data)
culture.combined@meta.data%>%
  ggplot(aes(percent.mt, colour=orig.ident))+
  stat_ecdf(geom="point")+
  theme_bw()

culture.combined@meta.data%>%
  group_by(orig.ident)%>%
  summarise(n=n(),
            median.mt=median(percent.mt),
            cells_mt25=sum(percent.mt<25),
            cells_mt20=sum(percent.mt<20),
            cells_mt15=sum(percent.mt<15),
            cells_mt10=sum(percent.mt<10),
            cells_mt05=sum(percent.mt<5))



# subsetting the data, 1st run as permissive as possible
culture.combined
culture.combined<-subset(culture.combined, subset = percent.mt<20 &
                           nFeature_RNA >700)
culture.combined

# ------- Normalizing the data -------
before_norm<-hist(colSums(as.matrix(culture.combined@assays$RNA@data)),
                  breaks = 100,
                  main = "Total expression before normalisation",
                  
                  xlab = "Sum of expression")
culture.combined <- NormalizeData(culture.combined, normalization.method = "LogNormalize", scale.factor = 10000)


after_norm<-hist(colSums(as.matrix(culture.combined@assays$RNA@data)),
                 breaks = 100,
                 main = "Total expression after normalisation",
                 xlab = "Sum of expression")

# ----- Identification of highly variable features (feature selection)

culture.combined <- FindVariableFeatures(culture.combined, selection.method = "vst", nfeatures = 2000)
culture.combined

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(culture.combined), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(culture.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


# ----- Scaling and PCA ------

all.genes <- rownames(culture.combined)

# also try with removing unwanted sources of heterogeneity, e.g. mt-contamination, cell-cycle etc
# culture.combined <- ScaleData(culture.combined, vars.to.regress = "percent.mt")
culture.combined <- ScaleData(culture.combined, features = all.genes)

# PCA
culture.combined <- RunPCA(culture.combined, features = VariableFeatures(object = culture.combined))
DimPlot(culture.combined, reduction = "pca")


# Determine dimensionality of the dataset
ElbowPlot(culture.combined)

DimHeatmap(culture.combined, dims = 1:6, cells = 500, balanced = TRUE)

# -------- Cluster the cells --------
# try with cluster 1:6 for PCA !
# also try resolution 0.3 - 0.9

set.seed(1)

culture.combined <- FindNeighbors(culture.combined, dims = 1:9)
culture.combined <- FindClusters(culture.combined, resolution = 0.5)

# UMAP/tSNE
# culture.combined<-RunUMAP(culture.combined, dims = 1:9)
# DimPlot(culture.combined, reduction = "umap",
#          pt.size = 1, label = TRUE)
# DimPlot(culture.combined, reduction = "umap",
#          pt.size = 1, group.by = "orig.ident")

culture.combined<-RunTSNE(culture.combined, dims = 1:9)
DimPlot(culture.combined, reduction = "tsne", pt.size = 1, label = TRUE)
DimPlot(culture.combined, reduction = "tsne", pt.size = 1,
        group.by = "orig.ident", label = TRUE)
DimPlot(culture.combined, reduction = "tsne", pt.size = 1,
        label = TRUE, split.by = "orig.ident")

grey1<-brewer.pal(n=9, "Greys")[4]
grey2<-brewer.pal(n=9, "Greys")[7]
oranges<-brewer.pal(n=9, "OrRd")[c(5, 7, 8)]
combi_tSNE_colors<-c(grey1, grey2, grey2, grey1, grey2, grey1, grey1, oranges[1], grey2, grey1, oranges[2], oranges[3])
DimPlot(culture.combined, reduction = "tsne", pt.size = 1,
        label = TRUE, cols = combi_tSNE_colors)
DimPlot(culture.combined, reduction = "tsne", pt.size = 1,
        label = TRUE, cols = combi_tSNE_colors, split.by = 'orig.ident')+
  scale_color_discrete(name="Cluster Id",
                     type = combi_tSNE_colors,
                     breaks=as.character(seq(0,11)),
                     labels=c("0 - enterocytes 1 (collagen)",
                              "1 - Lgr5+ CBCs (matrigel)",
                              "2 - TA cells (matrigel)",
                              "3 - enterocytes 2 (collagen)",
                              "4 - late enterocytes (matrigel)",
                              "5 - Lgr5 (collagen)",
                              "6 - early enterocytes (collagen)",
                              "7 - Paneth and goblet cells",
                              "8 - early enterocytes (matrigel)",
                              "9 - Ly6a+/Clu+ cells (collagen)",
                              "10 - enteroendocrine cells",
                              "11 - Tuft cells"))+
xlab("tSNE 1")+ylab("tSNE 2")




cluster.names<-c("0 - enterocytes 1 (collagen)",
                 "1 - Lgr5 (matrigel)",
                 "2 - TA cells (matrigel)",
                 "3 - enterocytes 2 (collagen)",
                 "4 - late enterocytes (matrigel)",
                 "5 - Lgr5 (collagen)",
                 "6 - enterocytes 3 (collagen)",
                 "7 - Paneth and goblet cells",
                 "8 - early enterocytes (matrigel)",
                 "9 - fetal cells (collagen)",
                 "10 - enteroendocrine cells",
                 "11 - Tuft cells")
