#Code by Alex Whitehead
# Feb 2021, contact @ alwhiteh@ucsd.edu

# for Renv installation of bioconductor packages, use renv::install("bioc::Biobase")
# If there is an old version of a package (from an old version of R) that needs to be removed
# you can use Renv::purge('RCurl')
library(Seurat) #for RNA
library(Signac) #for ATAC
library(dplyr)
library(Matrix)
library(scran)
library(ggplot2)
library(sctransform)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)
library(renv)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
#To restore Rdata use:
#load("Olson_Global_Objects.RData")

#renv::init()
renv::snapshot()

# for sc RNA seq you need a directory for each sample that contains the following files
#barcodes.tsv
#features.tsv (previously called genes.tsv)
#matrix.mtx



# scRNA -------------------------------------------------------------------

# Since we have a lot of samples, we will create a vector of paths
rna_1_1MI.data <- Read10X("/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scRNA/P1_1MI")  
rna_1_1S.data <- Read10X("/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scRNA/P1_1S")
rna_1_3MI.data <- Read10X("/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scRNA/P1_3MI")
rna_1_3S.data <- Read10X("/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scRNA/P1_3S")
rna_8_1MI.data <- Read10X("/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scRNA/P8_1MI")
rna_8_1S.data <- Read10X("/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scRNA/P8_1S")
rna_8_3MI.data <- Read10X("/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scRNA/P8_3MI")
rna_8_3S.data <- Read10X("/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scRNA/P8_3S")
 
#Create Seurat Objects for each one
rna_1_1MI<-CreateSeuratObject(rna_1_1MI.data, min.cells = 3, min.features = 200, project = "1_1MI")
rna_1_1S<-CreateSeuratObject(rna_1_1S.data, min.cells = 3, min.features = 200,project = "1_1S")
rna_1_3MI<-CreateSeuratObject(rna_1_3MI.data,min.cells = 3, min.features = 200,project = "1_3MI")
rna_1_3S<-CreateSeuratObject(rna_1_3S.data,min.cells = 3, min.features = 200,project = "1_3S")
rna_8_1MI<-CreateSeuratObject(rna_8_1MI.data,min.cells = 3, min.features = 200,project = "8_1MI")
rna_8_1S<-CreateSeuratObject(rna_8_1S.data,min.cells = 3, min.features = 200,project = "8_1S")
rna_8_3MI<-CreateSeuratObject(rna_8_3MI.data,min.cells = 3, min.features = 200,project = "8_3MI")
rna_8_3S<-CreateSeuratObject(rna_8_3S.data,min.cells = 3, min.features = 200,project = "8_3S")

# Merge the datasets
All_rna <- merge(rna_1_1MI, y = c(rna_1_1S,rna_1_3MI,rna_1_3S,rna_8_1MI,rna_8_1S,rna_8_3MI,rna_8_3S), 
                 add.cell.ids = c("1_1MI","1_1S","1_3MI","1_3S","8_1MI","8_1S","8_3MI","8_3S"),
                 project = "All_MFs")

#We will remove the individual seurat objects to save memory
rm(list=ls(pattern="rna_"))

#Remove apoptotic cells with a lot of mitochondrial genes
All_rna[["percent.mt"]] <- PercentageFeatureSet(object = All_rna, pattern = "^mt-")
head(x = All_rna@meta.data, 5)
All_rna #15147 features across 6503 cells
#Show how many genes, reads, and dead cells there are before filtering
VlnPlot(object = All_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = .05)
#Dead cells should have less reads, bigger cells (more RNA) should have genes (features)
pdf("ViolinAllRNABeforeMito.pdf")
plot1 <- FeatureScatter(object = All_rna, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(object = All_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2
dev.off()
#Get rid of cells outside of read BCs and dead cells with a lot of mito genes
All_rna2 <- subset(x = All_rna, subset = nFeature_RNA > 1000 & nFeature_RNA < 5800 & percent.mt < 18)
All_rna #19016 genes across 18431 cells
All_rna2 #19016 genes across 15739 cells
table(All_rna@meta.data$orig.ident) #1806 Cont MFs and 4697 MI MFs
table(All_rna2@meta.data$orig.ident) #1806 Cont MFs and 4697 MI MFs


#Normalize
All_rna2 <- NormalizeData(object = All_rna2, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = TRUE)
pdf("ViolinAllMFsAfterMito.pdf")
plot3 <- FeatureScatter(object = All_rna2, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot4 <- FeatureScatter(object = All_rna2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot3 + plot4 #compare to earlier plots and see that correlation with size is stronger
dev.off()
#Find Variable Features and Normalize
All_rna2 <- FindVariableFeatures(All_rna2)
length(x = All_rna2@assays$RNA@var.features) # just top 2000 diff exp
all.genesRNA <- rownames(x = All_rna2)

top10DiffGenes <- head(x = VariableFeatures(object = All_rna2), 10)
pdf("Top10DiffGenesMItoSS.pdf", width = 10, height = 6)
plot1 <- VariableFeaturePlot(object = All_rna2)
plot2 <- LabelPoints(plot = plot1, points = top10DiffGenes, repel = TRUE)
plot1 + plot2
dev.off()

#Time to Scale the Data!
All_rna2 <- SCTransform(All_rna2, vars.to.regress = "percent.mt", verbose = TRUE) #Correcting for nUMI is default YES

#Run UMAP
All_rna2 <- RunPCA(All_rna2, verbose = FALSE)
All_rna2 <- RunUMAP(All_rna2, dims = 1:30, verbose = FALSE)
# Make clusters
All_rna2 <- FindNeighbors(All_rna2, dims = 1:30, verbose = FALSE) 
All_rna2 <- FindClusters(All_rna2, verbose = FALSE)

pdf("All_Samples_RNA_Dimplot.pdf")
DimPlot(All_rna2, label = TRUE, group.by = "orig.ident") 
dev.off()

pdf("RNA_Dimplot_by_Cluster.pdf")
DimPlot(All_rna2, label = TRUE, group.by = "ident") + NoLegend()
dev.off()

#Find Markers for each cluster
All_rna2.markers <- FindAllMarkers(All_rna2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # this takes a while
#All_rna2.markers %>% group_by(cluster) %>% slice_max(order_by = "avg_logFC", n = 1) -> Top1Markers
All_rna2.markers %>% group_by(cluster) %>% slice_head(n = 2) -> Top2Markers

#export markers and top markers as .csv files
write.csv(Top2Markers, file = "top2markersbycluster.csv")
write.csv(All_rna2.markers, file = "allmarkersbycluster_small_clust.csv")

all_genes <- as.data.frame(all.genesRNA)
write.csv(all_genes, file = "allgenes.csv")

#Verify Specificity of Markers on UMAP plot
#These are the top ones 
pdf("All_MFs2_TopMarkeronUMAP.pdf")
for (i in 1:length(Top2Markers$gene)){
print(FeaturePlot(All_rna2, features = Top2Markers$gene[i], label = TRUE))}
dev.off()

# This will give you individual plots
# for (i in 1:length(Top2Markers$gene)){
#   temp_plot <- FeaturePlot(All_rna2, features = Top2Markers$gene[i], label = TRUE)
#   ggsave(temp_plot, file=paste0(Top2Markers$gene[i],"featureplot.png"), width = 14, height = 10, units = "in")
# }


# ID Cell Types by Really Granular Clustering -----------------------------

# 
# #The clusters we got were too small, so we should try again with less granular clustering
#All_rna2 <- FindNeighbors(All_rna2, dims = 1:30, verbose = FALSE)  
All_rna2 <- FindClusters(All_rna2, verbose = FALSE, resolution = 0.1)
All_rna2.markers <- FindAllMarkers(All_rna2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
All_rna2.markers %>% group_by(cluster) %>% slice_head(n = 2) -> Top2Markers
# 
# pdf("All_MFs2_TopMarkeronUMAP.pdf")
# for (i in 1:length(Top2Markers$gene)){
#   print(FeaturePlot(All_rna2, features = Top2Markers$gene[i], label = TRUE))}
# dev.off()
# 

DimPlot(All_rna2, label = TRUE)
# # This tell us what the major cell types are 
# # Cluster 0 = CF = DCN, Col1A1
# # Cluster 1 = Endothelial Cells = CD36, FABP4
# # Cluster 2 = SMCs = Acta2, Rgs5
# # Cluster 3 = Macrophage/Monocyte = Lyz2, Spp1
# # Cluster 4 = Epicardial Cells =  C3 + MSLN
# # Cluster 5 = ? = Bace2, Plvap
# # Cluster 6 = Ccr7, Cytip - resting T cells?
# # Cluster 7 = some other macrophge?
# # Cluster 8 = Myoz2 = CM!
# # Cluster 9 = Alas2, Hbb-bh1 = RBCs - according to author
# #

# Now we need to combine some of the smaller clusters / remove some


DefaultAssay(All_rna2) <- 'RNA' # Use the non-normalized values 
All_rna <- FindVariableFeatures(
   object = olson_rna,
   nfeatures = 10000)
 
saveRDS(All_rna2, file = "olson_rna_granular.RDS")
all_rna_gran <- readRDS(file = "olson_rna_granular.RDS")

all_rna_gran <- RenameIdents(all_rna_gran, `0` = "CF", `1` = "EC", `2` = "SMC",'3'= "MonoMac", 
                         `4` = "EpiC", `5` = "???", `6` = "T", `7` = "MonoMac2",
                         `8` = "CM", `9` = "RBC") 
# Find RNA clusters -------------------------------------------------------

# Here are what the clusters are:
# Cluster 0 = CF
# Cluster 1 = Endo
# Cluster 2 = Endo
# Cluster 3 = CF
# Cluster 4 = SMC
# Cluster 5 = Endo
# Cluster 6 = MonoMac
# Cluster 7 = MonoMac
# Cluster 8 = CF
# Cluster 9 = SMC
# Cluster 10 = CF
# Cluster 11 = Endo
# Cluster 12 = CF
# Cluster 13 = Endo 
# Cluster 14 = CF 
# Cluster 15 = Endo  (VWF, Plvap)
# Cluster 16 = Epicardial
# Cluster 17 = MonoMac 
# Cluster 18 = CFs (Tnc, Timp1)
# Cluster 19 = Cardiomyocytes
# Cluster 20 = B cells
# Cluster 21 = Pericyte/ SMC
# Cluster 22 = T cells
# Cluster 23 = RBCs

#Let's remove the RBC cluster
DimPlot(All_rna2)
All_rna2 <- subset(All_rna2, idents = "23", invert = TRUE)
#Let's rename the clusters to that they are by cell type
table(Idents(All_rna2))

All_rna2 <- RenameIdents(All_rna2, `0` = "CF1", `1` = "Endo1", `2` = "Endo2",'3'= "CF2", 
                        `4` = "SMC1", `5` = "Endo3", `6` = "MonoMac1", `7` = "MonoMac2",
                        `8` = "CF3", `9` = "SMC2", '10'="CF4",'11'="Endo4",'12'="CF5",
                        '13' = "Endo5",'14'="CF6",'15'="Endo6",'16'="Epicardial",'17'="MonoMac3",
                        '18'="CF7",'19'="Cardiomyocyte",'20'="B Cells",'21'="SMC3/Peri",'22'="T Cells")
#Now we will reorder the clusters for easy graphing
mylevels <- c("CF1","CF2","CF3","CF4","CF5","CF6","CF7","SMC1","SMC2","SMC3/Peri","Endo1",
                      "Endo2","Endo3","Endo4","Endo5","Endo6","MonoMac1","MonoMac2","MonoMac3",
                      "Cardiomyocyte","B Cells","T Cells", "Epicardial")
All_rna2@active.ident <- factor (x = All_rna2@active.ident, levels = mylevels )
All_rna2$SubCellTypes <- Idents(All_rna2)
table(Idents(All_rna2))

# We will copy orig.ident as TreatmentGroup 
table(All_rna2$orig.ident)
Idents(All_rna2) <- "orig.ident"
All_rna2$TreatmentGroup <- Idents(All_rna2)
DimPlot(All_rna2)


# Save the Seurat Object instead of all global objects to save space
saveRDS(All_rna2, file = "scRNA_data.RDS")
All_rna2 <- readRDS("scRNA_data.RDS")


#Let's create a heatmap of markers
Idents(All_rna2) <- "SubCellTypes"
cluster.averages <- AverageExpression(All_rna2, return.seurat = T)
dimension1Features <- unlist(TopFeatures(All_rna2[["pca"]],balanced = TRUE, dim =1))
dimension2Features <- unlist(TopFeatures(All_rna2[["pca"]],balanced = TRUE, dim =2))
dimension3Features <- unlist(TopFeatures(All_rna2[["pca"]],balanced = TRUE, dim =3))
dimensionalFeatures <- c(dimension1Features,dimension2Features,dimension3Features)

pdf("cluster_heatmaps.pdf")
DoHeatmap(cluster.averages, features = c(dimensionalFeatures[1:40],"Tnnt2","Cd79a","Trbc1","Msln"),
          size=3, draw.lines = F)
DoHeatmap(cluster.averages, features = unlist(Top2Markers[1:46,7]),
          size=3, draw.lines = F)
dev.off()

FeaturePlot(All_rna2, features = c("En1","Ogn"),  label = T)

Idents(All_rna2)<- "predicted.id"
All_rna2 <- RenameIdents(All_rna2,'MonoMac3' = "Cardiomyocyte2")
FeaturePlot(All_rna2, features = c("Ccl2"),  label = T)


# Compare RNA-seq Values with Dick SC Data -----------------------------------------------------------------
#Let's create a subsetted seurat object that only contains mono/mac
mac_rna <- subset(All_rna2, idents = c("MonoMac1","MonoMac2","MonoMac3"))
DimPlot(mac_rna)

#Save seurat object
saveRDS(mac_rna, file = "mac_rna.rds")
mac_rna <- readRDS("mac_rna.rds")

# Try to integrate with Harmony instead

#Let's read in the sc RNA seq from Dick et al and see how the cells cluster in aggregate.
dick <- readRDS(file='/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/Dick_Human_SSMac/Alex_Analysis_Dick_SS/All_MFs2.rds')
dick$author <- "dick"
mac_rna$author <- "olson"
agg.list <- list(dick, mac_rna)
agg.list <- lapply(X = agg.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = agg.list)

mac.anchors <- FindIntegrationAnchors(object.list = agg.list,anchor.features = features)
agg_mac_rna <- IntegrateData(anchorset = mac.anchors)

DefaultAssay(agg_mac_rna) <- "integrated"
agg_mac_rna <- ScaleData(agg_mac_rna, verbose = FALSE)
agg_mac_rna <- RunPCA(agg_mac_rna, npcs = 30, verbose = FALSE)
agg_mac_rna <- RunUMAP(agg_mac_rna, dims = 1:30, verbose = FALSE, reduction = "pca")
agg_mac_rna <- FindNeighbors(agg_mac_rna, dims = 1:30, verbose = FALSE, reduction = "pca") 
agg_mac_rna <- FindClusters(agg_mac_rna, verbose = FALSE, resolution = 0.2)

DimPlot(agg_mac_rna)
p1 <- DimPlot(agg_mac_rna, reduction = "umap", group.by = "author")
p2 <- DimPlot(agg_mac_rna, reduction = "umap",label = TRUE, repel = TRUE)
p1 + p2


FeaturePlot(agg_mac_rna, feature = c("Ccr2","Cxcr4"), label = T, blend = T)
FeaturePlot(agg_mac_rna, feature = c("Timd4"), label = T)
agg_mac_rna.markers <- FindAllMarkers(agg_mac_rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # this takes a while
agg_mac_rna.markers %>% group_by(cluster) %>% slice_head(n = 2) -> Top2MacAggMarkers

# Cluster 0 - Antigen Presnting Macs (MHC2 hi)
# Cluster 1 - Cathepsin D and Gpnmb Hi Macs (Likely CXCR4 hi)
# Cluster 2 - YS Macs
# Cluster 3 - Mainly in dick, Ly6C hi 
# Cluster 4 - Dividing cells - survivin, stathmin, ki67

# We will remove the contaminating cells
# Cluster 5 - Only in Dick - Activated? Degranulating?  Mmp9 + CD244 
# Cluster 6- Endothelial contamination? Mainly in Dick  - iNOS creating,
# Cluster 7 - B Cells -  CD79a, Ly6d- Ms4A1
# Cluster 8  - Fibrocyte?

DefaultAssay(agg_mac_rna) <- "RNA"
#test <- FindConservedMarkers(agg_mac_rna, grouping.var = "author", verbose = FALSE)
#head(nk.markers)

agg_mac_rna <- RenameIdents(agg_mac_rna, `0` = "APC Macs", `1` = "CXCR4", `2` = "YS Macs", 
                                `3` = "Ly6C Hi Macs", `4` = "Dividing", `5` = "Degranulating", `6` = "Endothelial?", `7` = "B Cells", 
                            `8` = "Fibrocyte?")
DimPlot(agg_mac_rna, label = TRUE)


table(Idents(agg_mac_rna))
agg_mac_rna2 <- subset(agg_mac_rna, idents = c("APC Macs","CXCR4","YS Macs","Ly6C Hi Macs"))

DimPlot(agg_mac_rna2)
FeaturePlot(agg_mac_rna2, features = c("Spp1", "Cxcr4","Ccl24"), label = T)
VlnPlot(agg_mac_rna2, features =c("Spp1", "Cxcr4","Ccl24","Ccr2"))

DefaultAssay(agg_mac_rna2) <- "integrated"
agg_mac_rna2 <- ScaleData(agg_mac_rna2, verbose = FALSE)
agg_mac_rna2 <- FindVariableFeatures(agg_mac_rna2, selection.method = "vst", nfeatures = 2000)
agg_mac_rna2 <- RunPCA(agg_mac_rna2, npcs = 30, verbose = FALSE)
agg_mac_rna2 <- RunUMAP(agg_mac_rna2, dims = 1:30, verbose = FALSE, reduction = "pca")
agg_mac_rna2 <- FindNeighbors(agg_mac_rna2, dims = 1:30, verbose = FALSE, reduction = "pca") 
agg_mac_rna2 <- FindClusters(agg_mac_rna2, verbose = FALSE, resolution = 0.15)

DefaultAssay(agg_mac_rna2) <- "RNA"
DimPlot(agg_mac_rna2, label = T, split.by = "author")

agg_mac_rna.markers <- FindAllMarkers(agg_mac_rna2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # this takes a while
agg_mac_rna.markers %>% group_by(cluster) %>% slice_head(n = 2) -> Top2MacAggMarkers




#DefaultAssay(agg_mac_rna) <- "RNA"
#DimPlot(agg_mac_rna, reduction = "umap", split.by = 'author', label = T)
#VlnPlot(agg_mac_rna, features = c("Timd4","Ccr2","Cxcr4"), pt.size = 0.01)
FeaturePlot(agg_mac_rna2, features = c("Spp1", "Cxcr4","Timd4","Ccr2"), label = T)
DefaultAssay(agg_mac_rna2) <- "SCT"
FeatureScatter(agg_mac_rna2, "Ccr2", "Cxcr4", group.by = "author", slot = "counts")
DefaultAssay(agg_mac_rna2) <- "integrated"
print(agg_mac_rna2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(agg_mac_rna2, dims = 1:2, reduction = "pca")
DimHeatmap(agg_mac_rna2, dims = 1, cells = 500, balanced = TRUE)
DefaultAssay(agg_mac_rna2) <- "SCT"
FeaturePlot(agg_mac_rna2, features = c("Spp1","Sepp1"), label = T)
FeatureScatter(agg_mac_rna2, "S100a11", "Igfbp4", group.by = "author")
FeatureScatter(agg_mac_rna2, "S100a11", "Fxyd5", group.by = "author")
FeaturePlot(agg_mac_rna2, features = c("S100a11","Igfbp4"), label = T)
FeatureScatter(agg_mac_rna2, "Cxcr4", "Ccr2", group.by = "author")
test <- subset(agg_mac_rna2, idents = c("0","2"))
FeatureScatter(test, "Cxcr4", "Ccl2", group.by = "author", slot = "counts")

#~~~~~~~~~ Create New Cluters Based on CCR2, Timd4, and CXCR4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Save existing ident as cluster.ident
# MI_MFs2$cluster_names = Idents(MI_MFs2)
# 
# chemocluster = FetchData(MI_MFs2,c("ident","orig.ident","Ccr2","Cxcr4","Timd4"))
# Ccr2Cells = WhichCells(object = MI_MFs2, expression = Ccr2 > 1)
# FeaturePlot(MI_MFs2, cells = Ccr2Cells, features = c("Ccr2"), label = TRUE)
# Cxcr4Cells = WhichCells(object = MI_MFs2, expression = Cxcr4 >1)
# length(table(Ccr2Cells)) #1783
# FeaturePlot(MI_MFs2, cells = Cxcr4Cells, features = c("Cxcr4"), label = TRUE)
# length(table(Cxcr4Cells)) #728
# Timd4Cells = WhichCells(object = MI_MFs2, expression = Timd4 >0.5)
# FeaturePlot(MI_MFs2, cells = Timd4Cells, features = c("Timd4"), label = TRUE)
# length(table(Timd4Cells)) #348
# MI_MFs2 <- SetIdent(object = MI_MFs2, value = "Other")
# MI_MFs2 <- SetIdent(object = MI_MFs2, cells = Ccr2Cells, value = "Ccr2Cells")
# MI_MFs2 <- SetIdent(object = MI_MFs2, cells = Cxcr4Cells, value = "Cxcr4Cells")
# MI_MFs2 <- SetIdent(object = MI_MFs2, cells = Timd4Cells, value = "Timd4Cells")
# pdf("Chemo_on_MI.pdf")
# FeaturePlot(MI_MFs2, features = c("Ccr2","Cxcr4","Timd4"), label = TRUE)
# dev.off()

# #Create new seurat with 60% of cells that bin niceley with expression measured
# ChemoSeurat <- subset(MI_MFs2, idents = c("Ccr2Cells","Cxcr4Cells","Timd4Cells"))
# # Reset MI_MFs2 Idents to clusters
# Idents(MI_MFs2) <- MI_MFs2$cluster_names
# #recluster and bin cells from new ChemoSeurat
# ChemoSeurat <- FindVariableFeatures(ChemoSeurat)
# 
# 
# #ChemoSeurat <- JackStraw(object = ChemoSeurat, num.replicate = 100) 
#ChemoSeurat <- ScoreJackStraw(ChemoSeurat, dims = 1:20)
#pdf("Chemo MFs Dimensions.pdf")
#ElbowPlot(object = ChemoSeurat)
#JackStrawPlot(ChemoSeurat, dims = 1:20)
#dev.off()
# 
# ChemoSeurat <- RunUMAP(ChemoSeurat, dims = 1:30, verbose = FALSE)
# 
# # With the new Seurat update, you don't need to SCTTransform again, BUT you can't do Jackstraw is subsetting
# #ChemoSeurat <- SCTransform(ChemoSeurat, vars.to.regress = "percent.mt", verbose = TRUE)
# #ChemoSeurat <- RunPCA(ChemoSeurat, verbose = FALSE)
# 
# ChemoSeurat <- FindNeighbors(ChemoSeurat, dims = 1:20, verbose = TRUE)
# ChemoSeurat <- FindClusters(ChemoSeurat, resolution = 0.8, verbose = TRUE)
# 
# pdf("ID_Chemo_Clusters.pdf")
# FeaturePlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4"), label = TRUE)
# VlnPlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4"))
# dev.off()
# 
# #Merge Clusters to obtain 3 main clusters: CCR2, CXCR4, TIMD4
# ChemoSeurat <- RenameIdents(object = ChemoSeurat, "4" = "Timd4","5"="Ccr2","2"="Ccr2",
#                             "9"="Ccr2","6"="Cxcr4","7"="Ccr2","8"="Ccr2","1"="Ccr2",
#                             "0"="Ccr2","3"="Cxcr4")
# #Groups 3 and 5 seemed to be double positive or ambiguous
# pdf("Chemo_New_IDs.pdf")
# FeaturePlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4"), label = TRUE)
# VlnPlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4"))
# DimPlot(ChemoSeurat, label = TRUE, group.by = "ident") + NoLegend()
# DimPlot(ChemoSeurat, label = TRUE, group.by = "orig.ident") + NoLegend()
# dev.off()
# 
# #Find DEGS Between the three clusters
# ChemoSeurat.markers <- FindAllMarkers(ChemoSeurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# ChemoSeurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) -> ChemoTop2Markers
# FeaturePlot(ChemoSeurat, features = c("Gpnmb","Timd4","Tgfbr1"), label = TRUE)
# 
# #Get CXCR4 Markers and then Run through GO Term Screen
# CXCR4_Markers <- ChemoSeurat.markers[c(204:328),(7)]
# write.csv(CXCR4_Markers, file = "CXCR4 Markers for GO.csv")
# 


# scATAC  -----------------------------------------------------------------
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)  
library(EnsDb.Mmusculus.v79)# here we will use EnsDb.Mmusculus.v79
library(ggplot2)
library(patchwork)
set.seed(1234)
library(hdf5r)

atac_counts <-Read10X_h5(filename = "/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scATAC/GSE153479_filtered_peak_bc_matrix.h5")


atac_metadata <- read.csv(file = "/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scATAC/GSE153479_singlecell.csv",
                          header = T, row.names =1)

chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(":","-"), genome = 'mm10',
                                    fragments = '/Volumes/TUNEZ/TSCC_Files/scratch/olson_sc_rna_atac/scATAC/GSE153479_fragments.tsv.gz',
                                    min.cell = 10, min.features = 200)

olson_sc_atac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = atac_metadata)

#we will remove chrom_assay and atac_counts to save memory
rm(chrom_assay,atac_counts)

olson_sc_atac[['peaks']]
granges(olson_sc_atac)

# Annotate Regions with Ensembl Database using UCSC names
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
Annotation(olson_sc_atac) <- annotations

# QC
#Compute Nucleosome Signal score per cell
olson_sc_atac <- NucleosomeSignal(object = olson_sc_atac)

# Compute TSS enrichment score per cell
olson_sc_atac <- TSSEnrichment(olson_sc_atac, fast = FALSE)

# add blacklist ration and fraction of reads in peaks
olson_sc_atac$pct_reads_in_peaks <- olson_sc_atac$peak_region_fragments / olson_sc_atac$passed_filters * 100
olson_sc_atac$blacklist_ratio <- olson_sc_atac$blacklist_region_fragments / olson_sc_atac$peak_region_fragments

#Visualize High and Low TSS Enrichment in dataset
olson_sc_atac$high.tss <- ifelse(olson_sc_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(olson_sc_atac, group.by = 'high.tss') + NoLegend()

# Visualize nucleosome enrichment of groups of cells
olson_sc_atac$nucleosome_group <- ifelse(olson_sc_atac$nucleosome_signal > 0.8, 'NS > 0.8', 'NS < 0.8')
FragmentHistogram(object = olson_sc_atac, group.by = 'nucleosome_group')

# Look at distributions of: 
pdf("sc_atac_QC.pdf", width = 15, height = 15)
VlnPlot(
  object = olson_sc_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.01,
  ncol = 5
)
dev.off()

# Remove outliers for these QC criteria - these parameters were used in the publication
olson_sc_atac <- subset(
  x = olson_sc_atac,
  subset = peak_region_fragments > 5000 &
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
olson_sc_atac

# We need to add the sample IDs from the barcode suffixes
#Name	    Sample	  Sample Index	Cell Barcode Suffix
#P1D3MI  	P1+3 dpi	SI-NA-G1           	-1
#P1D3Sham	P1+3 dps	SI-NA-H1          	-2
#P8D3MI  	P8+3 dpi	SI-NA-A2           	-3
#P8D3Sham	P8+3 dps	SI-NA-B2	          -4
olson_sc_atac@meta.data$condition <- rownames(olson_sc_atac@meta.data)
library(stringr)
olson_sc_atac@meta.data$condition <- substring(olson_sc_atac@meta.data$condition, 18,18)
olson_sc_atac@meta.data$condition <- factor(olson_sc_atac@meta.data$condition, 
                                            levels = c("1","2","3","4"), 
                                            labels = c("atac_1_3MI","atac_1_3S","atac_8_3MI","atac_8_3S"))
olson_sc_atac@meta.data$condition
Idents(olson_sc_atac) <- olson_sc_atac@meta.data$condition
olson_sc_atac$sample <- Idents(olson_sc_atac)

#normalize via TF-IDF and do single variable decomposition (together this is called LSI)
olson_sc_atac <- RunTFIDF(olson_sc_atac)
olson_sc_atac <- FindTopFeatures(olson_sc_atac, min.cutoff = 'q0')
olson_sc_atac <- RunSVD(olson_sc_atac)

# We need to see if there is a strong correlation between sequencing depth and biological variance
# Generally, the first principal component in sc-ATAC is sequencing variablility, not biological differences
DepthCor(olson_sc_atac)
# Yes, the first principal component is highly correlated with sequencing depth

# Perform Dimension Reduction
olson_sc_atac <- RunUMAP(object = olson_sc_atac, reduction = 'lsi', dims = 2:30)
olson_sc_atac <- FindNeighbors(object = olson_sc_atac, reduction = 'lsi', dims = 2:30)
olson_sc_atac <- FindClusters(object = olson_sc_atac, verbose = FALSE, algorithm = 3)
DimPlot(object = olson_sc_atac, label = TRUE) + NoLegend()

# assign genes to regions based on intersecting 2kb upstream of gene start sites
# You can also do this using Cicero (check out later)
gene.activities <- GeneActivity(olson_sc_atac)
olson_sc_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
# Add Gene Activity Matrix to Seurat ATAC matrix and normalize based on sequencing depth
olson_sc_atac <- NormalizeData(
  object = olson_sc_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(olson_sc_atac$nCount_RNA)
)

# Change the default assay of the ATAC seq object to display promoter intersections
olson_sc_atac$author <- 'Olson'
DefaultAssay(olson_sc_atac) <- 'RNA'
DefaultAssay(olson_sc_atac) <- 'peaks'

# Here we will save objects so that they can be loaded later
saveRDS(olson_sc_atac, file = "sc_atac.RDS")
readRDS(file = "sc_atac.RDS")


# Idk what features to use here, check back after integration with scRNA
FeaturePlot(
  object = olson_sc_atac,
  features = c('Col1a2', 'Cd34'),
  pt.size = 0.1,
  max.cutoff = 'q95')

DefaultAssay(olson_sc_atac) <- 'peaks'
sc_atac.markers <- FindAllMarkers(olson_sc_atac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # this takes a while
#All_rna2.markers %>% group_by(cluster) %>% slice_max(order_by = "avg_logFC", n = 1) -> Top1Markers
sc_atac.markers %>% group_by(cluster) %>% slice_head(n = 2) -> Top2atacMarkers
Top2atacMarkers_genes <- ClosestFeature(olson_sc_atac, regions = Top2atacMarkers$gene)




pdf("ATAC_features.pdf", width = 20, height = 60)
FeaturePlot(
  object = olson_sc_atac,
  features = Top2atacMarkers_genes$query_region,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  label = T)
dev.off()

pdf("ATAC_features_Publication.pdf")
FeaturePlot(
  object = olson_sc_atac,
  features = c("chr11-94940531-94951370","chr5-122092239-122095914",
               "chr2-26582216-26598835","chr18-61657229-61661148",
               "chr18-61107669-61109014","chr11-103174811-103178549",
               "chr2-33129447-33132423"), 
  # Col1A1(CF), Myl2(CM), Agpat2(EC),IL17b(SMC), Cst1r(Mac), Fmnl1 (Lymphocyte), Garnl3 (Epi)
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  label = T)
VlnPlot(olson_sc_atac, features = c("chr11-94940531-94951370","chr5-122092239-122095914",
                                    "chr2-26582216-26598835","chr18-61657229-61661148",
                                    "chr18-61107669-61109014","chr11-103174811-103178549",
                                    "chr2-33129447-33132423"))
dev.off()


pdf("test.pdf", width = 15, height = 40)
FeaturePlot(
  object = olson_rna,
  features = c("Myh7","Fgfr2","Pnlip","Phldb2","Rapgef2","Ano2","Gm8730","A930011G23Rik","Pid1","Mir142hg","Traf3ip3"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  label = T)
dev.off()

Idents(olson_sc_atac) <- 'seurat_clusters'
pdf("test.pdf", width = 15, height = 40)
FeaturePlot(
  object = olson_rna,
  features = c("Nkx2-5", "Ets1","Tcf21","Ebf1"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  label = T)
dev.off()

DefaultAssay(olson_sc_atac) <- 'rna'
Idents(olson_sc_atac) <- 'sample'
Idents(olson_sc_atac) <- 'seurat_clusters'

pdf("coverage.pdf")
CoveragePlot(
  object = olson_sc_atac,
  region = "Ets1",
  extend.upstream = 40000,
  extend.downstream = 20000
)
dev.off()
# We need to define ATAC clusters

DimPlot(olson_sc_atac, label = T)

atac_marker_anno <- ClosestFeature(olson_sc_atac, regions = sc_atac.markers$gene)
colnames(sc_atac.markers)[7] <- "query_region"
atac_markers <- inner_join(sc_atac.markers, atac_marker_anno, by = "query_region")

# Cluster 0 - CM 1 
# Cluster 1 - EC 1 
# Cluster 2 - CM 2 
# Cluster 3 - CF 1
# Cluster 4 - CF 2
# Cluster 5 - CM 3
# Cluster 6 - SMC 
# Cluster 7 - CF 3
# Cluster 8 - IDK
# Cluster 9 - MonoMac/Lymphocyte 1 
# Cluster 10 - CM 4
# Cluster 11 - EC 2
# Cluster 12 - CF 4
# Cluster 13 - IDK
# Cluster 14 - MonoMac/Lymphocyte 2
# Cluster 15 - IDK

atac_markers$cluster <- as.factor(atac_markers$cluster)
#test <- dplyr::filter(atac_markers,  p_val_adj < 0.1, cluster ==15)
#print(head(test$gene_id, 100))

Idents(olson_sc_atac) <- 'seurat_clusters'
UMAPPlot(olson_sc_atac, label = T)

DefaultAssay(olson_sc_atac) <- 'RNA'
FeaturePlot(olson_sc_atac, features = c("Tgfb1"))

#Save the SC_atac object for now, will update later
saveRDS(olson_sc_atac, file = "scATAC_data_nolab.RDS")

# Integrate scATAC and scRNA ----------------------------------------------



# Make sure scRNA's active assay is RNA (and not SCT)
All_rna2 <- readRDS("scRNA_data.RDS")
olson_rna <- All_rna2
#olson_sc_atac <- readRDS(file = "sc_atac_nolab.RDS")
#olson_rna_gran <- readRDS(file = "olson_rna_granular.RDS")


# This code chunk integrates granular clusters
# DefaultAssay(olson_rna_gran) <- 'RNA' # Use the non-normalized values 
# olson_rna <- FindVariableFeatures(
#   object = olson_rna_gran,
#   nfeatures = 10000)
# 
# DefaultAssay(olson_sc_atac) <- 'RNA'
# 
# transfer.anchors <- FindTransferAnchors(
#   reference = olson_rna_gran,
#   query = olson_sc_atac,
#   reduction = 'cca',
#   dims = 1:40)
# 
# predicted.labels <- TransferData(
#   anchorset = transfer.anchors,
#   refdata = olson_rna_gran$seurat_clusters,
#   weight.reduction = olson_sc_atac[['lsi']],
#   dims = 2:30)
# 
# olson_sc_atac <- AddMetaData(object = olson_sc_atac, metadata = predicted.labels)

# This code chunk integrates based on Smaller Clustering
DefaultAssay(olson_rna) <- 'RNA' # Use the non-normalized values 
Idents(olson_rna) <- 'SubCellTypes'
olson_rna <- FindVariableFeatures(
  object = olson_rna,
  nfeatures = 10000)

DefaultAssay(olson_sc_atac) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = olson_rna,
  query = olson_sc_atac,
  reduction = 'cca',
  dims = 1:40)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = olson_rna$SubCellTypes,
  weight.reduction = olson_sc_atac[['lsi']],
  dims = 2:30)

olson_sc_atac <- AddMetaData(object = olson_sc_atac, metadata = predicted.labels)

plot1 <- DimPlot(
  object = olson_rna,
  group.by = 'SubCellTypes',
  label = TRUE,
  repel = TRUE) + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = olson_sc_atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + ggtitle('scATAC-seq')

plot1 + plot2
ggsave("test.png", width = 15, height = 8)


saveRDS(olson_sc_atac, file = "scRNA_data.RDS")

plot1 <- DimPlot(object = olson_sc_atac,
    group.by = 'rna_clusters',
    label = TRUE,
    repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot2 <- DimPlot(object = olson_rna,
              group.by = 'SubCellTypes',
              label = TRUE,
              repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot1 + plot2
ggsave("Dimplots_rna_and_atac.png")

Idents(olson_sc_atac) <- "predicted.id"
DefaultAssay(olson_sc_atac) <- 'RNA'
#smc1 <- WhichCells(olson_sc_atac,idents = "SMC1" )
#smc2 <- WhichCells(olson_sc_atac,idents = "SMC2" )
#monomac3 <- WhichCells(olson_sc_atac,idents = "MonoMac3" )
#DimPlot(olson_sc_atac, cells.highlight= list(monomac3, smc1), cols.highlight = c("darkblue", "darkred"), cols= "grey")
DimPlot(olson_sc_atac, label = T)


# The predicted cluster IDs do not match perfectly, so we will have to manually correct some:

# Most CFs can be ID'ed by Col1a1
# CF1 = CF (FGF18)
# CF2 = CF (Fbln, Loxl)
# CF3 = CF (matrix)
# CF4 = CF (matrix)
# CF5 = CF (matrix)
# CF6 = CF (Fbln, Loxl)
# CF7 = CF (matrix)

# Most Endo's can be traced by Dll4 
# Endo1 = Endo 
# Endo2 = Endo 
# Endo3 = Endo 
# Endo4 = Endo 
# Endo5 = Endo 
# Endo6 = Endo 

# Most SMCs can be ID'ed by ???
# SMC1 = Unsure
# SMC2 = Classic Mature SMCs (Myh11, Rgs5)

# Cardiomyocyte = Cardiomyocyte (myosin heavy chains)
# B Cells = B cells (Cd3d)
# MonoMac1 = YS Containing Macs (Ccl24)
# MonoMac2 = Spp1 Making Macs, Ccr2 making macs? Mixed bag
# MonoMac3 = Cardiomyocytes (Lots of CM genes)
# Epicardial = Epicardial (Wt1, C3)


FeaturePlot(olson_sc_atac, features = c("Myh11","Cd3d"), label = T)


VlnPlot(olson_sc_atac, features = c("Spp1","Ccl24","Ccr2","Cxcr4"))

olson_sc_atac.rnamarkers <- FindAllMarkers(olson_sc_atac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # this takes a while
write.csv(olson_sc_atac.rnamarkers, "olson_sc_atac_rna_markers.csv")
pdf("ATAC_Markers_Cluster_ID.pdf")
FeaturePlot(olson_sc_atac, features = c("Dll4","Col1a1", "Socs2"),pt.size = 0.1, label = T)
dev.off()
FeaturePlot(olson_rna, features = c("Cxcr4","Ccr2"), label = T)

# After we fix the cluster IDs, then we need to save the object again
olson_sc_atac <- RenameIdents(olson_sc_atac,  'CF1' = "CF_a", 'CF2' = "CF_b",
                              'CF3' = "CF_c", 'CF4' = "CF_d", 'CF5' = "CF_e",
                              'CF6' = "CF_f",'CF7' = "CF_g", 'Endo1' = "Endo_a",
                              'Endo2' = "Endo_b",'Endo3' = "Endo_c",'Endo4' = "Endo_d",
                              'Endo5' = "Endo_e", 'Endo6' = "Endo_f",'SMC1' = "SMC_a",
                              'SMC2' = "SMC_b", 'B Cells' = "B Cells", 'MonoMac1' = "MonoMac_a",
                              "MonoMac2" = "MonoMac_b",'Cardiomyocyte' = "CM_a", 'MonoMac3' = "CM_b",
                              'Epicardial' = "EpiC")
pdf("DimPlot_RenamedATACClust.pdf")
DimPlot(olson_sc_atac, label = T)
dev.off()

olson_sc_atac.rnamarkers <- FindAllMarkers(olson_sc_atac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # this takes a while
write.csv(olson_sc_atac.rnamarkers, "olson_sc_atac_rna_markers.csv")
saveRDS(olson_sc_atac, file = "scRNA_data_renamed.RDS")

FeaturePlot(olson_sc_atac, features = c("Has2"), split.by = "sample", label = T)
plot1 <- VlnPlot(olson_sc_atac, features = c("Has2"), split.by = "sample")

plot1
ggsave("has2.png")
#olson_atac_mac <- subset(olson_sc_atac, idents = c("MonoMac1","MonoMac2"))
#olson_atac_mac.markers <- FindAllMarkers(olson_atac_mac, only.pos = TRUE, min.pct = 0.8, logfc.threshold = 0.5) # this takes a while
#FeaturePlot(olson_atac_mac, features = c("Cxcr4","Ccr2"), label = T)

DimPlot(olson_sc_atac)
olson_sc_atac$rna_clusters <- Idents(olson_sc_atac)
DefaultAssay(olson_sc_atac) <- 'peaks'
plot1 <- FeaturePlot(olson_sc_atac, features = "Col1a2")


saveRDS(olson_sc_atac, file = "scRNA_data_renamed.RDS")

plot1 <- CoveragePlot(object = olson_sc_atac, region = c("Col1a2"))
ggsave("Col1A1_Cov_Plot.png")


olson_sc_atac <- readRDS("scRNA_data_renamed.RDS")
FeaturePlot(olson_sc_atac, features = c("En1"), label = T)


# CF only analysis --------------------------------------------------------



CFs_atac <- subset(olson_sc_atac, idents =c("CF_a","CF_b","CF_c","CF_d","CF_e","CF_f","CF_g"))
DimPlot(CFs_atac)
CFs_atac <- RunUMAP(object = CFs_atac, reduction = 'lsi', dims = 2:30)
CFs_atac <- FindNeighbors(object = CFs_atac, reduction = 'lsi', dims = 2:30)
CFs_atac <- FindClusters(object = CFs_atac, verbose = FALSE, algorithm = 3)
DimPlot(object = CFs_atac, label = TRUE) 
DefaultAssay(CFs_atac) <- 'RNA'
CF_atac_rna_markers <- FindAllMarkers(CFs_atac)
DefaultAssay(CFs_atac) <- 'peaks'
CF_atac_peak_markers <- FindAllMarkers(CFs_atac)

DimPlot(object = CFs_atac, label = TRUE, split.by = "sample") 
FeaturePlot(CFs_atac, features = c("Col1a2"), split.by = "sample")
VlnPlot(CFs_atac, features = c("Col1a2","Col1a1","Tnc","Cd44"))
FeaturePlot(CFs_atac, features = c("Col1a2","Col1a1","Tnc","Cd44"), label = T)

#Based on the new clustering of CFs, most clusters are either unique to P1 or P8
# Good: 1,4, 5?
# Bad: 2,3, 0?

# Should do motif finding for clusters 2 + 3 


DefaultAssay(CFs_atac) <- 'peaks'


pfm <- getMatrixSet(x = JASPAR2020,opts = list(species = 9606, all_versions = FALSE))

CFs_atac <- AddMotifs(
  object = CFs_atac,
  genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

da_peaks <- FindMarkers(
  object = CFs_atac,
  ident.1 = '3',
  ident.2 = '4',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

enriched.motifs <- FindMotifs(
  object = CFs_atac,
  features = top.da.peak
)

MotifPlot(
  object = CFs_atac,
  motifs = head(rownames(enriched.motifs))
)

CFs_atac <- RunChromVAR(
  object = CFs_atac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(CFs_atac) <- 'chromvar'

# look at the activity of Mef2c
plot1 <- DimPlot(CFs_atac, split.by = "sample")
plot2 <- FeaturePlot(
  object = CFs_atac,
  features = "MA1144.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  split.by = "sample",
  label = T
)
plot1 + plot2
plot1
plot2
ggsave("AP-1_Motifs_CF_sc_atac.png")

saveRDS(olson_sc_atac, file = "scRNA_data_renamed.RDS")


# Since AP-1 is the motif unique to p8 MI CFs, let's find a similar cluster in the RNA-seq
# And see what AP-1 proteins are upregulated between those clusters

DimPlot(olson_rna, label = T)
CFs_rna <- subset(olson_rna, idents = c("CF1","CF2","CF3","CF4","CF5","CF6","CF7"))
# Scale then run PCA/UMAP
CFs_rna <- NormalizeData(object = CFs_rna, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = TRUE)
CFs_rna <- FindVariableFeatures(CFs_rna)
CFs_rna <- SCTransform(CFs_rna, vars.to.regress = "percent.mt", verbose = TRUE) #Correcting for nUMI is default YES
CFs_rna <- RunPCA(CFs_rna, verbose = FALSE)
CFs_rna <- RunUMAP(CFs_rna, dims = 1:30, verbose = FALSE)
CFs_rna <- FindNeighbors(CFs_rna, dims = 1:30, verbose = FALSE) 
CFs_rna <- FindClusters(CFs_rna, verbose = FALSE, resolution = 0.1)
DimPlot(CFs_rna, label = T)
CF_rna_markers <-  FindAllMarkers(CFs_rna, only.pos = F, min.pct = 0.5, logfc.threshold = 0.25)

DimPlot(object = CFs_rna, label = TRUE, split.by = "TreatmentGroup") 
ggsave("CF_Dimplot_Split.png")

# Seems like most CFs are lost with age
# Good - 0,1,3,4
# Bad - 2 (emerges with MI only in p8)

# Check AP-1 family protein expression - see if some are limited to Cluster 2
#FeaturePlot(CFs_rna, features = c("Jun","Fos","Atf3","Fosl2","Thbs1"), label = T)
# Check subfamilies of AP-1
# Check Fos subfamily
VlnPlot(CFs_rna, features = c("Fos","Fosb",
                              "Fosl1","Fosl2"))  # The Fosls are Fra 1 + 2
# Fosl2 looks promising?

# Check ATF family
VlnPlot(CFs_rna, features = c("Atf1","Atf2","Atf3","Atf4",
                              "Atf5","Atf6","Atf6b",
                              "Atf7","Batf","Batf2",
                              "Batf3","Jdp2","Atfx"))
# Check Jun Family
VlnPlot(CFs_rna, features = c("Jun","Junb","Jund"))

# Check MAF Family
VlnPlot(CFs_rna, features = c("Maf","Mafa","Mafb","Maff",
                              "Mafg","Mafk","Nrl"))

# Let's look at the predicted RNA binding of Fosl2 in the atac data
VlnPlot(CFs_atac, features = c("Fos","Fosb",
                              "Fosl1","Fosl2"), pt.size = 0.01)
# Fosl1 has a slightly higher candle, but nothing super convincing
FeaturePlot(CFs_atac, features = c("Junb"), label = T)

# Let's put the in a heatmap for only this cluster
DoHeatmap(CFs_rna, features = c("Atf1","Atf2","Atf3","Atf4",
                                "Atf5","Atf6","Atf6b",
                                "Atf7","Batf","Batf2",
                                "Batf3","Jdp2","Jun","Junb","Jund",
                                "Fos","Fosb",
                                "Fosl1","Fosl2","Maf","Mafa","Mafb","Maff",
                                "Mafg","Mafk","Nrl"))


#Let's check chromatin remodeling proteins to see if they are driving a unique AP-1 response
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4932839/
# Let's check histone acetylases
VlnPlot(CFs_rna, features = c("Kat2a","Kat2b","Kat5",
                                  "Kat6a","Kat6b","Kat8","Kat7",
                                  "Ep300","Crebbp","Hat1",
                                  "Atf2","Ncoa1","Taf1"))
DimPlot(CFs_rna, split.by = "TreatmentGroup")
# none look promising

# Let's check histone lysine methyltransferases
VlnPlot(CFs_rna, features = c("Kmt2a","Kmt2d","Kmt2c",
                              "Kmt2b","Setd1a","Set1l",
                              "Kat7","Nsd1","Whsc1",
                              "Whsc1l1","Ash1l","Setd2",
                              "Set2l"))
        # can't find other names for Set1L or Set2L
# There are a lot so here's another plot
VlnPlot(CFs_rna, features = c("Ezh1","Ezh2","Smyd1",
                              "Smyd2","Smyd3","Smyd4",
                              "Smyd5","Setmar","Setdb1",
                              "Setdb2","Ehmt1","Ehmt2",
                              "Suv39h1","Suv39h2"))
VlnPlot(CFs_rna, features = c("Pdrm1","Pdrm2","Pdrm4",
                              "Pdrm5","Pdrm6","Pdrm7",
                              "Pdrm8","Pdrm9","Pdrm9",
                              "Pdrm10","Pdrm11","Pdrm12",
                              "Pdrm13","Pdrm14","Pdrm15",
                              "Pdrm16"))
      #Not sure what ANY of these gene names should be 
VlnPlot(CFs_rna, features = c("Suv420h1","Suv420h2","Kmt2e",
                              "Setd3","Setd4","Setd5",
                              "Setd6","Setd7","Setd8",
                              "Dot1l"))


#Let's check SWI/SNF complexes aka BAF/PBAF Complexes
VlnPlot(CFs_rna, features = c("Smarca4","Arid1a","Arid1b",
                              "Smarcb1","Phf10","Dpf1",
                              "Dpf3","Dpf2",
                              "Actl6a","Actl6b","Smarce1",
                              "Smarcd1","Smarcd2","Smarcd3",
                              "Smarcc1","Smarcc2","Ss18","Pbrm1",
                              "Brd7"))

# Now let's check ISWI Family Chromatin Remodelers
VlnPlot(CFs_rna, features = c("Smarca5","Baz1a","Smarca1",
                              "Rsf1","Cecr2","Baz1a",
                              "Chrac1","Pole3",
                              "Bptf","Rbbp7","Rbbp4",
                              "Baz2a","Baz1b","Dek","Ercc6",
                              "Nmi","Sf3b1","Mybbp1a","Ddx21"))

# Now let's check CHD family chromatin remodeling complexes
VlnPlot(CFs_rna, features = c("Chd3","Chd4","Hdac1",
                              "Hdac2","Mta1","Mta2",
                              "Mta3",
                              "Actl6a","Actl6b"))
                              
FeaturePlot(CFs_rna, features = c("Hdac1",
                                  "Hdac2","Mta1","Mta2",
                                  "Mta3"), label = T)

VlnPlot(CFs_rna, features = c("Smarce1",
                              "Gatad2a","Gatad2b","Mbd2",
                              "Mbd3","Chd5","Hdac2","Mta3",
                              "Chd1","Chd2","Chd6","Chd7",
                              "Chd8","Chd9"))

# Now let's check INO80 Family chromatin remodeling genes
VlnPlot(CFs_rna, features = c("Ino80","Ruvbl1","Ruvbl2",
                              "Mcrs1","Tfpt",
                              "Yy1","Ino80c","Ino80b",
                              "Uchl5","Nfrkb","Ino80e",
                              "Ep400","Trrap","Kat5",
                              "Brd8","Epc1","Epc2",
                              "Ncoa1","Yeats4","Dmap1","Ing3",
                              "Meaf6","Morf4l1"))
      # Not sure if Ies6, Yl2 exists in humans
VlnPlot(CFs_rna, features = c("Morf4l2","Srcap","Arp6",
                              "Mcrs1",
                              "Hit1"))
    # not sure what Arp6 or Hit1 should be 


#None of the chromatin remodeling proteins fit the trend.
# However! AP-1 is able to open chromatin by itself!!
# https://rupress.org/jem/article/217/1/e20182009/132593/AP-1-activity-induced-by-co-stimulation-is

# Can we check for the transcription of AP-1 stimulators?
VlnPlot(CFs_rna, features = c("Mapk8","Mapk14","Mapk1")) # This is JNK, p38, Erk

VlnPlot(CFs_rna, features = c("Tcf21","Nppa","Hif1a"))
VlnPlot(CFs_atac, features = c("Tcf21","Nppa","Hif1a"))


# Let's look at some downstream AP-1 genes: (from supplmenetal Fig 2 of https://www.nature.com/articles/s41467-018-08236-0)
FeaturePlot(CFs_rna, features = c("Ctgf", "Egr1","Egr2","Mafb","Nfkb2","E2f8"), label = T)
FeaturePlot(CFs_rna, features = c("Thap1", "Stat6","Hic1","Nfatc1"), label = T)
DimPlot(CFs_rna, label = T, split.by = "TreatmentGroup")

#What about NF-kB downstream genes?
FeaturePlot(CFs_rna, features = c("Nr4a2","Nfkbia"), label = T)



# Are the other matricellular proteins made by other cell types in response to MI?
FeaturePlot(olson_rna, features = c("Thbs1","Thbs2","Ogn","Sparc","Spp1", "Dcn","Lum"), label = T)
# Spp1 is mainly from macrophages
# Thbs1 + Thbs2 + Ogn + Dcn + Lum = mainly CFs
# SPARC is from all non-hematopoetic cells

FeaturePlot(olson_rna, features = c("Sparc"), label = T, split.by = "TreatmentGroup")
VlnPlot(olson_rna, features = c("Sparc", "Dcn","Lum"))
# Sparc seems to be expressed across the board


# What about ACE?
FeaturePlot(olson_rna, features = c("Ace"), label = T)
FeaturePlot(CFs_rna, features = c("Ace"), label = T)
# YES, ACE is produced primarily in non-regenerative CFs


# What about HA related signaling
FeaturePlot(olson_rna, features = c("Has1", "Has2", "Has3","Cd44","Hmmr"), label = T)
# Look closer at CF clusters
FeaturePlot(CFs_rna, features = c("Has1", "Has2", "Has3","Cd44","Hmmr"), label = T)
VlnPlot(CFs_rna, features = c("Has1", "Has2", "Has3","Cd44","Hmmr"))
# Not sure what to make of this




# Endo Only Analysis ------------------------------------------------------



## Let's look at endothelial cells that are making ACE and CXCL12
DimPlot(olson_rna, label = T)
Endo_rna <- subset(olson_rna, idents = c("Endo1","Endo2","Endo3","Endo4","Endo5","Endo6"))
# Scale then run PCA/UMAP
Endo_rna <- NormalizeData(object = Endo_rna, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = TRUE)
Endo_rna <- FindVariableFeatures(Endo_rna)
Endo_rna <- SCTransform(Endo_rna, vars.to.regress = "percent.mt", verbose = TRUE) #Correcting for nUMI is default YES
Endo_rna <- RunPCA(Endo_rna, verbose = FALSE)
Endo_rna <- RunUMAP(Endo_rna, dims = 1:30, verbose = FALSE)
Endo_rna <- FindNeighbors(Endo_rna, dims = 1:30, verbose = FALSE) 
Endo_rna <- FindClusters(Endo_rna, verbose = FALSE, resolution = 0.5)
DimPlot(Endo_rna, label = T)
Endo_rna_markers <-  FindAllMarkers(Endo_rna, only.pos = F, min.pct = 0.5, logfc.threshold = 0.25)

DimPlot(object = Endo_rna, label = TRUE, split.by = "TreatmentGroup") 
ggsave("Endo_Dimplot_Split.png")

FeaturePlot(Endo_rna, features = c("Cxcl12","Ace"), label = T, split.by = "TreatmentGroup")
FeaturePlot(olson_rna, features = c("Tlr2","Tlr4"), label = T)
 
test <- WhichCells(Endo_rna, idents = 1)
test <- subset(Endo_rna, idents = c("1"))
Idents(test) <- "orig.ident"
table(Idents(test))
#ggsave("test.png", width = 40, height = 10)

# To understand the sc data in the context of the bulk data, are there just more endothelial cells (cxcl12+) and CFs (matricell+)
# relative to other cells after MI? Is this a game of proliferation regulation?
DimPlot(olson_rna)
# Let's rever to our cell-type clusters to see how composition changes
#DefaultAssay(olson_rna) <- 'SCT'
All_rna2 <- FindClusters(All_rna2, verbose = FALSE, resolution = 0.1)
DimPlot(All_rna2)

# # This tell us what the major cell types are 
# # Cluster 0 = CF = DCN, Col1A1
# # Cluster 1 = Endothelial Cells = CD36, FABP4
# # Cluster 2 = SMCs = Acta2, Rgs5
# # Cluster 3 = Macrophage/Monocyte = Lyz2, Spp1
# # Cluster 4 = Epicardial Cells =  C3 + MSLN
# # Cluster 5 = More Endothelial
# # Cluster 6 = Lymphocyte - Ccr7, Cytip
# # Cluster 7 = Macrophage 2 
# # Cluster 8 = Myoz2 = CM!

FeaturePlot(All_rna2, features = c("Ccr7","Myoz2"), label = T)
All_rna2 <- RenameIdents(All_rna2, `0` = "CF", `1` = "Endo", `2` = "SMC",'3'= "Mac", 
                         `4` = "EpiC", `5` = "Endo2", `6` = "Lymphocyte", `7` = "Mac2",
                         `8` = "CM")
All_rna2$CellTypes <- Idents(All_rna2)
table(Idents(All_rna2))

test <- subset(All_rna2, idents = c("CF"))
Idents(test) <- "orig.ident"
table(Idents(test))
Idents(All_rna2) <- "orig.ident"
table(Idents(All_rna2))

# This shows us proportion for CFs if you do division
# CFs:
# 1_1MI  1_1S 1_3MI  1_3S 8_1MI  8_1S 8_3MI  8_3S 
# 846  1111  1419  1014   279   303   592    98 

# Total 
# 1_1MI  1_1S 1_3MI  1_3S 8_1MI  8_1S 8_3MI  8_3S 
# 2226  2673  3471  2612  1015  1427  1672   541 

# Pct CF of total by sample is: 
# 1_1MI  1_1S 1_3MI  1_3S 8_1MI  8_1S 8_3MI  8_3S 
# 0.38  0.42  0.41   0.39  0.27  0.21  0.35  0.18

# Percentage decreases in nonregen MI, so unsure why Matricellular genes go up in bulk?

