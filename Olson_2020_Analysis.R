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
#To restore Rdata use:
#load("Olson_Global_Objects.RData")

#renv::init()
renv::snapshot()

# for sc RNA seq you need a directory for each sample that contains the following files
#barcodes.tsv
#features.tsv (previously called genes.tsv)
#matrix.mtx

# Since we have a lot of samples, we will create a vector of paths
rna_1_1MI.data <- Read10X("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scRNA/P1_1MI")  
rna_1_1S.data <- Read10X("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scRNA/P1_1S")
rna_1_3MI.data <- Read10X("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scRNA/P1_3MI")
rna_1_3S.data <- Read10X("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scRNA/P1_3S")
rna_8_1MI.data <- Read10X("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scRNA/P8_1MI")
rna_8_1S.data <- Read10X("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scRNA/P8_1S")
rna_8_3MI.data <- Read10X("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scRNA/P8_3MI")
rna_8_3S.data <- Read10X("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scRNA/P8_3S")
 
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
# Makeclusterss
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
write.csv(All_rna2.markers, file = "allmarkersbycluster.csv")

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


#The clusters we got were too small, so we should try again with less granular clustering
All_rna2 <- FindClusters(All_rna2, verbose = FALSE, resolution = 0.1)
All_rna2.markers <- FindAllMarkers(All_rna2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
All_rna2.markers %>% group_by(cluster) %>% slice_head(n = 2) -> Top2Markers

pdf("All_MFs2_TopMarkeronUMAP.pdf")
for (i in 1:length(Top2Markers$gene)){
  print(FeaturePlot(All_rna2, features = Top2Markers$gene[i], label = TRUE))}
dev.off()

# This tell us what the major cell types are 
# Cluster 0 = CF = DCN, Col1A1
# Cluster 1 = Endothelial Cells = CD36, FABP4
# Cluster 2 = SMCs = Acta2, Rgs5
# Cluster 3 = Macrophage/Monocyte = Lyz2, Spp1
# Cluster 4 = Epicardial Cells =  C3 + MSLN
# Cluster 5 = ? = Bace2, Plvap
# Cluster 6 = Ccr7, Cytip - resting T cells?
# Cluster 7 = some other macrophge?
# Cluster 8 = Myoz2 = CM!
# Cluster 9 = Alas2, Hbb-bh1 = 
#

#Use generalized lit markers to try to assign IDs to clusters
# Traditional Monocytes - CD14 + Lyz, Authors use ACE
# B cells - MS4A1
# NK Cells - GNLY, NKG7
# DCs - FCER1A, CST3
# Platelets - PPBP
# Naive CD4+ T - IL7R, CCR7
# Memory CD4+ T - IL7R, S100A4
# CD 8 T - CD8A
# Neutrophils - Ly6G

#~~~~~Clusters from the Paper~~~~~~~~~~~
# MHC2 Cluster - MMP12
# Monocyte Cluster - Ace
# Isg Cluster - Ifit1
# Timd4 Cluter - Timd4, Lyve1, Folr2
# CCR2 Cluster - least well defined, has Ccr2, Plac8, Il1b, C1qa


FeaturePlot(All_MFs2, features =c("Ace","Ifit1","Timd4","Mmp12","Cx3cr1","Igf1","H2-Eb1"), label = TRUE)
VlnPlot(All_MFs2, features = c("Ace","Ifit1","Timd4","Mmp12","H2-Eb1","Cd79a","Fabp4","Adgre1"),
        pt.size = 0)
VlnPlot(All_MFs2, features = c("Ace","Cd14","Cd68",'Cx3cr1',"Cd209a", "H2-Eb1"),
        pt.size = 0)
FeaturePlot(All_MFs2, features = c("Xcr1","Il1b"), label= TRUE)

#ID which cells are dividing 
pdf("Cells in Cell Cycle")
VlnPlot(All_MFs2, features = c("Xcr1","Irf8","Cdca3","Top2a"), pt.size = .5)
dev.off()

pdf("test")
FeaturePlot(All_MFs2, features = c("Ccl5","Fscn1","Ccr7"), label = TRUE)
VlnPlot(All_MFs2, features = c("Ccl5","Fscn1","Ccr7","Batf3",""))
dev.off()
# FROM this we know that:
# cluster 10 = Monocytes
# cluster 11 = Isg cluster
# cluster 1 = Timd4 cluster
# cluster __ = PROBABLY MHC2 because it has high H2-Eb1 and also MMP12
# cluster 12 is probably B cell contamination due to Cd79a
# cluster 7 = cDC because CD209a 
# cluster 13 = stromal cells referenced in paper


#Make a heatmap
pdf("PCHeatmapAll_MFs2.pdf", width = 20, height = 20)
DimHeatmap(object = All_MFs2, dims = 1:20,cells = 500, balanced = TRUE)
dev.off()

#Determine Dimensionality and make Jackstraw and ElbowPlots
#All_MFs2 <- JackStraw(object = All_MFs2, num.replicate = 100) 
#pdf("All_MFs2_ElbowPlot.pdf")
#ElbowPlot(object = All_MFs2)
#dev.off()

#All_MFs2 <- ScoreJackStraw(All_MFs2, dims = 1:20)
#pdf("All_MFs2_JackstrawPlot.pdf")
#JackStrawPlot(All_MFs2, dims = 1:20)
#dev.off()

#Make Clusters alternatively with pca instead of UMAP, not used
# All_MFs2 <- FindNeighbors(All_MFs2, reduction = "pca", dims = 1:20, verbose = TRUE, force.recalc = TRUE)
# All_MFs2 <- FindClusters(object = All_MFs2, resolution = 0.8)
# All_MFs2 <- RunTSNE(object = All_MFs2, dims = 1:20, resolution = 0.8)
# DimPlot(object = All_MFs2, reduction = 'tsne', group.by = 'orig.ident',label = 'TRUE')
# 
# # find all markers of cluster Control vs MI (over all clusters), not used in publication
# control_vs_MI_markers <- FindAllMarkers(All_MFs2, min.pct = 0.25, test.use = "MAST")
# DoHeatmap(All_MFs2)
# test <- (FetchData(object = All_MFs2,vars = c("ident","orig.ident","Klf4")))
# head(test)
# FeaturePlot(object = All_MFs2, features = c("Ly6c1","Timd4","H2-Eb1","Ccr2","Irf7","Ace","Lyve1","Plac8","Ifit1","Ifit3"), label = "TRUE")
# #get markers for each cluster
# control_vs_MI_markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20
# DoHeatmap(All_MFs2, features = top20$gene)
# VlnPlot(All_MFs2, features = c( "Retnla","Lyve1", "Timd4", "Cd163", "Folr2", "Klf2", 
#                                                  "Pf4", "F13a1", "Cd36", "Ccr2", "Fcgr1", "Adgre1", "C1qa",
#                                                  "Spp1", "Ms4a7"))



#Save seurat object
#saveRDS (All_MFs2, file = "All_MFs2.rds")

#Save all genes for easy lookup
all_genes <- as.data.frame(all.genesMF)

#Subset MI and Control sets of cells and then recalculate variable genes, scale, run PCA, and then cluster
All_MFs2$CellType <- Idents(All_MFs2)
Idents(All_MFs2) <- "orig.ident"
Cont_MFs2 <- subset(All_MFs2, idents = "Cont_MFs")
MI_MFs2 <- subset(All_MFs2, idents = "MI_MFs")
Cont_MFs2 <- FindVariableFeatures(Cont_MFs2)
MI_MFs2 <- FindVariableFeatures(MI_MFs2)
Cont_MFs2 <- SCTransform(Cont_MFs2, vars.to.regress = "percent.mt", verbose = TRUE)
MI_MFs2 <- SCTransform(MI_MFs2, vars.to.regress = "percent.mt", verbose = TRUE)
Cont_MFs2 <- RunPCA(Cont_MFs2, verbose = FALSE)
#Cont_MFs2 <- JackStraw(object = Cont_MFs2, num.replicate = 100) 
#Cont_MFs2 <- ScoreJackStraw(Cont_MFs2, dims = 1:20)
#pdf("Control MFs Dimensions.pdf")
#ElbowPlot(object = All_MFs2)
#JackStrawPlot(Cont_MFs2, dims = 1:20)
#dev.off()
Cont_MFs2 <- RunUMAP(Cont_MFs2, dims = 1:30, verbose = FALSE)
Cont_MFs2 <- RunTSNE(Cont_MFs2, dims = 1:30, verbose = FALSE)
Cont_MFs2 <- FindNeighbors(Cont_MFs2, dims = 1:30, verbose = FALSE)
Cont_MFs2 <- FindClusters(Cont_MFs2, verbose = FALSE)
MI_MFs2 <- RunPCA(MI_MFs2, verbose = FALSE)
MI_MFs2 <- RunUMAP(MI_MFs2, dims = 1:30, verbose = FALSE)
#MI_MFs2 <- JackStraw(object = MI_MFs2, num.replicate = 100) 
#MI_MFs2 <- ScoreJackStraw(MI_MFs2, dims = 1:20)
#pdf("MI MFs Dimensions.pdf")
#ElbowPlot(object = MI_MFs2)
#JackStrawPlot(MI_MFs2, dims = 1:20)
#dev.off()
MI_MFs2 <- FindNeighbors(MI_MFs2, dims = 1:30, verbose = FALSE)
MI_MFs2 <- FindClusters(MI_MFs2, verbose = FALSE)

#Visualize New Clusters
pdf("Cont and MI Clusters.pdf")
DimPlot(Cont_MFs2, label = TRUE, group.by = "ident") + NoLegend()
DimPlot(MI_MFs2, label = TRUE, group.by = "ident") + NoLegend()
dev.off()

#Find Markers for Each Cluster under each condition
Cont_MFs2.markers <- FindAllMarkers(Cont_MFs2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Cont_MFs2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) -> ContMFsTop2Markers

MI_MFs2.markers <- FindAllMarkers(MI_MFs2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MI_MFs2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) -> MIMFsTop2Markers


#ID Control Clusters 
# 0 = Timd4
# 1 = MHC2 - top2 are H2-Eb1, Cd74
# 2 = CCR2 - C5ar1, Ccl12
# 3 = Likely Same as Timd4 group - needs to be combined, but slightly older in pseudotime?
# 4 = Monocytes - has Ace, Plac8
# 5 = Proliferating Cells Cdca3, Top2a
# 6 = Probably APC group - has Fn1
# 7 = B cells - Cd79a
# 8 = Isg Cells - Ifit1
# 9 = Stromal Cells - Fap4

FeaturePlot(Cont_MFs2, features = c("Itgax","Fn1"), label = TRUE)

# Merge Clusters 0 and 3 since they are both the Timd4+, Lyve1+ cells
Cont_MFs2 <- RenameIdents(object = Cont_MFs2, '3' ='0')

#Change the cluster names 
Cont_MFs2 <- RenameIdents(Cont_MFs2, `0` = "Timd4", `1` = "MHC2", `2` = "CCR2", 
                          `4` = "Monocyte", `5` = "Prolif", `6` = "APC", `7` = "B Cells",
                          `8` = "Isg", `9` = "Stromal")
pdf("Cont MFs with Cell Types")
DimPlot(Cont_MFs2, label = TRUE)
dev.off()

#saveRDS (Cont_MFs2, file = "Cont_MFs2.rds")

# ID MI Clusters
# 0 = 
# 1 = 
# 2 = 
# 3 = MHC2 - Cd74, H2-Eb1
# 4 = Monocyte - Ms4a7, Fcrls, Cd14
# 5 = Timd4 - Folr2, F13a1, TGFBR2 - works with Science 2012 paper
# 6 = APCs? - Arg1, Spp1 (also could be M1ish cells?), Thbsp1
# 7 = Proliferating Cells - Stmn1, Hist1h2ap
# 8 = Isg - Ifit3, Isg15
# 9 = DCs - Cd209a, Napsa
# 10 = 
# 11 = 
# 12 = Stromal - Sparc, Fabp4

FeaturePlot(MI_MFs2, features = c("Timd4", "Trem2"), label = TRUE)
DotPlot(MI_MFs2, features = c("Ly6c2","Cd72", "Sirpb1a", "Dpp4","Ly6c1","Spp1"))
DotPlot(Cont_MFs2, features = c("Ly6c2","Cd72", "Sirpb1a", "Dpp4","Ly6c1","Spp1"))
DotPlot(All_MFs2, features = c("Ly6c2","Cd72", "Sirpb1a", "Dpp4","Ly6c1","Spp1"))

MI_MFs2 <- RenameIdents(MI_MFs2, `0` = "MI 11", `1` = "MI 10", `2` = "2",'3'= "MHC2", 
                        `4` = "Monocyte", `5` = "Timd4", `6` = "APC/DCs", `7` = "Prolif",
                        `8` = "Isg", `9` = "DCs/MI7", '10'="10",'11'="11",'12'="Stromal")

#saveRDS (MI_MFs2, file = "MI_MFs2.rds")
# Examples of Plot Types
# DotPlot(Cont_MFs2, features = c("Timd4", "H2-Eb1"))
# RidgePlot(Cont_MFs2, features = c("Timd4","H2-Eb1"))
# # CellScatter(Cont_MFs2,"Cont_AAACGGGGTCATGCCG", "Cont_AAAGATGAGATATGCA", features = c("Ace", "Cd14","Gapdh"))
# pdf("Biased Search for Rene", width = 20, height = 20)
# FeaturePlot(All_MFs2, features = c("Lmna","Syne1","Syne2", "Sun1","Sun2","Yap1","Wwtr1","Piezo1"), label = TRUE, ncol=2)
# FeaturePlot(Cont_MFs2, features = c("Lmna","Syne1","Syne2", "Sun1","Sun2","Yap1","Wwtr1","Piezo1"), label = TRUE)
# FeaturePlot(MI_MFs2, features = c("Lmna","Syne1","Syne2", "Sun1","Sun2","Yap1","Wwtr1","Piezo1"), label = TRUE)
# VlnPlot(All_MFs2, features = c("Lmna","Spp1","Timd4","Mpeg1"))
# dev.off()


# Show Timd4 cluster from Cont in the DimPlot of All Cells
ContTimd4Cells <- WhichCells(Cont_MFs2, idents = "Timd4")
MITimd4Cells <- WhichCells(MI_MFs2, idents = "Timd4")
DimPlot(All_MFs2, cells.highlight = c("ContTimd4Cells",""), label = TRUE )
FeaturePlot(MI_MFs2, features = c("Tnf","Timd4","Tgfb1","TdTomato","Il1b"), label = TRUE)
DoHeatmap(All_MFs2, features = c("Tnf","Timd4","Tgfb1","TdTomato","Il11ra1"))
VlnPlot(MI_MFs2, features = "Il1b")

# Primary Differential Cytokine is Il1b after MI
pdf("Cytokines from All and MI.pdf", width = 15, height = 15)
FeaturePlot(All_MFs2, features = c("Il1b", "Il4", "Il11","Il10","Il6","Tnf","Tgfb1"), label = TRUE)
FeaturePlot(MI_MFs2, features = c("Il1b", "Il4", "Il11","Il10","Il6","Tnf","Tgfb1"), label = TRUE)
VlnPlot(MI_MFs2, features = c("Il1b","Tgfb1"))
dev.off()

#Primary Differential Cytokine Receptor is Il10ra after MI
pdf("Cytokine Receptors from MI.pdf", width = 10, height = 15)
FeaturePlot(MI_MFs2, features = c("Il1r1", "Il1r2","Il6ra","Il10ra","Il10rb"), label = TRUE)
dev.off()

# Cd44/RHAMM/MMP3/MMP9 Expression after MI
pdf("Cd44 and HA signals.pdf", width = 10, height = 10)
FeaturePlot(MI_MFs2, features = c("Cd44", "Hmmr","Mmp3","Mmp9"), label = TRUE)
dev.off()

#ECM Genes after MI
pdf("ECM Mac Genes after MI.pdf", width = 10, height = 15)
FeaturePlot(MI_MFs2, features = c("Fn1", "Sparc","Lox","Loxl1","Loxl2","Loxl3"), label = TRUE)
dev.off()

# Chemokines after MI
pdf("Chemokines after MI.pdf", width = 15, height = 15)
FeaturePlot(MI_MFs2, features = c("Ccl2", "Cxcl1","C5ar1","C5ar2","Ccl5","Ccl7","Cxcl10","Cxcl12"), label = TRUE)
dev.off()

# Chemokine Receptors after MI
pdf("Chemokine Receptors after MI.pdf", width = 15, height = 15)
FeaturePlot(MI_MFs2, features = c("Ccr2", "Ccr5","Cxcr4","Cxcr3","TdTomato","Gpr35","Ccr1"), label = TRUE)
dev.off()

# Chemokine Receptors after MI
pdf("Chemokine Receptors after MI Chemo.pdf", width = 15, height = 15)
FeaturePlot(MI_MFs2, features = c("Ccr2", "Ccr5","Cxcr4","Cxcr3","Cxcr2","TdTomato","Gpr35","Ccr1"), label = TRUE)
dev.off()

# Hypoxia Genes
pdf("Hypoxia Genes after MI.pdf")
FeaturePlot(MI_MFs2, features = c("Hif1a","Nppa"), label = TRUE)
dev.off()

# Hypoxia Genes
pdf("SOCS Genes after MI.pdf", width = 15, height = 15)
FeaturePlot(MI_MFs2, features = c("Socs1","Socs2","Socs3","Socs4","Socs5","Socs6","Socs7","Cish","Stat1","Stat3"), label = TRUE)
#FeaturePlot(ChemoSeurat, features = c("Socs1","Socs2","Socs3","Socs4","Socs5","Socs6","Socs7","Cish","Stat1","Stat3"), label = TRUE)
dev.off()


#~~~~~~~~~ Create New Cluters Based on CCR2, Timd4, and CXCR4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Save existing ident as cluster.ident
MI_MFs2$cluster_names = Idents(MI_MFs2)

chemocluster = FetchData(MI_MFs2,c("ident","orig.ident","Ccr2","Cxcr4","Timd4"))
Ccr2Cells = WhichCells(object = MI_MFs2, expression = Ccr2 > 1)
FeaturePlot(MI_MFs2, cells = Ccr2Cells, features = c("Ccr2"), label = TRUE)
Cxcr4Cells = WhichCells(object = MI_MFs2, expression = Cxcr4 >1)
length(table(Ccr2Cells)) #1783
FeaturePlot(MI_MFs2, cells = Cxcr4Cells, features = c("Cxcr4"), label = TRUE)
length(table(Cxcr4Cells)) #728
Timd4Cells = WhichCells(object = MI_MFs2, expression = Timd4 >0.5)
FeaturePlot(MI_MFs2, cells = Timd4Cells, features = c("Timd4"), label = TRUE)
length(table(Timd4Cells)) #348
MI_MFs2 <- SetIdent(object = MI_MFs2, value = "Other")
MI_MFs2 <- SetIdent(object = MI_MFs2, cells = Ccr2Cells, value = "Ccr2Cells")
MI_MFs2 <- SetIdent(object = MI_MFs2, cells = Cxcr4Cells, value = "Cxcr4Cells")
MI_MFs2 <- SetIdent(object = MI_MFs2, cells = Timd4Cells, value = "Timd4Cells")
pdf("Chemo_on_MI.pdf")
FeaturePlot(MI_MFs2, features = c("Ccr2","Cxcr4","Timd4"), label = TRUE)
dev.off()

#Create new seurat with 60% of cells that bin niceley with expression measured
ChemoSeurat <- subset(MI_MFs2, idents = c("Ccr2Cells","Cxcr4Cells","Timd4Cells"))
# Reset MI_MFs2 Idents to clusters
Idents(MI_MFs2) <- MI_MFs2$cluster_names
#recluster and bin cells from new ChemoSeurat
ChemoSeurat <- FindVariableFeatures(ChemoSeurat)


#ChemoSeurat <- JackStraw(object = ChemoSeurat, num.replicate = 100) 
#ChemoSeurat <- ScoreJackStraw(ChemoSeurat, dims = 1:20)
#pdf("Chemo MFs Dimensions.pdf")
#ElbowPlot(object = ChemoSeurat)
#JackStrawPlot(ChemoSeurat, dims = 1:20)
#dev.off()

ChemoSeurat <- RunUMAP(ChemoSeurat, dims = 1:30, verbose = FALSE)

# With the new Seurat update, you don't need to SCTTransform again, BUT you can't do Jackstraw is subsetting
#ChemoSeurat <- SCTransform(ChemoSeurat, vars.to.regress = "percent.mt", verbose = TRUE)
#ChemoSeurat <- RunPCA(ChemoSeurat, verbose = FALSE)

ChemoSeurat <- FindNeighbors(ChemoSeurat, dims = 1:20, verbose = TRUE)
ChemoSeurat <- FindClusters(ChemoSeurat, resolution = 0.8, verbose = TRUE)

pdf("ID_Chemo_Clusters.pdf")
FeaturePlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4"), label = TRUE)
VlnPlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4"))
dev.off()

#Merge Clusters to obtain 3 main clusters: CCR2, CXCR4, TIMD4
ChemoSeurat <- RenameIdents(object = ChemoSeurat, "4" = "Timd4","5"="Ccr2","2"="Ccr2",
                            "9"="Ccr2","6"="Cxcr4","7"="Ccr2","8"="Ccr2","1"="Ccr2",
                            "0"="Ccr2","3"="Cxcr4")
#Groups 3 and 5 seemed to be double positive or ambiguous
pdf("Chemo_New_IDs.pdf")
FeaturePlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4"), label = TRUE)
VlnPlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4"))
DimPlot(ChemoSeurat, label = TRUE, group.by = "ident") + NoLegend()
DimPlot(ChemoSeurat, label = TRUE, group.by = "orig.ident") + NoLegend()
dev.off()

#Find DEGS Between the three clusters
ChemoSeurat.markers <- FindAllMarkers(ChemoSeurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ChemoSeurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) -> ChemoTop2Markers
FeaturePlot(ChemoSeurat, features = c("Gpnmb","Timd4","Tgfbr1"), label = TRUE)

#Get CXCR4 Markers and then Run through GO Term Screen
CXCR4_Markers <- ChemoSeurat.markers[c(204:328),(7)]
write.csv(CXCR4_Markers, file = "CXCR4 Markers for GO.csv")

pdf("ChemoPlots.pdf", width = 20, height = 20)
FeaturePlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4","Syngr1","Olfml3","Tgfbr1","S1pr5","Hif1a"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(ChemoSeurat, features = c("Ccr2","Cxcr4","Timd4","Ctsd","Ctsb","Ctsz","Ctsl","Ccl24"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(ChemoSeurat, features = c("Il1b", "Il4", "Il11","Il10","Il6","Tnf","Tgfb1","Pdgfa","Pdgfb","Il1rn","Igf2bp3"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(ChemoSeurat, features = c("Tnc", "Sparc", "Spp1","Postn","Thbs1","Thbs2","Ogn","Hspg2","Robo2","Robo3","Robo4","Ace"),
            label = TRUE, cols = c("lightblue","red"))
VlnPlot(ChemoSeurat, features = c("Spp1", "Col1A1","Col4a4","Igf1","Tgfbr1","Ace"), slot = "data")
dev.off()


pdf("MI_Cont_Chemo_Overlaps.pdf")
FeaturePlot(Cont_MFs2, features = c("Ccr2","Cxcr4","Timd4"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(MI_MFs2, features = c("Ccr2","Cxcr4","Timd4"),
            label = TRUE, cols = c("lightblue","red"))
DimPlot(All_MFs2, label=TRUE, group.by = "orig.ident")
FeaturePlot(All_MFs2, features = c("Ccr2","Cxcr4","Timd4"),
            label = TRUE, cols = c("lightblue","red"))
dev.off()

pdf("Dick_Plots_for_Adam.pdf")
FeaturePlot(MI_MFs2, features = c("Ccr2"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(MI_MFs2, features = c("Cxcr4"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(MI_MFs2, features = c("Timd4"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(MI_MFs2, features = c("Ccr5"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(MI_MFs2, features = c("Ccr1"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(MI_MFs2, features = c("Socs3"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(ChemoSeurat, features = c("Ccl24"),
            label = TRUE, cols = c("lightblue","red"))
FeaturePlot(ChemoSeurat, features = c("Spp1"),
            label = TRUE, cols = c("lightblue","red"))
dev.off()

# Here we will save all the objects so that they can be loaded later
save.image(file = "Dick_Global_Objects.RData")
load("Dick_Global_Objects.Rdata")
