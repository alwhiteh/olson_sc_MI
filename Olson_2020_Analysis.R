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



# scRNA -------------------------------------------------------------------

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
#All_rna2 <- readRDS(file = "scRNA_data.RDS")

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


# Compare RNA-seq Values with Dick SC Data -----------------------------------------------------------------
#Let's create a subsetted seurat object that only contains mono/mac
mac_rna <- subset(All_rna2, idents = c("MonoMac1","MonoMac2","MonoMac3"))
DimPlot(mac_rna)


#Run UMAP, make clusters
mac_rna <- SCTransform(mac_rna, vars.to.regress = "percent.mt", verbose = TRUE) #Correcting for nUMI is default YES
mac_rna <- RunPCA(mac_rna, verbose = FALSE)
mac_rna <- RunUMAP(mac_rna, dims = 1:30, verbose = FALSE)
mac_rna <- FindNeighbors(mac_rna, dims = 1:30, verbose = FALSE) 
mac_rna <- FindClusters(mac_rna, verbose = FALSE, resolution = 0.1)
mac_rna.markers <- FindAllMarkers(mac_rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # this takes a while
mac_rna.markers %>% group_by(cluster) %>% slice_head(n = 2) -> Top2MacMarkers

pdf("Macs_TopMarkeronUMAP.pdf")
for (i in 1:length(Top2MacMarkers$gene)){
  print(FeaturePlot(mac_rna, features = Top2MacMarkers$gene[i], label = TRUE))}
dev.off()

Idents(mac_rna) <- "orig.ident"
DimPlot(mac_rna)
mylevels <- c("MI","Sham","MI","Sham","MI","Sham","MI","Sham")
Idents(mac_rna) <- mylevels
mac_rna$infarct <- Idents(mac_rna)
DimPlot(mac_rna)
mylevels <- c("1","1","1","1","8","8","8","8")
Idents(mac_rna) <- mylevels
mac_rna$age <- Idents(mac_rna)
DimPlot(mac_rna)
mylevels <- c("1","1","3","3","1","1","3","3")
Idents(mac_rna) <- mylevels
mac_rna$timepoint <- Idents(mac_rna)
DimPlot(mac_rna)
#age, infarction, and time alone do not yield clusters, but together they do (i.e. orig.ident)
Idents(mac_rna) <- "SubCellTypes"
FeaturePlot(mac_rna, features = c("Cxcr4","Ccr2","Timd4"), label = T)

# Check CD14, CD16, and CD11b expression
FeaturePlot(mac_rna, features = c("Cd14","Fcgr3","Itgam"), label = T)
#Maybe cluster 1 is monocytes? 

#Save seurat object
saveRDS(mac_rna, file = "mac_rna.rds")


# Try to integrate with Harmony instead

#Let's read in the sc RNA seq from Dick et al and see how the cells cluster in aggregate.
dick <- readRDS(file='/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/Dick_Human_SSMac/Alex_Analysis_Dick_SS/All_MFs2.rds')
#chemo_seurat <- readRDS(file = '')
#agg_mac_rna <-  merge(mac_rna, y = dick,
#                      add.cell.ids = c("olson","dick"),
#                      project = "Aggregate_data")

agg.list <- list(dick, mac_rna)
agg.list <- lapply(X = agg.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

mac.anchors <- FindIntegrationAnchors(object.list = agg.list,dims = 1:20)
agg_mac_rna <- IntegrateData(anchorset = mac.anchors, dims = 1:20)

#DefaultAssay(agg_mac_rna) <- "integrated"
#agg_mac_rna <- ScaleData(agg_mac_rna, verbose = FALSE)
agg_mac_rna <- FindVariableFeatures(agg_mac_rna, selection.method = "vst", nfeatures = 2000)
agg_mac_rna <- RunPCA(agg_mac_rna, npcs = 30, verbose = FALSE)
agg_mac_rna <- RunUMAP(agg_mac_rna, dims = 1:30, verbose = FALSE)
agg_mac_rna <- FindNeighbors(agg_mac_rna, dims = 1:30, verbose = FALSE) 
agg_mac_rna <- FindClusters(agg_mac_rna, verbose = FALSE, resolution = 0.1)
agg_mac_rna.markers <- FindAllMarkers(agg_mac_rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # this takes a while
agg_mac_rna.markers %>% group_by(cluster) %>% slice_head(n = 2) -> Top2MacAggMarkers

Idents(agg_mac_rna) <- "orig.ident"
DimPlot(agg_mac_rna)




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
#test <- H5File$new("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scATAC/GSE153479_filtered_peak_bc_matrix.h5")

atac_counts <-Read10X_h5(filename = "/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scATAC/GSE153479_filtered_peak_bc_matrix.h5")


atac_metadata <- read.csv(file = "/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scATAC/GSE153479_singlecell.csv",
                          header = T, row.names =1)

chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(":","-"), genome = 'mm10',
                                    fragments = '/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/olson_ss_mi/olson_sc_MI/scATAC/GSE153479_fragments.tsv.gz',
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

pdf("test.pdf", width = 15, height = 40)
FeaturePlot(
  object = olson_rna,
  features = c("Nkx2-5", "Ets1","Tcf21","Ebf1"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  label = T)
dev.off()

DefaultAssay(olson_sc_atac) <- 'peaks'

pdf("coverage.pdf")
CoveragePlot(
  object = olson_sc_atac,
  region = "Cd8a",
  extend.upstream = 40000,
  extend.downstream = 20000
)
dev.off()
# We need to define ATAC clusters







# Integrate scATAC and scRNA ----------------------------------------------

# I don't think the original authors actually did this... 
# Might explain why the cells aren't overlapping very well

# We may need to use really broad clustering (i.e. cell type level) in order for this to work


# Make sure scRNA's active assay is RNA (and not SCT)
#olson_rna <- readRDS(file = "scRNA_data.RDS")
#olson_rna <- All_rna2
#olson_sc_atac <- readRDS(file = "sc_atac.RDS")
olson_rna_gran <- readRDS(file = "scRNA_granular.RDS")

DefaultAssay(olson_rna_gran) <- 'RNA' # Use the non-normalized values 
#olson_rna_gran <- FindVariableFeatures(
#  object = olson_rna,
#  nfeatures = 10000)

transfer.anchors <- FindTransferAnchors(
  reference = olson_rna_gran,
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
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = olson_sc_atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

plot1 <- DimPlot(object = olson_sc_atac,
    group.by = 'predicted.id',
    label = TRUE,
    repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot2 <- DimPlot(object = olson_sc_atac,
              group.by = 'predicted.id',
              label = TRUE,
              repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')







