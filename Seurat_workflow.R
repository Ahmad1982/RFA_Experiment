#########################################################################################
### Seurat workflow for all cells.						       
### Control example and RFA example were aggregated using the cellranger aggr pipeline.
#########################################################################################

library(Seurat)
# Load aggregated data
rfa.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")
rfa <- CreateSeuratObject(counts = rfa.data, project = "RFA_10x", min.cells = 3, min.features = 200)

# Add group information
rfa$group <- "RFA"
rfa@meta.data[1:10664,4] <- "CON"
rfa@meta.data[10663:10666,]
sce <- rfa

# QC and selecting cells for further analysis
sce <- PercentageFeatureSet(sce, pattern = "^Mt-", col.name = "percent.mt")
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce <- subset(sce, subset = nFeature_RNA > 400 & nFeature_RNA < 5500)

# Apply SCTransform normalization
sce <- SCTransform(sce, vars.to.regress = "percent.mt", verbose = FALSE)

# Perform dimensionality reduction by PCA and UMAP embedding
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims = 1:40)
sce <- FindNeighbors(sce, dims = 1:40)
sce <- FindClusters(sce)

# Visualization
DimPlot(sce, reduction = "umap", label = TRUE)

FeaturePlot(sce, cols = c("lightgrey", "red"), features = c("Foxp3", "Cd3d", "Cd8a", "Cd8b1", "Cd4", "S100a9", "G0s2", "Cxcr2", "Ncr1"), pt.size = 0.2, ncol = 3)

# Save data
save(sce, file = "10x_rfa_aggr_sce.umap.rds")

#########################################################################################
### Seurat workflow for lymphoid population.					       
### Clusters expressing Cd3d were extracted from aggregated samples.		       
#########################################################################################

library(Seurat)

# Load data
load("10x_rfa_aggr_sce.umap.rds")

# extract Cd3d positive clusters
sce <- subset(sce, idents = c("7", "8", "9", "10", "12", "14"))
VlnPlot(sce, features = "Cd3e", split.by = "group")
FeaturePlot(sce, cols = c("lightgrey", "red"),  features = "Cd3e")

# Redo dimensionality reduction by UMAP embedding
sce <- RunUMAP(sce, dims = 1:40)
sce <- FindNeighbors(sce, reduction = "pca", dims = 1:40)
sce <- FindClusters(sce, resolution = 1)

# Visualization
DimPlot(sce, reduction = "umap", label = TRUE)
FeaturePlot(sce, cols = c("lightgrey", "red"), features = c("Foxp3", "Cd3d", "Cd8a", "Cd8b1", "Cd4", "Mki67"))

# Save data
save(sce, file = "10x_rfa_Cd3e_aggr_sce.umap.rds")