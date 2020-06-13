################################################################################
### Monocle workflow for CD8 cluster.
################################################################################

library(monocle)
library(Seurat)

# Load data
load("10x_rfa_Cd3e_aggr_sce.umap.rds")

# extract CD8 clusters
sce <- subset(sce, idents = c("0", "2", "3", "5", "7", "8", "9", "10", "14"))

# Store data in a CellDataSet object
data <- as(as.matrix(sce@assays$RNA@counts), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = sce@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds@expressionFamily
head(pData(cds))
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"

# Add cluster annotation
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),c("0" = 'Cd8_s1',"2" = 'Cd8_s2',"3" = 'Cd8_s3',"5" = 'Cd8_s4',"7" = 'Cd8_s5',"8" = 'Cd8_s6',"9" = 'Cd8_s7',"10" = 'Cd8_s8',"14" = 'Cd8_s9'))
# cell_type_color <- c("Fibroblasts" = "#E088B8","NK" = "#46C7EF","Macrophage" = "#EFAD1E","EC" = "#8CB3DF")
table(pData(cds)$cell_type2)

# trajectory progress
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~cell_type2")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

# reduce data dimensionality
ds <- reduceDimension(
  cds,
  max_components = 2,
  method = 'DDRTree')

# order cells along the trajectory
cds <- orderCells(cds)
save(cds, file = "Cd8_monocle2.rds")
load("Cd8_monocle2.rds")
plot_cell_trajectory(cds, color_by = "cell_type2")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State")
GM_state <- function(a){
  if (length(unique(pData(a)$State)) > 1){
    T0_counts <- table(pData(a)$State, pData(a)$cell_type2)[,"Cd8_s4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds_1 <- orderCells(cds, root_state = GM_state(cds))
plot_cell_trajectory(cds_1, color_by = "cell_type2")+ scale_color_manual(values=c("#F8766C", "#FF7F00", "#B2A102", "#FFC800", "#00C087", "#44B500", "#01BCD6", "#01B3F2","#F166E8"))
plot_cell_trajectory(cds_1, color_by = "Pseudotime")


################################################################################
### Monocle workflow for macrophage cluster.
################################################################################

library(monocle)
library(Seurat)

# Load data
load("10x_rfa_aggr_sce.umap.rds")

# extract macrophage clusters
sce <- subset(sce, idents = c("0", "1", "2", "3", "13"))

# Store data in a CellDataSet object
data <- as(as.matrix(sce@assays$RNA@counts), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = sce@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds@expressionFamily
head(pData(cds))
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"

# Add cluster annotation
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),c("0" = 'Mac_s1',"1" = 'Mac_s2',"2" = 'Mac_s3',"3" = 'Mac_s4',"13" = 'Mac_s5'))
cell_type_color <- c("Mac_s1" = "#F8766C","Mac_s2" = "#AAAAAA","Mac_s3" = "#FF7F00","Mac_s4" = "#B2A102","Mac_s5" = "#D277FF")
table(pData(cds)$cell_type2)

# trajectory progress
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~cell_type2")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

#disp_table <- dispersionTable(cds)
#ordering_genes <- subset(disp_table, mean_expression >= 0.1)
#cds <- setOrderingFilter(cds, ordering_genes)

# reduce data dimensionality
ds <- reduceDimension(
  cds,
  max_components = 2,
  method = 'DDRTree')

# order cells along the trajectory
cds <- orderCells(cds)
colnames(pData(cds))
plot_cell_trajectory(cds, color_by = "cell_type2")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State")
GM_state <- function(a){
  if (length(unique(pData(a)$State)) > 1){
    T0_counts <- table(pData(a)$State, pData(a)$cell_type2)[,"Mac_s5"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds_1 <- orderCells(cds, root_state = GM_state(cds))
plot_cell_trajectory(cds_1, color_by = "cell_type2", cell_size = 1.5)+ scale_color_manual(values=c("#F8766C", "#AAAAAA", "#FF7F00", "#B2A102", "#D277FF"))
plot_cell_trajectory(cds_1, color_by = "Pseudotime")