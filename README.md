---
title: "NNT-dependent Cardiometabolic HFpEF - snRNA Sequencing Analysis"
# author: "Mark E. Pepin, MD, PhD, MS"
output:
  html_document:
    code_folding: hide
    keep_md: yes
    toc: yes
    toc_float: yes
header-includes:
- \usepackage{booktabs}
- \usepackage{longtable}
- \usepackage{array}
- \usepackage{multirow}
- \usepackage[table]{xcolor}
- \usepackage{wrapfig}
- \usepackage{float}
- \usepackage{colortbl}
- \usepackage{pdflscape}
- \usepackage{tabu}
- \usepackage{threeparttable}x
mainfont: Times
fontsize: 10pt
always_allow_html: yes
editor_options: 
  markdown: 
    wrap: 72
---



**Code Authors**: Mark E. Pepin, MD, PhD, MS **Contact**:
[mepepin\@bwh.harvard.edu](mailto:mepepin@bwh.harvard.edu){.email}\
**Institution**: Brigham and Women's Hospital \| Broad Institute of
Harvard and MIT\
**Location**: Boston, MA, USA

# Sample Pre-processing

Cellranger output data were filtered to remove of ambient RNA using
cellbender in full running mode. The output files, which constitute
feature-barcode matrices and cluster assignments, were then imported
into the R (4.3.1) environment using the Seurat package (5.0.1) for
downstream preprocessing, annotation, and differential expression. To
improve data quality, cells with fewer than 200 detected genes, as well
as those containing more than 5% mitochondrial gene content, were
excluded from downstream analysis.


``` r
start_time <- Sys.time()
library(Seurat)
Index <- openxlsx::read.xlsx("../1_Input/M020_Sample.Info_NvJ_mep.xlsx")
snRNA.dat_1 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[1],"_", Index$Background[1], "_", Index$MouseID[1], "/filtered_feature_bc_matrix/"))
snRNA.dat_2 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[2],"_", Index$Background[2], "_", Index$MouseID[2], "/filtered_feature_bc_matrix/"))
snRNA.dat_3 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[3],"_", Index$Background[3], "_", Index$MouseID[3], "/filtered_feature_bc_matrix/"))
snRNA.dat_7 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[4],"_", Index$Background[4], "_", Index$MouseID[4], "/filtered_feature_bc_matrix/"))
snRNA.dat_8 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[5],"_", Index$Background[5], "_", Index$MouseID[5], "/filtered_feature_bc_matrix/"))
snRNA.dat_9 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[6],"_", Index$Background[6], "_", Index$MouseID[6], "/filtered_feature_bc_matrix/"))
snRNA.dat_13 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[7],"_", Index$Background[7], "_", Index$MouseID[7], "/filtered_feature_bc_matrix/"))
snRNA.dat_14 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[8],"_", Index$Background[8], "_", Index$MouseID[8], "/filtered_feature_bc_matrix/"))
snRNA.dat_15 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[9],"_", Index$Background[9], "_", Index$MouseID[9], "/filtered_feature_bc_matrix/"))
snRNA.dat_19 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[10],"_", Index$Background[10], "_", Index$MouseID[10], "/filtered_feature_bc_matrix/"))
snRNA.dat_20 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[11],"_", Index$Background[11], "_", Index$MouseID[11], "/filtered_feature_bc_matrix/"))
snRNA.dat_21 <- Read10X(data.dir = paste0("../1_Input/snRNA/", Index$Treatment[12],"_", Index$Background[12], "_", Index$MouseID[12], "/filtered_feature_bc_matrix/"))
# Initialize the Seurat object with the raw (non-normalized data).
snRNA_seurat_1 <- CreateSeuratObject(counts = snRNA.dat_1, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_1[["SampleID"]] <- Index$MouseID[1]
snRNA_seurat_1[["Background"]] <- Index$Background[1]
snRNA_seurat_1[["Treatment"]] <- Index$Treatment[1]
snRNA_seurat_2 <- CreateSeuratObject(counts = snRNA.dat_2, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_2[["SampleID"]] <- Index$MouseID[2]
snRNA_seurat_2[["Background"]] <- Index$Background[2]
snRNA_seurat_2[["Treatment"]] <- Index$Treatment[2]
snRNA_seurat_3 <- CreateSeuratObject(counts = snRNA.dat_3, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_3[["SampleID"]] <- Index$MouseID[3]
snRNA_seurat_3[["Background"]] <- Index$Background[3]
snRNA_seurat_3[["Treatment"]] <- Index$Treatment[3]
snRNA_seurat_7 <- CreateSeuratObject(counts = snRNA.dat_7, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_7[["SampleID"]] <- Index$MouseID[4]
snRNA_seurat_7[["Background"]] <- Index$Background[4]
snRNA_seurat_7[["Treatment"]] <- Index$Treatment[4]
snRNA_seurat_8 <- CreateSeuratObject(counts = snRNA.dat_8, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_8[["SampleID"]] <- Index$MouseID[5]
snRNA_seurat_8[["Background"]] <- Index$Background[5]
snRNA_seurat_8[["Treatment"]] <- Index$Treatment[5]
snRNA_seurat_9 <- CreateSeuratObject(counts = snRNA.dat_9, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_9[["SampleID"]] <- Index$MouseID[6]
snRNA_seurat_9[["Background"]] <- Index$Background[6]
snRNA_seurat_9[["Treatment"]] <- Index$Treatment[6]
snRNA_seurat_13 <- CreateSeuratObject(counts = snRNA.dat_13, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_13[["SampleID"]] <- Index$MouseID[7]
snRNA_seurat_13[["Background"]] <- Index$Background[7]
snRNA_seurat_13[["Treatment"]] <- Index$Treatment[7]
snRNA_seurat_14 <- CreateSeuratObject(counts = snRNA.dat_14, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_14[["SampleID"]] <- Index$MouseID[8]
snRNA_seurat_14[["Background"]] <- Index$Background[8]
snRNA_seurat_14[["Treatment"]] <- Index$Treatment[8]
snRNA_seurat_15 <- CreateSeuratObject(counts = snRNA.dat_15, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_15[["SampleID"]] <- Index$MouseID[9]
snRNA_seurat_15[["Background"]] <- Index$Background[9]
snRNA_seurat_15[["Treatment"]] <- Index$Treatment[9]
snRNA_seurat_19 <- CreateSeuratObject(counts = snRNA.dat_19, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_19[["SampleID"]] <- Index$MouseID[10]
snRNA_seurat_19[["Background"]] <- Index$Background[10]
snRNA_seurat_19[["Treatment"]] <- Index$Treatment[10]
snRNA_seurat_20 <- CreateSeuratObject(counts = snRNA.dat_20, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_20[["SampleID"]] <- Index$MouseID[11]
snRNA_seurat_20[["Background"]] <- Index$Background[11]
snRNA_seurat_20[["Treatment"]] <- Index$Treatment[11]
snRNA_seurat_21 <- CreateSeuratObject(counts = snRNA.dat_21, project = "snRNA_NNT", min.cells = 3, min.features = 200)
snRNA_seurat_21[["SampleID"]] <- Index$MouseID[12]
snRNA_seurat_21[["Background"]] <- Index$Background[12]
snRNA_seurat_21[["Treatment"]] <- Index$Treatment[12]
# Initialize the Seurat object with the raw (non-normalized data).
snRNA_seurat_1[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_1, assay = "RNA", pattern = "^mt-")
snRNA_seurat_1 <- subset(snRNA_seurat_1, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_2[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_2, assay = "RNA", pattern = "^mt-")
snRNA_seurat_2 <- subset(snRNA_seurat_2, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_3[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_3, assay = "RNA", pattern = "^mt-")
snRNA_seurat_3 <- subset(snRNA_seurat_3, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_7[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_7, assay = "RNA", pattern = "^mt-")
snRNA_seurat_7 <- subset(snRNA_seurat_7, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_8[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_8, assay = "RNA", pattern = "^mt-")
snRNA_seurat_8 <- subset(snRNA_seurat_8, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_9[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_9, assay = "RNA", pattern = "^mt-")
snRNA_seurat_9 <- subset(snRNA_seurat_9, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_13[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_13, assay = "RNA", pattern = "^mt-")
snRNA_seurat_13 <- subset(snRNA_seurat_13, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_14[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_14, assay = "RNA", pattern = "^mt-")
snRNA_seurat_14 <- subset(snRNA_seurat_14, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_15[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_15, assay = "RNA", pattern = "^mt-")
snRNA_seurat_15 <- subset(snRNA_seurat_15, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_19[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_19, assay = "RNA", pattern = "^mt-")
snRNA_seurat_19 <- subset(snRNA_seurat_19, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_20[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_20, assay = "RNA", pattern = "^mt-")
snRNA_seurat_20 <- subset(snRNA_seurat_20, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_seurat_21[["percent.mt"]] <- PercentageFeatureSet(snRNA_seurat_21, assay = "RNA", pattern = "^mt-")
snRNA_seurat_21 <- subset(snRNA_seurat_21, subset = nFeature_RNA > 200 & percent.mt < 5)
snRNA_list<-list(snRNA_seurat_1, snRNA_seurat_2, snRNA_seurat_3, snRNA_seurat_7, snRNA_seurat_8, snRNA_seurat_9, snRNA_seurat_13, snRNA_seurat_14, snRNA_seurat_15, snRNA_seurat_19, snRNA_seurat_20, snRNA_seurat_21)
# normalize and identify variable features for each dataset independently
TAA.list <- lapply(X = snRNA_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
    x <- ScaleData(x)
    x <- SCTransform(x)
})
# Find most variable features across samples to integrate
integ_features <- SelectIntegrationFeatures(object.list = TAA.list, nfeatures = 3000)
merged_seurat <- merge(x = TAA.list[[1]],
		       y = TAA.list[2:length(TAA.list)],
		       merge.data = TRUE)
DefaultAssay(merged_seurat) <- "SCT"
# Manually set variable features of merged Seurat object
VariableFeatures(merged_seurat) <- integ_features
# Calculate PCs using manually set variable features
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)
library(harmony)
harmonized_seurat <- RunHarmony(merged_seurat, 
				group.by.vars = c("SampleID"), 
				reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:40)
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))
harmonized_seurat <- PrepSCTFindMarkers(harmonized_seurat, assay = "SCT", verbose = TRUE)
# markers <- FindAllMarkers(
#   object = harmonized_seurat,
#   assay = "SCT",
#   verbose = F
# )
# saveRDS(harmonized_seurat, file = "../1_Input/NNT_Integration_snRNA.rds")
```

# Unbiased Cell-Type Identification

Normalization and scaling were performed using the "LogNormalize" method
with a scaling factor of 10,000. Highly variable genes were identified
using the "vst" method, and the 10,000 genes with the highest variance
were retained for subsequent analysis. Unsupervised clustering via
Uniform Manifold Approximation and Projection (UMAP) was performed on
the scaled and variable gene expression data, and the top 30 principal
components were used to define Louvian clustering with the
"FindClusters" function (resolution = 0.2). Uniform Manifold
Approximation and Projection (UMAP) was employed for dimensionality
reduction and visualization of the clustered data. Cellular origins of
each sequenced nucleus were individually estimated using singleR
(2.6.0), which computes a Pearson correlation between the nuclear
transcriptome and those of purified cell types. We used the cardiac cell
atlas developed by Litviňuková et al.19 as the curated reference
dataset.


``` r
# Load necessary libraries
library(Seurat)
library(reticulate)
# Define the file paths
h5ad_file <- "../1_Input/Annotation/Global_lognormalised.h5ad"
seurat_rds_file <- "../1_Input/Annotation/Global_lognormalised_seurat.rds"
# Step 1: Set up Python environment and import anndata
# Adjust this if needed to use a specific Python environment
use_condaenv("myenv", required = TRUE)
# Load reticulate and anndata
anndata <- import("anndata")
# Step 2: Read the .h5ad file using anndata
ad <- anndata$read_h5ad(h5ad_file)
# Step 3: Extract relevant components from the AnnData object
# Extract the count matrix
counts <- t(ad$X)  # Transpose to match Seurat's row/column orientation
# Extract cell metadata (if available)
cell_metadata <- ad$obs
cell_metadata <- as.data.frame(cell_metadata)
# Extract feature (gene) metadata (if available)
feature_metadata <- ad$var
feature_metadata <- as.data.frame(feature_metadata)
# Ensure the row names of the count matrix and cell metadata are aligned
rownames(cell_metadata) <- colnames(counts)
# Step 4: Create a Seurat object from the extracted data
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = cell_metadata)
# Step 5: Save the Seurat object as an .rds file
saveRDS(seurat_obj, file = seurat_rds_file)
# Print a message to indicate completion
cat("Seurat object saved as:", seurat_rds_file, "\n")

###
library(Seurat)
library(scCustomize)
library(Nebulosa)
library(ggplot2)
library(ggtrace)
library(ggrepel)
snRNA.combined<-readRDS(file = "../1_Input/NNT_Integration_snRNA.rds")
Idents(snRNA.combined) <- "SCT_snn_res.0.2" # Change the identity of the clusters to cell types
UMAP_Treatment<-DimPlot(snRNA.combined,  label = T, split.by = "Treatment") + NoLegend()
UMAP_Background<-DimPlot(snRNA.combined,  label = T, split.by = "Background")
UMAP_Treatment + UMAP_Background
pdf(file = "../2_Output/UMAP_Clusters_snnres_0.2.pdf", height = 5, width = 9)
UMAP_Treatment
UMAP_Background
dev.off()
# Unbiased cluster identification
DEGs_Clusters<-FindAllMarkers(snRNA.combined, assay = "SCT")
write.csv(DEGs_Clusters, "../2_Output/DEGs_Clusters.csv")
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral2"))(paletteLength)
top5_markers <- Extract_Top_Markers(marker_dataframe = DEGs_Clusters, num_genes = 5, named_vector = FALSE,
    make_unique = TRUE)
pdf(file = "../2_Output/DotPlot_Clusters_top5DEGs.pdf", height = 10, width = 7)
Clustered_DotPlot(seurat_object = snRNA.combined, 
                  features = top5_markers, 
                  k = 15, 
                  colors_use_exp = myColor,
                  colors_use_idents = NA,
                  cluster_ident = F)
dev.off()
top30_markers <- Extract_Top_Markers(marker_dataframe = DEGs_Clusters, num_genes = 5, named_vector = FALSE,
    make_unique = TRUE)
pdf(file = "../2_Output/Heatmap_Clusters.pdf", height = 15, width = 10)
DoHeatmap(snRNA.combined, features = top30_markers, size = 3, disp.min = -2, disp.max = 2) + scale_fill_gradientn(colors = c("dodgerblue4", "white", "coral2"))
dev.off()

library(scRNAseq)
library(dplyr)
library(Seurat)
hESCs <- as.SingleCellExperiment(readRDS(file = "../1_Input/NNT_Integration_snRNA.rds"))
# reference dataset (formatted above)
sceM <- readRDS("../1_Input/Annotation/Global_lognormalised_seurat.rds")
sceM <- subset(sceM, subset = cell_or_nuclei == "Nuclei" & region=="LV" & modality=="snRNA")
# Count the number of each cell type
sceM@meta.data %>% 
  group_by(cell_type) %>%
  summarise(number_nuclei = length(cell_type))
# down-sample the reference dataset
metadata <- sceM@meta.data ## Extract metadata with cell type information
metadata$cell_names <- colnames(sceM)
random_sample <- metadata %>%
  group_by(!!sym("cell_type")) %>%
  sample_n(size = min(500, n()))  # Handles cases where there are fewer than 100 cells in a cell type
random_sample %>% # report the cell counts for down-sampled reference dataset
  group_by(cell_type) %>%
  summarise(number_nuclei = length(cell_type))
selected_cells <- random_sample$cell_names ## Extract cell names from the sample
seurat_sampled <- subset(sceM, cells = selected_cells) # Subset the Seurat object
#
sceM <- as.SingleCellExperiment(seurat_sampled)
rownames(sceM) <- stringr::str_to_title(rownames(sceM))
library(scuttle)
sceM <- logNormCounts(sceM)
sceM <- sceM[,!is.na(sceM$cell_type)]
# Annotation with singleR
library(SingleR)
pred.hesSingleRpred.hesc <- SingleR(test = hESCs, ref = sceM, assay.type.test="counts",
    labels = sceM$cell_type)
hESC_seur <- as.Seurat(hESCs, counts = "counts", data = NULL)
hESC_seur[["SingleR.labels"]] <- pred.hesSingleRpred.hesc$labels
# 
# # Or if `method="cluster"` was used:
# snRNA.combined[["SingleR.cluster.labels"]] <- 
#         pred.hesSingleRpred.hesc$labels[match(snRNA.combined[[]][["my.input.clusters"]], rownames(pred.hesSingleRpred.hesc))]
snRNA.combined[["SingleR.labels"]] <- pred.hesSingleRpred.hesc$labels
Idents(snRNA.combined) <- snRNA.combined$SingleR.labels
UMAP_Treatment<-DimPlot(snRNA.combined,  label = T, split.by = "Treatment") + NoLegend()
UMAP_Background<-DimPlot(snRNA.combined,  label = T, split.by = "Background")
UMAP_Treatment + UMAP_Background
pdf(file = "../2_Output/UMAP_Clusters_unbiased.annotation.pdf", height = 5, width = 9)
UMAP_Treatment
UMAP_Background
dev.off()
saveRDS(snRNA.combined, , file = "../1_Input/snRNA_unbiased.Annnotation_snRNA.rds")
# Annotation performance with SingleR
## heatmap of annotation scores
pdf("../2_Output/annotation.heatmap.pdf")
plotScoreHeatmap(pred.hesSingleRpred.hesc)
dev.off()
## Delta Distribution
plotDeltaDistribution(pred.hesSingleRpred.hesc, ncol = 3)
```

# Gene Markers to Validate Cell Type Identification


``` r
snRNA.combined <- readRDS(file = "../1_Input/snRNA_unbiased.Annnotation_snRNA.rds")
library(Seurat)
library(Nebulosa)
library(ggpubr)
EC <- c("Vwf", "Pecam1", "Cdh5")
Neuronal <- c("Plp1", "Nrxn1", "Nrxn3")
Fibroblast <- c("Dcn", "Gsn", "Pdgfra")
Adipocyte <- c("Gpam", "Fasn", "Lep") #, 
Mesothelial <- c("Msln", "Wt1") #, "Bnc1"
Pericyte <- c("Rgs5", "Abcc9", "Kcnj8")
SMC <- c("Myh11", "Acta2", "Tagln")
ACM <- c("Nppa","Myl4") #"Myl7",
VCM <- c("Myh7", "Myl2", "Fhl2")
Lymphoid <- c("Cd5")
Myeloid <- c("Cd14", "C1qa", "Cd68")
# # Plot density function
EC_density<-plot_density(snRNA.combined, EC, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Neuronal_density<-plot_density(snRNA.combined, Neuronal, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast_density<-plot_density(snRNA.combined, Fibroblast, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Adipocyte_density<-plot_density(snRNA.combined, Adipocyte, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Mesothelial_density<-plot_density(snRNA.combined, Mesothelial, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Pericyte_density<-plot_density(snRNA.combined, Pericyte, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
SMC_density<-plot_density(snRNA.combined, SMC, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
ACM_density<-plot_density(snRNA.combined, ACM, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
VCM_density<-plot_density(snRNA.combined, VCM, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Lymphoid_density<-plot_density(snRNA.combined, Lymphoid, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Myeloid_density<-plot_density(snRNA.combined, Myeloid, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
# 
pdf(file = "../2_Output/Cardiac_Density_plots.pdf", height = 3, width = 3)
EC_density
Neuronal_density
Fibroblast_density
Adipocyte_density
Mesothelial_density
SMC_density
ACM_density
VCM_density
Lymphoid_density
Myeloid_density
dev.off()

VCM_dens <- plot_density(snRNA.combined, "Myh7", pal = "magma") + theme(legend.position="none")
VSMC_dens
ACM_dens <- plot_density(snRNA.combined, "Nppa", pal = "magma") + theme(legend.position="none")
EC_dens <- plot_density(snRNA.combined, "Cdh5", pal = "magma") + theme(legend.position="none")
VSMC_dens <- plot_density(snRNA.combined, "Acta2", pal = "magma") + theme(legend.position="none")
Fibroblast_dens <- plot_density(snRNA.combined, "Pdgfra", pal = "magma") + theme(legend.position="none")
Adipocyte_dens <- plot_density(snRNA.combined, "Lep", pal = "magma") + theme(legend.position="none")
Neuronal_dens <- plot_density(snRNA.combined, "Plp1", pal = "magma") + theme(legend.position="none")
Lymphoid_dens <- plot_density(snRNA.combined, "Cd5", pal = "magma") + theme(legend.position="none")
Myeloid_dens <- plot_density(snRNA.combined, "Cd14", pal = "magma") + theme(legend.position="none")
# Cluster Identification using gene markers
pdf(file = "../2_Output/CellType_FeaturePlots.pdf", width = 11, height = 11)
  ggarrange(VCM_dens,
            ACM_dens,
            EC_dens,
            VSMC_dens,
            Fibroblast_dens,
            Adipocyte_dens,
            Neuronal_dens,
            Lymphoid_dens,
            Myeloid_dens,
            ncol = 3, nrow = 3)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = Marker_list, ncol = 3)
dev.off()
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = c(Mesothelial[1],Lymphoid[1], Myeloid[1]), split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = EC, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = Neuronal, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = Fibroblast, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = Mesothelial, split.by = "Treatment", ncol = 2)
# # FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = Pericyte, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = SMC, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = ACM, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = VCM, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = Lymphoid, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = Myeloid, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, features = Adipocyte, split.by = "Treatment", ncol = 2)
# FeaturePlot(snRNA.combined, reduction = "umap", label = T, split.by = "Background", features = "Nnt", ncol = 2)
dev.off()
## Use ridgeplots to identify bimodal gene marker distributions (enriched clusters)
pdf(file = "../2_Output/Celltype_RidgePlots.pdf", height = 10, width = 15)
RidgePlot(snRNA.combined,features = CM_genes, ncol = 2)
RidgePlot(snRNA.combined, features = EC, ncol = 2)
RidgePlot(snRNA.combined, features = Neuronal, ncol = 2)
RidgePlot(snRNA.combined, features = Fibroblast, ncol = 2)
RidgePlot(snRNA.combined, features = Mesothelial, ncol = 2)
RidgePlot(snRNA.combined, features = Pericyte, ncol = 2)
RidgePlot(snRNA.combined, features = SMC, ncol = 2)
RidgePlot(snRNA.combined, features = ACM, ncol = 2)
RidgePlot(snRNA.combined, features = VCM, ncol = 2)
RidgePlot(snRNA.combined, features = Lymphoid, ncol = 2)
RidgePlot(snRNA.combined, features = Myeloid, ncol = 2)
dev.off()
saveRDS(snRNA.combined, file = "../1_Input/snRNA_Clustering_snRNA.rds")
```

# Annotating Cells

Based on the cell identification above, clusters were renamed and saved
for further sub-cluster analysis.


``` r
library(Seurat)
snRNA.combined <- readRDS(file = "../1_Input/snRNA_Clustering_snRNA.rds")
Idents(snRNA.combined) <- snRNA.combined$SCT_snn_res.0.2
snRNA.combined <- RenameIdents(snRNA.combined,
             `0` = "EC",
             `1` = "Cardiomyocyte",
             `2` = "Fibroblast",
             `3` = "Mural_Cell",
             `4` = "Myeloid",
             `5` = "Fibroblast",
             `6` = "Fibroblast",
             `7` = "Mast_Cell",
             `8` = "Cardiomyocyte",
             `9` = "Mural_Cell",
             `10` = "Cardiomyocyte",
             `11` = "Lymphoid",
             `12` = "Mural_Cell",
             `13` = "Neural_Cell",
             `14` = "Cardiomyocyte")
snRNA.combined$CellType <- snRNA.combined@active.ident
snRNA.combined$CellType <- factor(snRNA.combined$CellType, levels = c("Cardiomyocyte", 
                                                                      "EC", 
                                                                      "Fibroblast",
                                                                      "Mural_Cell",
                                                                      "Lymphoid",
                                                                      "Myeloid",
                                                                      "Mast_Cell",
                                                                      "Neural_Cell"))
snRNA.combined <- SetIdent(snRNA.combined, value = "CellType")
UMAP_CellTypes<-DimPlot(snRNA.combined,  label = T, repel = T, split.by = "Treatment",label.size = 4,
                        cols = c("coral2", 
                                 "wheat", 
                                 "steelblue4",
                                 "deepskyblue3",
                                 "azure4",
                                 "goldenrod2",
                                 "tan2",
                                 "darkcyan")) + NoLegend()
UMAP_CellTypes
UMAP_CellTypes_bkg<-DimPlot(snRNA.combined,  label = T, repel = T, split.by = "Background",label.size = 4,
                        cols = c("coral2", 
                                 "wheat", 
                                 "steelblue4",
                                 "deepskyblue3",
                                 "azure4",
                                 "goldenrod2",
                                 "tan2",
                                 "darkcyan")) + NoLegend()
pdf(file = "../2_Output/UMAP_Cell.Clusters.pdf", height = 5, width = 7)
UMAP_CellTypes
UMAP_CellTypes_bkg
dev.off()
pdf(file = "../2_Output/UMAP_backgroud.pdf", height = 5, width = 7)
UMAP_CellTypes_bkg
dev.off()

Philipp <- c("Slc7a10", "Pycard","Il1b", "Nlrp3", "Nlrx1")
# Overlay these gene markers onto the UMAP to identify clusters
pdf("Philipp_Inflammasome.pdf")
FeaturePlot(snRNA.combined, reduction = "umap", label = T,features = Philipp, ncol = 2)
plot_density(snRNA.combined, Philipp, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
dev.off()

# Save file
saveRDS(snRNA.combined, file = "../1_Input/NNT_Labelling_snRNA.rds")
```

## UMAP and Cellular Proportions


``` r
# Cell-type Specific Differential Expression
library(ggplot2)
library(dplyr)
library(Seurat)
TAA.combined <- readRDS(file = "../1_Input/NNT_Labelling_snRNA.rds")
TAA.combined$CellType <- factor(TAA.combined$CellType, levels = c("Cardiomyocyte", 
                                                                      "EC", 
                                                                      "Fibroblast",
                                                                      "Mural_Cell",
                                                                      "Lymphoid",
                                                                      "Myeloid",
                                                                      "Mast_Cell",
                                                                      "Neural_Cell"))
# Show Integration Success between cohorts/datasets
p1 <- DimPlot(TAA.combined, reduction = "umap", group.by = "Treatment") + ggtitle(NULL)# Chou vs. Pepin
p1
# Show distribution of clustering across conditions
p2 <- DimPlot(TAA.combined, reduction = "umap", group.by = "SCT_snn_res.0.2", split.by = "Treatment", label = T, repel = TRUE) + NoLegend() + ggtitle(NULL)
p2
pdf(file = "../2_Output/UMAP_Clusters.pdf", height = 4, width = 7)
p2
dev.off()

p3 <- DimPlot(TAA.combined, reduction = "umap", group.by = "CellType", split.by = "Treatment", label = T, repel = TRUE,
              cols = c("coral2", 
                                 "wheat", 
                                 "steelblue4",
                                 "deepskyblue3",
                                 "azure4",
                                 "goldenrod2",
                                 "tan2",
                                 "darkcyan")) + NoLegend() + ggtitle(NULL)
pdf(file = "../2_Output/UMAP_Treatment.pdf", height = 4, width = 7)
p3
dev.off()
#Figure 2C -  Proportional Graph
TAA.combined$Group <- paste0(TAA.combined$Background, "_", TAA.combined$Treatment)
library(dittoSeq)
pdf(file = "../2_Output/Proportional.Bar_Treatment.pdf")
dittoBarPlot(
    object = TAA.combined,
    var = "CellType",
    group.by = "Group")+ ggtitle(NULL)
dev.off()
DittoPLOT <- dittoBarPlot(
    object = TAA.combined,
    var = "CellType",
    group.by = "SampleID",
    split.by = "Group",
    data.out = T)
write.csv(DittoPLOT$data, "../2_Output/Proportional.Cells.csv")
BarPlot <- dittoBarPlot(
    object = TAA.combined,
    var = "CellType",
    group.by = "SampleID",
    split.by = "Treatment", 
    data.out = T)
# Sample-specific proportional plot
pdf(file = "../2_Output/Proportional.Bar_Groups.pdf", height = 4, width = 6)
dittoBarPlot(
    object = TAA.combined,
    var = "CellType",
    group.by = "Background",
    split.by = "Treatment",
    retain.factor.levels = T,
    color.panel = c(Cardiomyocyte = "coral2", 
                    EC = "wheat", 
                    Fibroblast = "steelblue4", 
                    Mural_Cell = "deepskyblue3",
                    Lymphoid = "azure4",
                    Myeloid = "goldenrod2", 
                    Mast_Cell = "tan2",
                    Neural_Cell = "darkcyan")
    ) + 
  facet_wrap(~factor(Treatment, levels = c("Ctrl", "HFpEF")), 
             ncol = 3, 
             scales = "free_x") +
  ggtitle(NULL) + 
  labs(x = NULL)
dev.off()
# Extract Data for export
BarPlot_sampledata <- BarPlot$data %>% 
  select(-count, -label.count.total.per.facet) %>% 
  tidyr::pivot_wider(., names_from = "grouping", values_from = "percent")
openxlsx::write.xlsx(BarPlot_sampledata, "../2_Output/BarPlot_CellTypes.xlsx")
# Visualize the features/genes
VizDimLoadings(TAA.combined, dims = 1:5, reduction = "pca")
pdf(file = "../2_Output/PC_Heatmaps.pdf")
DimHeatmap(TAA.combined, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
```

# Loop: Cell Type-specific and Substrain-dependent Effects of HFD+L-NAME relative to Ctrl

To identify disease- and cell type-specific differentially-expressed
genes within the N vs J mice, iterative comparison across cell types was
performed. Differential gene expression analysis was carried out using
psuedobulk-based quantification of gene expression for each cell type
using the DESeq2 (1.44.0) algorithm within the R (4.4.1) statistical
computing environment.20,21 Gene-set enrichment analysis for all
comparisons was performed using the Elsevier curated pathway database
within EnrichR (3.2). Trajectory analysis of left ventricular
cardiomyocytes (“Cardiomyocyte”) populations was performed using
Monocle3 (1.4.22) in the R computing environment.


``` r
ann_colorInvestVec<-c(N_Ctrl="lightblue2", N_HFpEF = "steelblue3", J_Ctrl ="grey", J_HFpEF = "black")
ann_GROUP <- list(Group = c(N_Ctrl="lightblue2", N_HFpEF = "steelblue3", J_Ctrl ="grey", J_HFpEF = "black"), 
              Treatment = c("Ctrl" = "lightblue2", "HFpEF" = "steelblue3"),
              Background = c("J" = "coral2", "N" = "gray"))
names(ann_colorInvestVec)<-as.factor(c("N_Ctrl", "N_HFpEF", "J_Ctrl", "J_HFpEF"))
dbs <- c("GWAS_Catalog_2023") # Enrichment database
cellType_time <- Sys.time()
library(ggpubr)
library(Seurat)
library(dplyr)
library(DESeq2)
library(openxlsx)
TAA.combined <- readRDS(file = "../1_Input/NNT_Labelling_snRNA.rds")
TAA.combined$SampleID <- paste0("C", TAA.combined$SampleID)
cell_types <- as.character(unique(TAA.combined$CellType))
backgrounds <- as.character(unique(TAA.combined$Background))
# Forloop
for (cell_type in cell_types) {
  print(cell_type)
  tryCatch({
    for (background in backgrounds) {
      print(background)
    ## Create Folder Structure
    ifelse(!dir.exists(file.path(paste0("../2_Output/", cell_type, "/"))), 
           dir.create(file.path(paste0("../2_Output/", cell_type, "/"))), 
           FALSE)
    ifelse(!dir.exists(file.path(paste0("../2_Output/", cell_type, "/", background,"/"))), 
           dir.create(file.path(paste0("../2_Output/", cell_type, "/", background,"/"))), 
                      FALSE)
# --------------------------------------
# HFD vs. Ctrl Diet
TAA_CM_J <- subset(TAA.combined, subset = Treatment %in% c("Ctrl", "HFpEF") & CellType==cell_type & Background==background)
Idents(TAA_CM_J) <- "Treatment"
Index<-read.xlsx("../1_Input/M020_Sample.Info_NvJ_mep.xlsx", rowNames = T) %>% dplyr::filter(Background == background)
TAA_pseudo <- AggregateExpression(TAA_CM_J, return.seurat = F, assays = "RNA", group.by = "SampleID")$RNA
TAA_pseudo <- as.data.frame(TAA_pseudo)
colnames(TAA_pseudo) <- rownames(Index) #############
countData<-dplyr::select(TAA_pseudo, all_of(rownames(Index))) %>% dplyr::filter(rowSums(.)/ncol(.)>1)
# write.csv(countData, paste0("../1_Input/", "Normalized.Counts_",cell_type, ".csv")) # need to run once because this is used downstream for heatmaps, etc
######### RUN DESeq2
dds<-DESeq2::DESeqDataSetFromMatrix(countData=countData, colData = Index, design= ~Treatment)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds<-DESeq(dds, test="Wald", fitType="parametric")
resdf<-as.data.frame(results(dds, format = "DataFrame")) %>% merge(., countData, by = 0)
resdf$external_gene_name<-resdf$Row.names
rownames(resdf) <- resdf$Row.names
# Write excel
openxlsx::write.xlsx(resdf, file = paste0("../2_Output/", cell_type, "/", background, "/", background, "_", cell_type, "_HFpEF.vs.Ctrl_DEGs.xlsx"), rowNames=T)
de_CM <- openxlsx::read.xlsx(paste0("../2_Output/", cell_type, "/", background, "/", background, "_", cell_type, "_HFpEF.vs.Ctrl_DEGs.xlsx"), rowNames=T) 
# Volcano Plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
results = mutate(de_CM, minuslogpvalue = -log(pvalue), log2FC=log2FoldChange)
results<-results %>% filter(pvalue!=0)
results$gene_name<-rownames(results)
results <- results %>% 
  mutate(., sig=ifelse(pvalue<0.0001 & log2FC>.3, 
                       "P < 0.0001 and Log(Fold-Change) > 0.3", 
                       ifelse(pvalue<0.0001 & log2FC< 0-.3,
                              "P < 0.0001 and Log(Fold-Change) < -0.3", 
                              "Not Sig")
                       )
         )
results$sig<-factor(results$sig, 
levels = c("P < 0.0001 and Log(Fold-Change) < -0.3",
  "Not Sig",
  "P < 0.0001 and Log(Fold-Change) > 0.3")
  )
max(results$minuslogpvalue, na.rm = TRUE)
max(results$log2FC, na.rm = TRUE)
min(results$log2FC, na.rm = TRUE)
p = ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  geom_point(aes(fill=sig, size = minuslogpvalue),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=0.3, size=.5, linetype="dashed") +
  geom_vline(xintercept=-0.3, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0-log(0.0001), size=.5, linetype="dashed") +
  labs(x=expression(Log[2](Fold-Change)), y=expression(-Log[10](P-value))) + 
  xlim(min(results$log2FC, na.rm = TRUE),max(results$log2FC, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(0, max(results$minuslogpvalue, na.rm = TRUE)), expand = c(0,0)) +
  # geom_hline(yintercept = 0, size = 1) + 
  # geom_vline(xintercept=0, size=0.5) +
  scale_fill_manual(values=c("dodgerblue2", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-0.3), 15, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>0.3), 15, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
pdf(file = paste0("../2_Output/", cell_type, "/", background, "/", background, "_Volcano_HFDLNAMEvCon.pdf"), height = 5.5, width = 5)
  print(
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-0.3), 10, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>0.3), 10, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
  )
dev.off()
## Figure 4B - GSEA
DEGs_UP <- results %>% filter(log2FoldChange > 0 & pvalue<0.05) #%>% top_n(n = 250, -padj)
DEGs_DOWN <- results %>% filter(log2FoldChange < 0 & pvalue<0.05) #%>% top_n(n = 250, -padj)
##Enrichr
library(enrichR)
library(stringr)
enriched_UP <- enrichr(rownames(DEGs_UP), dbs)
enrich_UP<-enriched_UP[[dbs]]
head(enrich_UP)
enriched_DOWN <- enrichr(rownames(DEGs_DOWN), dbs)
enrich_DOWN<-enriched_DOWN[[dbs]]
head(enrich_DOWN)
write.csv(enrich_DOWN, paste0("../2_Output/", cell_type, "/", background, "/", background, "_HFD_PathwayDEGs_DOWN.csv"))
write.csv(enrich_UP, paste0("../2_Output/", cell_type, "/", background, "/", background, "_PathwayDEGs_UP.csv"))
#### Heatmap
    MDS_data<-read.csv(paste0("../1_Input/Normalized.Counts_", cell_type, ".csv"), row.names = 1)
    names(MDS_data) <- sub('^X', 'C', names(MDS_data))
    Index<-read.xlsx("../1_Input/M020_Sample.Info_NvJ_mep.xlsx") %>% mutate(Group = paste0(Background,"_",Treatment))
    rownames(Index) <- Index$MouseID
    Index$MouseID <- as.character(Index$MouseID)
    #
    Index_hm <- Index %>% dplyr::select(Treatment, Background)
    results_p05<-results %>% filter(pvalue<0.05)
    HM_data <- subset(MDS_data, rownames(MDS_data) %in% results_p05$Row.names) %>% data.matrix()
    paletteLength <- 100
    Index_hm$Treatment <- factor(Index_hm$Treatment, levels = c("Ctrl", "HFpEF"))
    Index_hm$Background <- factor(Index_hm$Background, levels = c("N", "J"))
    hm_COLOR = list(Background = c(N="gray", J = "coral2"), Treatment = c(Ctrl ="lightblue2", HFpEF = "steelblue"))
myColor <- colorRampPalette(c("dodgerblue4", "white", "gold2"))(paletteLength)
pheatmap::pheatmap(HM_data, scale="row",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         cutree_cols = 2,
         cutree_rows = 2,
         angle_col = 45,
         fontsize_col = 8,
         color = myColor,
         annotation_colors = hm_COLOR,
         show_rownames = FALSE,
         border_color = NA,
         annotation_col = Index_hm,
         filename=paste0("../2_Output/", cell_type,  "/", background, "/", cell_type, "_", background,"_Heatmap_Normcount.pdf"))
    }
###### N vs. J - HFpEF Only (for each cell-type)
dbs <- c("GWAS_Catalog_2023") # Enrichment database
TAA_CM_HFpEF <- subset(TAA.combined, subset = Treatment == "HFpEF" & CellType==cell_type)
Idents(TAA_CM_HFpEF) <- "Background"
Index_HFpEF<-read.xlsx("../1_Input/M020_Sample.Info_NvJ_mep.xlsx", rowNames = T) %>% dplyr::filter(Treatment == "HFpEF")
TAA_pseudo_HFpEF <- AggregateExpression(TAA_CM_HFpEF, return.seurat = F, assays = "RNA", group.by = "SampleID")$RNA
TAA_pseudo_HFpEF <- as.data.frame(TAA_pseudo_HFpEF)
names(TAA_pseudo_HFpEF) <- sub('^g', '', names(TAA_pseudo_HFpEF))
countData_HFpEF<-dplyr::select(TAA_pseudo_HFpEF, all_of(rownames(Index_HFpEF))) %>% dplyr::filter(rowSums(.)/ncol(.)>1)
dds_HFpEF <- DESeq2::DESeqDataSetFromMatrix(countData=countData_HFpEF, colData = Index_HFpEF, design= ~Background)
dds_HFpEF <- dds_HFpEF[ rowSums(counts(dds_HFpEF)) > 1, ]
dds_HFpEF<-DESeq(dds_HFpEF, test="Wald", fitType="parametric")
resdf_HFpEF<-as.data.frame(results(dds_HFpEF, format = "DataFrame")) %>% merge(., countData_HFpEF, by = 0)
resdf_HFpEF$external_gene_name<-resdf_HFpEF$Row.names
write.xlsx(resdf_HFpEF, file = paste0("HFpEF_",cell_type,"_DEGs.xlsx"))
# resdf_HFpEF <- read.xlsx(paste0("HFpEF_",cell_type,"_DEGs.xlsx"), rowNames = T)
# dbs <- c("GWAS_Catalog_2023") # Enrichment database
resdf_HFpEF <- resdf_HFpEF %>% filter(padj < 0.05) #%>% top_n(n = 500, -pvalue)
enriched_HFpEF <- enrichr(toupper(resdf_HFpEF$external_gene_name), dbs)
enrich_HFpEF<-enriched_HFpEF[[dbs]] %>% arrange(P.value)
head(enrich_HFpEF)
write.csv(enrich_HFpEF, paste0("HFpEF_", cell_type,"_GWAS.csv"))

###### N vs. J - Ctrl Only (for each cell-type)
TAA_CM_Ctrl <- subset(TAA.combined, subset = Treatment == "Ctrl" & CellType==cell_type)
Idents(TAA_CM_Ctrl) <- "Background"
Index_Ctrl<-read.xlsx("../1_Input/M020_Sample.Info_NvJ_mep.xlsx", rowNames = T) %>% dplyr::filter(Treatment == "Ctrl")
TAA_pseudo_Ctrl <- AggregateExpression(TAA_CM_Ctrl, return.seurat = F, assays = "RNA", group.by = "SampleID")$RNA
TAA_pseudo_Ctrl <- as.data.frame(TAA_pseudo_Ctrl)
names(TAA_pseudo_Ctrl) <- sub('^g', '', names(TAA_pseudo_Ctrl))
countData_Ctrl<-dplyr::select(TAA_pseudo_Ctrl, all_of(rownames(Index_Ctrl))) %>% dplyr::filter(rowSums(.)/ncol(.)>1)
dds_Ctrl <- DESeq2::DESeqDataSetFromMatrix(countData=countData_Ctrl, colData = Index_Ctrl, design= ~Background)
dds_Ctrl <- dds_Ctrl[ rowSums(counts(dds_Ctrl)) > 1, ]
dds_Ctrl<-DESeq(dds_Ctrl, test="Wald", fitType="parametric")
resdf_Ctrl<-as.data.frame(results(dds_Ctrl, format = "DataFrame")) %>% merge(., countData_Ctrl, by = 0)
resdf_Ctrl$external_gene_name<-resdf_Ctrl$Row.names
write.xlsx(resdf_Ctrl, file = paste0("Ctrl_", cell_type,"_DEGs.xlsx"))
# resdf_Ctrl <- read.xlsx(paste0("Ctrl_", cell_type, "_DEGs.xlsx"), rowNames = T)
dbs <- c("GWAS_Catalog_2023") # Enrichment database
resdf_Ctrl <- resdf_Ctrl %>% filter(padj < 0.05) #%>% top_n(n = 500, -pvalue)
enriched_Ctrl <- enrichr(toupper(resdf_Ctrl$external_gene_name), dbs)
enrich_Ctrl<-enriched_Ctrl[[dbs]] #%>% arrange(P.value)
head(enrich_Ctrl)
write.csv(enrich_Ctrl, paste0("Ctrl_", cell_type, "_GWAS.csv"))
########### VENN DIAGRAM
library(ggVennDiagram)
x<-list(Ctrl = resdf_Ctrl$external_gene_name, HFpEF = resdf_HFpEF$external_gene_name)
library(VennDiagram)
venn.diagram(x, 
             fill = c("darkcyan", "grey"),
             alpha = c(0.75, 0.75),
             lty = 'blank',
             main.cex = 1.5,
             sub.cex = 2,
             main = cell_type,
             filename = paste0("../2_Output/", cell_type, "/", cell_type, "_VENN_N.v.J.svg"),
             imagetype = "svg",
             na = "remove",
             disable.logging = T,
             width = 5,
             height = 5,
             units = "in")

##################################################################
    # Unsupervised Analysis
library(limma)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
#Filter normalized counts to remove outliers
    MDS_data<-read.csv(paste0("../1_Input/Normalized.Counts_", cell_type, ".csv"), row.names = 1)
    names(MDS_data) <- sub('^X', 'C', names(MDS_data))
    Index<-read.xlsx("../1_Input/M020_Sample.Info_NvJ_mep.xlsx") %>% mutate(Group = paste0(Background,"_",Treatment))
    rownames(Index) <- Index$MouseID
    Index$MouseID <- as.character(Index$MouseID)
    #
    Index_hm <- Index %>% dplyr::select(Treatment, Background)
MDS_data<-dplyr::select(MDS_data, all_of(rownames(Index)))
# MDS in ggplot2
Ntop = 1000
library(magrittr)
library(dplyr)
library(ggpubr)
library(matrixStats)
library("ggrepel")
library(wesanderson)
MDS.set<-as.data.frame(MDS_data)
RowVar<-rowVars(data.matrix(MDS.set)) #calculate variances for each row (vector)
MDS.set<-as.data.frame(cbind(MDS.set, RowVar)) #Add to the MDS.set dataset
MDS_matrix<-MDS.set %>% arrange(desc(RowVar)) %>% top_n(Ntop,RowVar) #Select top N rows by variance
# Compute MDS
mds <- MDS_matrix %>% dplyr::select(-RowVar) %>% t(.) %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
rownames(mds)<-rownames(Index)
mds$MouseID<-rownames(mds)
mds<-dplyr::inner_join(mds, Index)
#K-means clustering#K-means Indexclustering
clust <- kmeans(mds[,1:2], 2)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(kmeans.2 = clust)
###
library(ggpubr)
library(cowplot) 
# Main plot
pmain <- ggplot(mds, aes(x = Dim.1, y = Dim.2, color = Group))+
  scale_color_manual(values = ann_colorInvestVec) +
  theme_classic()+
  theme(panel.background = element_rect("white", colour = "black", size=2), 
      # panel.grid.major = element_line(colour = "gray75", size=.75), 
      # panel.grid.minor = element_line(colour = "gray90", size=0.4),
      legend.position="bottom",
      legend.key=element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face="bold")) +
  geom_hline(yintercept = 0, size = .5, linetype = 2) +
  geom_vline(xintercept=0, size=.5, linetype = 2) +
  geom_point(aes(size = 4))+ #Add points for each sample
  # stat_ellipse()+ # create elliptical shapes
  # geom_text_repel(data=mds, aes(label=Sample_ID), show.legend  = F) + #label the samples
  labs(x="Principal Component 1", 
       y="Principal Component 2")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = mds, aes(x = Dim.1, fill = Group),
              alpha = 0.7, size = 0.5)+
  scale_fill_manual(values = ann_colorInvestVec)
# Marginal densities along y axis
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+ #must set coord_flip = true if using coord_flip() below
  geom_density(data = mds, aes(x = Dim.2, fill = Group),
                alpha = 0.7, size = 0.5)+
  scale_fill_manual(values = ann_colorInvestVec)+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
pdf(file=paste0("../2_Output/", cell_type, "/", cell_type, "_MDS.pdf"), height = 5, width = 5, onefile = F)
print(ggdraw(p2))
dev.off()
######################### PCA Analysis
#Plot Features of the PCA
library(dplyr)
library(plotly)
##Import the data to be used for PCA
rownames(Index)<-Index$MouseID
#transpose the dataset (required for PCA)
data.pca<-t(MDS_data)
data.pca<-as.data.frame(data.pca)
##merge the file
data.pca_Final<-merge(Index, data.pca, by=0)
rownames(data.pca_Final)<-data.pca_Final$Row.names
pca.comp<-prcomp(data.pca_Final[,(ncol(Index)+2):ncol(data.pca_Final)])

pcaCharts=function(x) {
    x.var <- x$sdev ^ 2
    x.pvar <- x.var/sum(x.var)
    par(mfrow=c(2,2))
    plot(x.pvar,xlab="Principal component",
         ylab="Proportion of variance", ylim=c(0,1), type='b')
    plot(cumsum(x.pvar),xlab="Principal component",
         ylab="Cumulative Proportion of variance",
         ylim=c(0,1),
         type='b')
    screeplot(x)
    screeplot(x,type="l")
    par(mfrow=c(1,1))
                  }
pcaCharts(pca.comp)
png(file=paste0("../2_Output/", cell_type, "/", background, "/", cell_type, "_" , background, "_PCA.Charts.png"))
pcaCharts(pca.comp)
dev.off()
######
library(dplyr)
library(pathview)
library(biomaRt)
library(openxlsx)
CM_N <- openxlsx::read.xlsx(paste0("../2_Output/", cell_type, "/", backgrounds[1], "/", backgrounds[1], "_", cell_type,"_HFD.v.CON_DEGs.xlsx"), rowNames=T)
CM_N$GeneName <- rownames(CM_N)
DE_CM_N <- CM_N %>% filter(pvalue<0.05)
CM_J <- openxlsx::read.xlsx(paste0("../2_Output/", cell_type, "/", backgrounds[2], "/", backgrounds[2], "_", cell_type,"_HFD.v.CON_DEGs.xlsx"), rowNames=T)
CM_J$GeneName <- rownames(CM_J)
DE_CM_J <- CM_J %>% filter(pvalue<0.05)
# J Only DEGs
N_UP<-dplyr::filter(DE_CM_N, log2FoldChange>0)
N_DOWN<-filter(DE_CM_N, log2FoldChange<0)
N_ONLY<-anti_join(DE_CM_N, DE_CM_J, by = "GeneName")
N_ONLY.UP<-N_ONLY %>% filter(log2FoldChange>0)
N_ONLY.DOWN<-N_ONLY %>% filter(log2FoldChange<0)
write.xlsx(N_ONLY, paste0("../2_Output/", cell_type, "/", cell_type, "_N_ONLY.xlsx"), overwrite = TRUE)
# DKO Only DEGs
J_UP<-filter(DE_CM_J, log2FoldChange>0)
J_DOWN<-filter(DE_CM_J, log2FoldChange<0)
J_ONLY<-anti_join(DE_CM_J, DE_CM_N, by = "GeneName")
J_ONLY.UP<-J_ONLY %>% filter(log2FoldChange>0)
J_ONLY.DOWN<-J_ONLY %>% filter(log2FoldChange<0)
write.xlsx(J_ONLY, paste0("../2_Output/", cell_type, "/", cell_type, "_J.ONLY.xlsx"), overwrite = TRUE)
# Overlapping DEGs
Conserved_DEGs<-inner_join(DE_CM_N, DE_CM_J, by = "GeneName") 
rownames(Conserved_DEGs)<-make.unique(Conserved_DEGs$GeneName, sep = ".")
Conserved_DEGs <- Conserved_DEGs %>% rename_all(~stringr::str_replace_all(.,c("\\.y"="_J", "\\.x"="_N")))
Conserved_Both.UP<-Conserved_DEGs %>% filter(log2FoldChange_N > 0, log2FoldChange_J > 0)
Conserved_Both.DOWN<-Conserved_DEGs %>% filter(log2FoldChange_N < 0, log2FoldChange_J < 0)
Conserved_Inverse<-Conserved_DEGs %>% filter((log2FoldChange_N>0 & log2FoldChange_J<0) | (log2FoldChange_N<0 & log2FoldChange_J>0))
write.xlsx(Conserved_DEGs, paste0("../2_Output/", cell_type, "/", cell_type, "_Overlapping_DEGs.xlsx"), overwrite = TRUE)
#Merge dataframe for IPA
Merged<-full_join(DE_CM_N, DE_CM_J, by = "GeneName")
rownames(Merged)<-make.unique(Merged$GeneName, sep = ".")
Merged <- Merged %>% rename_all(~stringr::str_replace_all(.,c("\\.y"="_J", "\\.x"="_N")))
write.xlsx(Merged, paste0("../2_Output/", cell_type, "/", cell_type, "_Merged_DEGs.xlsx"), overwrite = TRUE)
########### VENN DIAGRAM
library(ggVennDiagram)
x<-list(N = DE_CM_N$GeneName, J = DE_CM_J$GeneName)
library(VennDiagram)
venn.diagram(x, fill = c("red", "grey"), alpha = c(0.75, 0.75), lty = 'blank', main.cex = 3, sub.cex = 2, sub = "HFD+L-NAME vs. Ctrl", main = cell_type,filename = paste0("../2_Output/", cell_type, "/", cell_type, "_DEGs.Overlap.svg"), imagetype = "svg", na = "remove", disable.logging = T, width = 8, height = 8, units = "in")
#Write excel worksheet
wb_DESeq<-createWorkbook()
#Unfiltered
  addWorksheet(wb_DESeq, "J_ONLY_p05")
  writeData(wb_DESeq, "J_ONLY_p05", J_ONLY, startCol = 1)
#P-value Significant (0.05)
  addWorksheet(wb_DESeq, "N_ONLY_p05")
  writeData(wb_DESeq, "N_ONLY_p05", N_ONLY, startCol = 1)
#Q-value Significant (0.05)
  addWorksheet(wb_DESeq, "CM_Conserved_DEGs")
  writeData(wb_DESeq, "CM_Conserved_DEGs", Conserved_DEGs, startCol = 1)
saveWorkbook(wb_DESeq, file = paste0("../2_Output/", cell_type, "/", cell_type, "_Venn.Diagram.xlsx"), overwrite = TRUE)
############################################
# Volcano Plot - Combined
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
ALL_DEGs <- full_join(DE_CM_N, DE_CM_J, by = "GeneName")
rownames(ALL_DEGs) <- ALL_DEGs$GeneName
ALL_DEGs <- ALL_DEGs %>% rename_all(~stringr::str_replace_all(.,c("\\.y"="_J", "\\.x"="_N")))
results <- ALL_DEGs %>% 
  mutate(., baseMean_mean=rowMeans(dplyr::select(., starts_with("baseMean"))), sig=ifelse(log2FoldChange_N>0 & log2FoldChange_J>0, 
                       "log2FoldChange_N > 0 log2FoldChange_J > 0", 
                       ifelse(log2FoldChange_N<0 & log2FoldChange_J<0,
                              "log2FoldChange_N < 0 and log2FoldChange_J < 0", 
                              "Not Sig")
                       )
         ) %>%
#  mutate(., avg_pct = rowMeans(dplyr::select(.,pct.1_N, pct.2_N, pct.1_J, pct.2_J))) %>%
  mutate(., Abs_FC = abs(log2FoldChange_J)+abs(log2FoldChange_N))
results$sig<-factor(results$sig, 
levels = c("log2FoldChange_N < 0 and log2FoldChange_J < 0",
  "Not Sig",
  "log2FoldChange_N > 0 log2FoldChange_J > 0")
  )
p = ggplot(results, aes(log2FoldChange_N, log2FoldChange_J)) + 
  theme_classic() +
  geom_point(aes(fill=sig, size = baseMean_mean),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=0, size=.5, linetype="dashed") +
  # geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0, size=.5, linetype="dashed") +
  # geom_hline(yintercept=-1, size=.5, linetype="dashed") +
  labs(x="HFD+L-NAME vs. Ctrl (N)", y="HFD+L-NAME vs. Ctrl (J)") + 
  xlim(min(results$log2FoldChange_N, na.rm = TRUE),max(results$log2FoldChange_N, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(min(results$log2FoldChange_J, na.rm = TRUE), max(results$log2FoldChange_J, na.rm = TRUE)), expand = c(0,0)) +
  scale_fill_manual(values=c("dodgerblue2", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, log2FoldChange_J< -.5 & log2FoldChange_N < -.5), 15, Abs_FC),
                  aes(label=GeneName)) +
  geom_text_repel(data=top_n(filter(results, log2FoldChange_J>.5 & log2FoldChange_N > .5), 15, Abs_FC), 
  aes(label=GeneName)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
pdf(file = paste0("../2_Output/", cell_type, "/", cell_type, "_Combined.pdf"), height = 5, width = 5)
  print(
    p+
  geom_text_repel(data=top_n(filter(results, log2FoldChange_J< -.5 & log2FoldChange_N < -.5), 15, Abs_FC),
                  aes(label=GeneName)) +
  geom_text_repel(data=top_n(filter(results, log2FoldChange_J>.5 & log2FoldChange_N > .5), 15, Abs_FC), 
  aes(label=GeneName)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
  )
dev.off()
###################
DEGs_UP <- results %>% filter(log2FoldChange_N > 0, log2FoldChange_J > 0, padj_N < 0.05, padj_J < 0.05)
DEGs_DOWN <- results %>% filter(log2FoldChange_N < -0, log2FoldChange_J < -0, padj_N < 0.05, padj_J < 0.05)
# Volcano Plot - N Only
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
ALL_DEGs <- anti_join(CM_N, CM_J, by = "GeneName")
rownames(ALL_DEGs) <- ALL_DEGs$GeneName
results <- ALL_DEGs %>% 
 mutate(., baseMean_mean=rowMeans(dplyr::select(., starts_with("baseMean"))), sig=ifelse(log2FoldChange>1, 
                       "log2FoldChange > 1", 
                       ifelse(log2FoldChange<0-1,
                              "log2FoldChange < -1", 
                              "Not Sig")
                       )
         ) %>%
  mutate(., Abs_FC = abs(log2FoldChange), minuslogpvalue=0-log(pvalue, 10))
results$sig<-factor(results$sig, levels = c("log2FoldChange < -1","Not Sig", "log2FoldChange > 1"))
p = ggplot(results, aes(log2FoldChange, minuslogpvalue)) + 
  theme_classic() +
  geom_point(aes(fill=sig, size = Abs_FC),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=0, size=.5, linetype="dashed") +
  # geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0, size=.5, linetype="dashed") +
  # geom_hline(yintercept=-1, size=.5, linetype="dashed") +
  labs(x="Log(Fold-Change)", y="-Log(P-value)") + 
  xlim(min(results$log2FoldChange, na.rm = TRUE),max(results$log2FoldChange, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(0, max(results$minuslogpvalue, na.rm = TRUE)), expand = c(0,0)) +
  scale_fill_manual(values=c("dodgerblue2", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, log2FoldChange< -1), 15, Abs_FC),
                  aes(label=GeneName)) +
  geom_text_repel(data=top_n(filter(results, log2FoldChange > 1), 15, Abs_FC), 
  aes(label=GeneName)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
pdf(file = paste0("../2_Output/", cell_type, "/", cell_type, "_Volcano_N.Only.pdf"), height = 7, width = 5)
  print(
  p+
  geom_text_repel(data=top_n(filter(results, log2FoldChange< -1), 15, Abs_FC),
                  aes(label=GeneName)) +
  geom_text_repel(data=top_n(filter(results, log2FoldChange > 1), 15, Abs_FC), 
  aes(label=GeneName)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
  )
dev.off()

# Candidate gene analysis
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(gtools)
library(openxlsx)
# Parameters
# GOI <- c("Slc25a51", "Naxe", "Nnt", "Sirt5", "Nudt13", "Sirt3", "Nadk2", "Nmnat3", "Naxd", "Kyat3", "Kmo", "Sirt4", "Noct")
ifelse(!dir.exists(file.path(paste0("../2_Output/", cell_type, "/Pathway"))), dir.create(file.path(paste0("../2_Output/", cell_type, "/Pathway"))), FALSE)
# choose OXPHOS genes
library("KEGGREST")
names <- keggGet("mmu00190")[[1]]$GENE
namesodd <-  names[seq(0,length(names),2)]
GOI <- gsub("\\;.*","",namesodd)
write.csv(GOI, file = "mmu01212.csv",quote = F, row.names = F)
GOI <- c(intersect(GOI, N_ONLY$GeneName), "Idh1", "Mdh1")
#Import Index file
Counts<-read.csv(paste0("../1_Input/Normalized.Counts_", cell_type, ".csv"), row.names = 1)
names(Counts) <- sub('^X', 'C', names(Counts))
colData<-read.xlsx("../1_Input/M020_Sample.Info_NvJ_mep.xlsx") %>% mutate(Group = paste0(Background,"_",Treatment))
rownames(colData) <- colData$MouseID
colData$MouseID <- as.character(colData$MouseID)
colData <- colData %>% mutate(Group = paste0(Background, "_", Treatment))
#
colData <- colData %>% dplyr::select(MouseID,Treatment, Background, Group)
Counts<-dplyr::select(Counts, all_of(rownames(colData)))

#Filter results by the gene vector
DEGs<-subset(Counts, rownames(Counts) %in% GOI)
tDEGs<-as.data.frame(t(DEGs))
## convert all genes to numeric (from factors)
asNumeric=function(x){as.numeric(as.character(x))}
factorsNumeric=function(d){modifyList(d, lapply(d[, sapply(d, is.character)], asNumeric))}
##
tDEGs<-factorsNumeric(tDEGs)
tDEGs$MouseID<-rownames(tDEGs)
colData.ex<-dplyr::inner_join(tDEGs, colData)
colData.ex$Group<-factor(colData.ex$Group, levels = c("J_Ctrl", "J_HFpEF", "N_Ctrl", "N_HFpEF"))
colData.ex<-dplyr::group_by(colData.ex, "Group") #define groups for statistics
write.xlsx(colData.ex, paste0("../2_Output/", cell_type, "/Pathway/Candidate_genes.xlsx"), overwrite = T)
## For loop creating a graph for each gene
plotlist = list()
p<-1
for (i in GOI){
  tryCatch({ # skip iteration if error results (i.e. gene not found)
g_plot<-ggboxplot(colData.ex, x = "Group", 
          y = i, 
          fill = "Group",
          add = "jitter"
          ) + 
  scale_fill_manual(values = ann_colorInvestVec) +
  # stat_compare_means(aes(group = Group),
  #                   comparisons = my_comparisons,
  #                   label = "p.signif",
  #                   bracket.nudge.y = 5
  #                   ) +
  theme(axis.text.x=element_text(size=rel(0.75), 
                                 angle = 45, hjust = 1), 
        axis.text.y=element_text(size=rel(0.75)), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold"), 
        legend.position="none") + # resize labels, remove legend
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))  # expand = expansion(mult = c(0, 0.1)) ### Y scale (to see the statistics)
pdf(file=paste0("../2_Output/", cell_type, "/Pathway/", i, "_Expression.pdf"), width = 6, height = 6)
print(g_plot)
dev.off()
plotlist[[i]] = g_plot
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # skip iteration if error results (i.e. gene not found)
}
t<-marrangeGrob(grobs = plotlist, legend, nrow=2, ncol=2)
ggsave(paste0("../2_Output/", cell_type, "/Pathway/_DEGs_Pathway.pdf"), t, width = 6, height = 7)
dev.off()
t
#Heatmap
rownames(colData.ex)<-colData.ex$MouseID
colnames(colData.ex)<-as.character(colnames(colData.ex))
hm_data.candidate<-colData.ex %>% dplyr::select(-MouseID, -Treatment, -Background, -Group)
hm_data.t<-t(hm_data.candidate[,-ncol(hm_data.candidate)])
colnames(hm_data.t)<-colData.ex$MouseID
hm_data.t<-data.matrix(hm_data.t)
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "gold2"))(paletteLength)
hm_COLOR = list(Group = c(N_Ctrl="black", N_HFpEF = "darkcyan", J_CON ="grey", J_HFD.LNAME = "coral2"))
Index_hm <- Index %>% dplyr::select(-MouseID)
pheatmap::pheatmap(hm_data.t, scale="row",
         cluster_cols = F,
         cluster_rows = T,
         angle_col = 45,
         fontsize_col = 8,
         color = myColor,
         show_rownames = T,
         border_color = NA,
         annotation_colors = ann_GROUP,
         annotation_col = Index, 
        filename = paste0("../2_Output/", cell_type, "/Pathway/_Heatmap.pdf"))
dev.off()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
cellType_endtime <- Sys.time()
cellType_endtime - cellType_time


Pathways <- data.frame()
for (cell_type in cell_types){
  dframe <- read.csv(paste0("../2_Output/", cell_type, "/", "Pathway.Enrich_", cell_type, "_Ctrl_N.v.J_GWAS.csv"), row.names = 1) %>% dplyr::select(-Old.P.value, -Old.Adjusted.P.value, -Odds.Ratio)
  if (nrow(dframe)>0){
  dframe$cell_type <- cell_type
  Pathways <- bind_rows(Pathways, dframe)
  }
  print(Pathways)
}
Pathways_CVD <- Pathways %>% filter(grepl("Blood Pressure|Coronary|Atrial|Heart|Card|Stroke|Cholesterol|Obesity", Term)) %>% filter(P.value<0.05)
write.csv(Pathways_CVD, "../2_Output/CVD.Pathway_GWAS.csv")
```

# Cardiomyocyte-Specific Cell State Transition

Trajectory analysis was employed to characterize the phenotypic shifts
between VSMC_A and VSMCB among non-diseased aortic control tissue to
better understand how VSMCs might transition from one functional state
to another over time or in response to specific pathologic contexts.


``` r
# Check if the folder exists, and create it if it doesn't
if (!dir.exists("../2_Output/Cardiomyocyte/Trajectory")) {dir.create("../2_Output/Cardiomyocyte/Trajectory")
  cat("Folder created:", "../2_Output/Cardiomyocyte/Trajectory", "\n")
} else {
  cat("Folder already exists:", "../2_Output/Cardiomyocyte/Trajectory", "\n")
}
#
ann_colors = list(Treatment = c(Ctrl = "black", HFpEF ="darkcyan"),
                  Background = c("J" = "coral2", "N" = "gray"),
                  CellType = c(Cardiomyocyte = "coral2", 
                    EC = "wheat", 
                    Fibroblast = "steelblue4", 
                    Mural_Cell = "deepskyblue3",
                    Lymphoid = "azure4",
                    Myeloid = "goldenrod2", 
                    Mast_Cell = "tan2",
                    Neural_Cell = "darkcyan"))
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggtrace)
Merged_scaled <- readRDS(file = "../1_Input/NNT_Labelling_snRNA.rds")
Merged_scaled@active.assay = "RNA"
Merged_scaled$Treatment <- factor(Merged_scaled$Treatment, levels = c("Ctrl", "HFpEF"))
Merged_scaled <- subset(Merged_scaled, subset = Treatment == c("HFpEF"))
Merged_scaled <- subset(Merged_scaled, subset = Background == "N") 
Merged_scaled <- subset(Merged_scaled, subset = CellType == "Cardiomyocyte")
Merged_scaled <- PrepSCTFindMarkers(Merged_scaled, assay = "SCT", verbose = TRUE)
# Re-run the clustering
hdat_Cardiomyocyte <- Merged_scaled
DefaultAssay(hdat_Cardiomyocyte ) <- "SCT"
hdat_Cardiomyocyte <- FindClusters(hdat_Cardiomyocyte, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0), graph.name = "SCT_snn")
hdat_Cardiomyocyte <- FindNeighbors(object = hdat_Cardiomyocyte, reduction = "harmony")
hdat_Cardiomyocyte <- RunUMAP(Merged_scaled, assay = "SCT", dims = 1:40)
# hdat_Cardiomyocyte <- RunPCA(Merged_scaled, assay = "SCT")
# hdat_Cardiomyocyte <- RunTSNE(Merged_scaled, assay = "SCT")
hdat_Cardiomyocyte<-SetIdent(hdat_Cardiomyocyte, value = "SCT_snn_res.0.2")
UMAP_Cardiomyocyte<-DimPlot(hdat_Cardiomyocyte,  label = T, reduction = "umap") + NoLegend()
# PCA_Cardiomyocyte <-DimPlot(hdat_Cardiomyocyte,  label = T, reduction = "pca") + NoLegend()
# TSNE_Cardiomyocyte <-DimPlot(hdat_Cardiomyocyte,  label = T, reduction = "tsne") + NoLegend() 
# TSNE_Cardiomyocyte
# PCA_Cardiomyocyte
UMAP_Cardiomyocyte
# convert to singlecell_dataset for downstream analysis
cds<-SeuratWrappers::as.cell_data_set(hdat_Cardiomyocyte)
cds <- estimate_size_factors(cds)
# cds <- preprocess_cds(cds, method = "PCA", num_dim = 50)
# Include gene names (not done by default by the seurat conversion)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
# # Run Monocle
# cds <- reduce_dimension(cds, 
#                         reduction_method = "UMAP", 
#                         preprocess_method = "PCA")
cds <- cluster_cells(cds, 
                     reduction_method = "UMAP",
                     resolution = 0.000002) #This step creates "partitions" that are used in the trajectory inference
plot_cells(cds, 
           reduction_method = "UMAP",
           show_trajectory_graph = FALSE, 
           color_cells_by = "partition") # this shows the partitions overlain on the UMAP
cds <- learn_graph(cds, use_partition = TRUE) # creating trajectory inference within each partition (if `use_partition = TRUE`)
# cds <- order_cells(cds) # Used when setting the nodes; if already known, use the next line
# names(which(cds@principal_graph_aux$UMAP$pseudotime==0))
root_cell <- "TTGAGTGAGAATTTGG-1_4" # names(which(cds@principal_graph_aux$UMAP$pseudotime==0))
cds<-order_cells(cds, root_cells = root_cell)
# Plot the pseudotime on UMAP
plot_cells(cds, 
           color_cells_by = "pseudotime",
           label_branch_points = FALSE,
           label_leaves = FALSE)
pdf(file = "../2_Output/Cardiomyocyte/Trajectory/Trajectory_Cardiomyocyte_UMAP_Partition.pdf", height = 4, width = 4)
plot_cells(cds,
           color_cells_by = "partition",
           show_trajectory_graph = F,
           graph_label_size = 3,
           cell_size = .5,
           label_roots = F,
           label_branch_points = FALSE,
           label_leaves = FALSE)
dev.off()
pdf(file = "../2_Output/Cardiomyocyte/Trajectory/Trajectory_Cardiomyocyte_UMAP_Pseudotime.pdf", height = 4, width = 5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = F,
           graph_label_size = 1,
           cell_size = .7,
           label_branch_points = FALSE,
           label_leaves = F)
dev.off()
# Identify pseudotime
modulated_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8) # Identify differentially-expressed genes with pseudotime
modulated_genes <- na.omit(modulated_genes) # remove NA's
modulated_genes <- modulated_genes %>% filter(modulated_genes$q_value < 0.05 & modulated_genes$status =="OK") # filter cds results down
modulated_genes <- modulated_genes[order(-modulated_genes$morans_test_statistic), ] # order by moran's test
modulated_genes <- top_n(modulated_genes, 1000, -q_value)
#### Create a heatmap of genes with similar pseudotime kinetics
genes <- row.names(subset(modulated_genes, q_value < 0.05))
openxlsx::write.xlsx(modulated_genes, "../2_Output/Cardiomyocyte/Trajectory/Pseudotime_DEGs.xlsx")
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
pt.matrix <- as.data.frame(exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))])
#
cell_names <- colnames(pt.matrix)
Index<-as.data.frame(cds@colData) %>% dplyr::select(CellType, Background, Treatment)
Index<-subset(Index, row.names(Index) %in% cell_names)
Index$CellType <- factor(Index$CellType, levels = c("Cardiomyocyte", 
                                                    "EC", 
                                                    "Fibroblast",
                                                    "Mural_Cell",
                                                    "Lymphoid",
                                                    "Myeloid",
                                                    "Mast_Cell",
                                                    "Neural_Cell"))
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=6)$y})) # Create a spline that smooths the pseudotime-based expression along 6 degrees of freedom.
filtered_matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
colnames(pt.matrix) <- cell_names
###########
paletteLength <- 20
myColor <- colorRampPalette(c("dodgerblue4", "white", "goldenrod1"))(paletteLength)
# genes_traj <- modulated_genes %>% top_n(., 5,  -p_value)
# ha = rowAnnotation(foo = anno_mark(at = which(rownames(modulated_genes)==rownames(genes_traj)), labels=rownames(modulated_genes),
#         labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
# ha = rowAnnotation(foo = anno_mark(at = which(rownames(modulated_genes)==TCF21_GENES), 
#                                    labels=TCF21_GENES,
#         labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
heatmap_DMC<-pheatmap::pheatmap(pt.matrix, scale="row", 
                      cluster_cols = F, 
                      cluster_rows = TRUE,
                      cutree_rows = 4,
                      fontsize_col = 8,
                      color = myColor,
                      annotation_col = Index,
                      annotation_colors = ann_colors,
                      show_colnames = F,
                      show_rownames = F,
                      # right_annotation = ha,
                      border_color = NA)
# pdf("../2_Output/Cardiomyocyte/Trajectory/Trajectory_Pseudotime_Heatmap.pdf", height = 10, width = 7)
# heatmap_DMC
# dev.off()
################################
hc <-heatmap_DMC$tree_row
lbl <- cutree(hc, 4)
cluster1<-which(lbl==1)
cluster2<-which(lbl==2)
cluster3<-which(lbl==3)
cluster4<-which(lbl==4)
Cluster1_data<-pt.matrix[cluster1,]
Cluster2_data<-pt.matrix[cluster2,]
Cluster3_data<-pt.matrix[cluster3,]
Cluster4_data<-pt.matrix[cluster4,]
Cluster1_GENES <- rownames(Cluster1_data)
Cluster2_GENES <- rownames(Cluster2_data)
Cluster3_GENES <- rownames(Cluster3_data)
Cluster4_GENES <- rownames(Cluster4_data)
##Enrichr
library(enrichR)
dbs <- c("WikiPathways_2019_Mouse")
enriched_1 <- enrichr(Cluster1_GENES, dbs)
enrich_1<-enriched_1[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_1)
enriched_2 <- enrichr(Cluster2_GENES, dbs)
enrich_2<-enriched_2[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_2)
enriched_3 <- enrichr(Cluster3_GENES, dbs)
enrich_3<-enriched_3[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_3)
enriched_4 <- enrichr(Cluster4_GENES, dbs)
enrich_4<-enriched_4[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_4)
library(openxlsx)
wb_DESeq<-createWorkbook()
  addWorksheet(wb_DESeq, "Cluster 1")
  writeData(wb_DESeq, "Cluster 1", enrich_1, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 2")
  writeData(wb_DESeq, "Cluster 2", enrich_2, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 3")
  writeData(wb_DESeq, "Cluster 3", enrich_3, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 4")
  writeData(wb_DESeq, "Cluster 4", enrich_4, startCol = 1)
saveWorkbook(wb_DESeq, file = paste0("../2_Output/Cardiomyocyte/Trajectory/Trajectory_Cardiomyocyte_Pathways.xlsx"), overwrite = TRUE)
library(stringr)
Cluster1_Names<-str_to_title(unlist(strsplit(enrich_1$Genes[1], ";", fixed = T)))
Cluster2_Names<-str_to_title(unlist(strsplit(enrich_2$Genes[1], ";", fixed = T)))
Cluster3_Names<-str_to_title(unlist(strsplit(enrich_3$Genes[1], ";", fixed = T)))
Cluster4_Names<-str_to_title(unlist(strsplit(enrich_4$Genes[1], ";", fixed = T)))
############################################################################
GENES_HM<-c(Cluster1_Names, Cluster2_Names, Cluster3_Names, Cluster4_Names, "Nnt")
pt.df<-as.data.frame(pt.matrix)
pt.df$gene_name<-rownames(pt.matrix)
ha = rowAnnotation(link = anno_mark(at = which(pt.df$gene_name %in% GENES_HM),
                   labels = as.character(pt.df[which(pt.df$gene_name %in% GENES_HM), "gene_name"]),
                   labels_gp = gpar(fontsize = 8),
                   padding = unit(1, "mm"))
                   )
heatmap_combination<-ComplexHeatmap::pheatmap(pt.matrix, scale="row",
                    cluster_cols = F,
                    cluster_rows = T,
                    cutree_rows = 4,
                    # cutree_cols = 3,
                     fontsize_col = 6,
                     color = myColor,
                    annotation_names_col = FALSE,
                    show_colnames = F,
                     show_rownames = F,
                     border_color = NA,
                    annotation_colors = ann_colors,
                    right_annotation = ha,
                    annotation_col = Index,
                    border = TRUE)
heatmap_combination
# pdf(file = "../2_Output/Cardiomyocyte/Trajectory/Trajectory_Heatmap.pdf", height = 7, width = 7)
# heatmap_combination
# dev.off()

## Examine specific genes within the module(s) of interest
library("viridis")
library(ggplot2)
library(ggpubr)
contractility_genes <- c("Mybpc3", "Camk2d", "Pln", "Nnt")
# Transitional_genes <- c("SETD5", "INSR", "OGT", "NOX4", "PLA2G7")
GENES <- contractility_genes
# GENES <- contractility_genes
# GENES <- modulated_genes %>% top_n(., 8,  -pvalue)
# GENES <- rownames(GENES)
p1 <- plot_cells(cds, 
           genes = GENES[1:5],
           label_cell_groups = F,
           label_roots = F,
           label_branch_points = F,
           label_leaves = F,
           show_trajectory_graph = F,
           min_expr = 1,
           alpha = 0.8,
           trajectory_graph_color = "grey28",
           ) +
  # geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme(
  axis.text = element_text(size = 6),    # Adjust the size as needed
  axis.title = element_text(size = 8),  # Adjust the size as needed
  legend.text = element_text(size = 4),  # Adjust the size as needed
  legend.title = element_text(size = 6),
  legend.key.size = unit(2, 'mm')) +
  scale_colour_viridis_c(option = "inferno")
pdf(file="../2_Output/UMAP_Pseudotime_Genes.pdf", height = 3, width = 3)
# ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="right")
p1
dev.off()
# Plot genes according to "pseudotime"
lineage_cds <- cds[rowData(cds)$gene_short_name %in% GENES, ] #colData(cds)$cell_type %in% c("VSMC")
lineage_cds<-order_cells(lineage_cds, root_cells = root_cell)
pdf(file="../2_Output/Pseudotime_curves.pdf", height = 9, width = 3)
monocle3::plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="ident",
                         min_expr=0)
#### Scvelo
# Find the common set of cells between the data matrix and metadata
common_cells <- intersect(colnames(hdat_Cardiomyocyte@assays$SCT@data), rownames(hdat_Cardiomyocyte@meta.data))
# Subset the Seurat object to include only these common cells
Merged_scaled <- subset(hdat_Cardiomyocyte, cells = common_cells)
# Now, check again to ensure that the number of cells in data and metadata match
dim(hdat_Cardiomyocyte@assays$SCT@data)  # Should match the rows of the metadata
dim(hdat_Cardiomyocyte@meta.data)
# Subset Seurat object to ensure metadata and data match
cell_names <- colnames(hdat_Cardiomyocyte@assays$SCT@counts)
hdat_Cardiomyocyte <- subset(hdat_Cardiomyocyte, cells = cell_names)
# Export the raw count matrix
counts_matrix <- as.matrix(t(hdat_Cardiomyocyte@assays$SCT@counts))
write.csv(counts_matrix, file = "../2_Output/Cardiomyocyte/Trajectory/python/counts_matrix.csv")
# Export the normalized data matrix
data_matrix <- as.matrix(t(hdat_Cardiomyocyte@assays$SCT@data))
write.csv(data_matrix, file = "../2_Output/Cardiomyocyte/Trajectory/python/normalized_data_matrix.csv")
# Export cell metadata
cell_metadata <- hdat_Cardiomyocyte@meta.data
write.csv(cell_metadata, file = "../2_Output/Cardiomyocyte/Trajectory/python/cell_metadata.csv")
# Export gene metadata (optional, could be just rownames)
gene_metadata <- data.frame(genes = rownames(hdat_Cardiomyocyte@assays$SCT@data))
write.csv(gene_metadata, file = "../2_Output/Cardiomyocyte/Trajectory/python/gene_metadata.csv")
```

# ShinyR App Development - ShinyCell Scaffold

The following script was used to prepare the unsupervised analysis for
deployment on a Shiny Server to make all data and analyses widely
accessible, tailored specifically for coding-naive basic scientists.
This is now deployed on
<https://markpepin.shinyapps.io/HFpEF_Transcriptomics/>


``` r
library(Seurat)
library(ShinyCell)
library(qs)
TAA.combined<-qread(file = "../1_Input/NNT_Labelling_snRNA.qs")
TAA.combined$Group <- paste(TAA.combined$Background, TAA.combined$Treatment, sep = "_")
# saveRDS(TAA.combined, "../1_Input/NNT_Labelling_snRNA.rds")
# Make a config file
scConf <- createConfig(
  obj = TAA.combined,                # Your Seurat object
  meta.to.include = c("CellType","Group", "Background", "Treatment", "SCT_snn_res.0.2", "SampleID"),  # Metadata columns to include
  legendCols = 4,                    # Number of columns in categorical metadata legends
  maxLevels = 50                     # Max number of levels for categorical metadata
)
# makeShinyApp(TAA.combined, scConf, gene.mapping = TRUE,
#              shiny.title = "ShinyCell Quick Start") 
makeShinyApp(
  obj = TAA.combined,               # Seurat object or file path to the dataset
  scConf = scConf,                  # Configuration data.table generated in Step 1
  default.dimred = c("UMAP_1", "UMAP_2"),
  shiny.footnotes = "Pepin et al. 2024",
  gex.assay = "SCT",                # Assay in Seurat object to use (e.g., "RNA", "integrated")
  gex.slot = "data",                # Slot in the Seurat assay to use ("data" is default for normalized expression)
  shiny.title = "NNT-dependent Cardiometabolic HFpEF: snRNA-Sequencing Analysis",  # Title for the Shiny app
  shiny.dir = "HFpEF_APP",       # Directory to create the app files
  enableSubset = TRUE,              # Enable subsetting cells functionality in the app
  defPtSiz = 1.25,                  # Default point size for single cells in plots
  default.gene1 = "Myh7",          # Primary default gene to show in feature plots
  default.gene2 = "Nppb",            # Secondary default gene for comparison
  default.multigene = c("Vwf", "Pecam1", "Cdh5", "Plp1", "Nrxn1", "Nrxn3", "Dcn", "Gsn", "Pdgfra", "Gpam", "Fasn", "Lep", "Msln", "Wt1", "Rgs5", "Abcc9", "Kcnj8", "Myh11", "Acta2", "Tagln", "Nppa","Myl4","Myh7", "Myl2", "Fhl2", "Cd5","Cd14", "C1qa", "Cd68")
)
```

# Supplemental Table: R Session Information


``` r
end_time <- Sys.time()
# execution_time <- end_time - start_time
sinfo<-devtools::session_info()
sinfo$platform
```

```
##  setting  value
##  version  R version 4.4.1 (2024-06-14)
##  os       macOS Sonoma 14.6.1
##  system   aarch64, darwin20
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/Los_Angeles
##  date     2024-10-11
##  pandoc   3.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
```

``` r
sinfo$packages %>% kable( 
                         align="c", 
                         longtable=T, 
                         booktabs=T,
                         caption="Packages and Required Dependencies") %>% 
    kable_styling(latex_options=c("striped", "repeat_header", "condensed"))
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Packages and Required Dependencies</caption>
 <thead>
  <tr>
   <th style="text-align:left;">  </th>
   <th style="text-align:center;"> package </th>
   <th style="text-align:center;"> ondiskversion </th>
   <th style="text-align:center;"> loadedversion </th>
   <th style="text-align:center;"> path </th>
   <th style="text-align:center;"> loadedpath </th>
   <th style="text-align:center;"> attached </th>
   <th style="text-align:center;"> is_base </th>
   <th style="text-align:center;"> date </th>
   <th style="text-align:center;"> source </th>
   <th style="text-align:center;"> md5ok </th>
   <th style="text-align:center;"> library </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> abind </td>
   <td style="text-align:center;"> abind </td>
   <td style="text-align:center;"> 1.4.8 </td>
   <td style="text-align:center;"> 1.4-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/abind </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/abind </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bslib </td>
   <td style="text-align:center;"> bslib </td>
   <td style="text-align:center;"> 0.8.0 </td>
   <td style="text-align:center;"> 0.8.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bslib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bslib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cachem </td>
   <td style="text-align:center;"> cachem </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cachem </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cachem </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cli </td>
   <td style="text-align:center;"> cli </td>
   <td style="text-align:center;"> 3.6.3 </td>
   <td style="text-align:center;"> 3.6.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cli </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cli </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cluster </td>
   <td style="text-align:center;"> cluster </td>
   <td style="text-align:center;"> 2.1.6 </td>
   <td style="text-align:center;"> 2.1.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cluster </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cluster </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> codetools </td>
   <td style="text-align:center;"> codetools </td>
   <td style="text-align:center;"> 0.2.20 </td>
   <td style="text-align:center;"> 0.2-20 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/codetools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/codetools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-31 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> colorspace </td>
   <td style="text-align:center;"> colorspace </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> 2.1-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/colorspace </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/colorspace </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cowplot </td>
   <td style="text-align:center;"> cowplot </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cowplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cowplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:center;"> data.table </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/data.table </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/data.table </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> deldir </td>
   <td style="text-align:center;"> deldir </td>
   <td style="text-align:center;"> 2.0.4 </td>
   <td style="text-align:center;"> 2.0-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/deldir </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/deldir </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> devtools </td>
   <td style="text-align:center;"> devtools </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/devtools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/devtools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-11 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> digest </td>
   <td style="text-align:center;"> digest </td>
   <td style="text-align:center;"> 0.6.37 </td>
   <td style="text-align:center;"> 0.6.37 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/digest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/digest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dotCall64 </td>
   <td style="text-align:center;"> dotCall64 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dotCall64 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dotCall64 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dplyr </td>
   <td style="text-align:center;"> dplyr </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dplyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dplyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ellipsis </td>
   <td style="text-align:center;"> ellipsis </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ellipsis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ellipsis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> evaluate </td>
   <td style="text-align:center;"> evaluate </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/evaluate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/evaluate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fansi </td>
   <td style="text-align:center;"> fansi </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fansi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fansi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> farver </td>
   <td style="text-align:center;"> farver </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/farver </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/farver </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastDummies </td>
   <td style="text-align:center;"> fastDummies </td>
   <td style="text-align:center;"> 1.7.4 </td>
   <td style="text-align:center;"> 1.7.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastDummies </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastDummies </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastmap </td>
   <td style="text-align:center;"> fastmap </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmap </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fitdistrplus </td>
   <td style="text-align:center;"> fitdistrplus </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fs </td>
   <td style="text-align:center;"> fs </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-25 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future </td>
   <td style="text-align:center;"> future </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future.apply </td>
   <td style="text-align:center;"> future.apply </td>
   <td style="text-align:center;"> 1.11.2 </td>
   <td style="text-align:center;"> 1.11.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future.apply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future.apply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> generics </td>
   <td style="text-align:center;"> generics </td>
   <td style="text-align:center;"> 0.1.3 </td>
   <td style="text-align:center;"> 0.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/generics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/generics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggplot2 </td>
   <td style="text-align:center;"> ggplot2 </td>
   <td style="text-align:center;"> 3.5.1 </td>
   <td style="text-align:center;"> 3.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggrepel </td>
   <td style="text-align:center;"> ggrepel </td>
   <td style="text-align:center;"> 0.9.6 </td>
   <td style="text-align:center;"> 0.9.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrepel </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrepel </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggridges </td>
   <td style="text-align:center;"> ggridges </td>
   <td style="text-align:center;"> 0.5.6 </td>
   <td style="text-align:center;"> 0.5.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggridges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggridges </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> globals </td>
   <td style="text-align:center;"> globals </td>
   <td style="text-align:center;"> 0.16.3 </td>
   <td style="text-align:center;"> 0.16.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/globals </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/globals </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> glue </td>
   <td style="text-align:center;"> glue </td>
   <td style="text-align:center;"> 1.8.0 </td>
   <td style="text-align:center;"> 1.8.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/glue </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/glue </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> goftest </td>
   <td style="text-align:center;"> goftest </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> 1.2-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/goftest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/goftest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-10-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gridExtra </td>
   <td style="text-align:center;"> gridExtra </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gridExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gridExtra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2017-09-09 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gtable </td>
   <td style="text-align:center;"> gtable </td>
   <td style="text-align:center;"> 0.3.5 </td>
   <td style="text-align:center;"> 0.3.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> harmony </td>
   <td style="text-align:center;"> harmony </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/harmony </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/harmony </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmltools </td>
   <td style="text-align:center;"> htmltools </td>
   <td style="text-align:center;"> 0.5.8.1 </td>
   <td style="text-align:center;"> 0.5.8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmltools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmltools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmlwidgets </td>
   <td style="text-align:center;"> htmlwidgets </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httpuv </td>
   <td style="text-align:center;"> httpuv </td>
   <td style="text-align:center;"> 1.6.15 </td>
   <td style="text-align:center;"> 1.6.15 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httpuv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httpuv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httr </td>
   <td style="text-align:center;"> httr </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ica </td>
   <td style="text-align:center;"> ica </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> 1.0-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ica </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ica </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> igraph </td>
   <td style="text-align:center;"> igraph </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/igraph </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/igraph </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> irlba </td>
   <td style="text-align:center;"> irlba </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/irlba </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/irlba </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-03 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jquerylib </td>
   <td style="text-align:center;"> jquerylib </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jquerylib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jquerylib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jsonlite </td>
   <td style="text-align:center;"> jsonlite </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jsonlite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jsonlite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> kableExtra </td>
   <td style="text-align:center;"> kableExtra </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/kableExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/kableExtra </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KernSmooth </td>
   <td style="text-align:center;"> KernSmooth </td>
   <td style="text-align:center;"> 2.23.24 </td>
   <td style="text-align:center;"> 2.23-24 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> knitr </td>
   <td style="text-align:center;"> knitr </td>
   <td style="text-align:center;"> 1.48 </td>
   <td style="text-align:center;"> 1.48 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/knitr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/knitr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> later </td>
   <td style="text-align:center;"> later </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/later </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/later </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lattice </td>
   <td style="text-align:center;"> lattice </td>
   <td style="text-align:center;"> 0.22.6 </td>
   <td style="text-align:center;"> 0.22-6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lattice </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lattice </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lazyeval </td>
   <td style="text-align:center;"> lazyeval </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lazyeval </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lazyeval </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> leiden </td>
   <td style="text-align:center;"> leiden </td>
   <td style="text-align:center;"> 0.4.3.1 </td>
   <td style="text-align:center;"> 0.4.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/leiden </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/leiden </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lifecycle </td>
   <td style="text-align:center;"> lifecycle </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lifecycle </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lifecycle </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> listenv </td>
   <td style="text-align:center;"> listenv </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/listenv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/listenv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lmtest </td>
   <td style="text-align:center;"> lmtest </td>
   <td style="text-align:center;"> 0.9.40 </td>
   <td style="text-align:center;"> 0.9-40 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lmtest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lmtest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> magrittr </td>
   <td style="text-align:center;"> magrittr </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/magrittr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/magrittr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS </td>
   <td style="text-align:center;"> MASS </td>
   <td style="text-align:center;"> 7.3.61 </td>
   <td style="text-align:center;"> 7.3-61 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MASS </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MASS </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Matrix </td>
   <td style="text-align:center;"> Matrix </td>
   <td style="text-align:center;"> 1.7.0 </td>
   <td style="text-align:center;"> 1.7-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Matrix </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Matrix </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> matrixStats </td>
   <td style="text-align:center;"> matrixStats </td>
   <td style="text-align:center;"> 1.4.1 </td>
   <td style="text-align:center;"> 1.4.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/matrixStats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/matrixStats </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> memoise </td>
   <td style="text-align:center;"> memoise </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/memoise </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/memoise </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mime </td>
   <td style="text-align:center;"> mime </td>
   <td style="text-align:center;"> 0.12 </td>
   <td style="text-align:center;"> 0.12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mime </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mime </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-09-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miniUI </td>
   <td style="text-align:center;"> miniUI </td>
   <td style="text-align:center;"> 0.1.1.1 </td>
   <td style="text-align:center;"> 0.1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/miniUI </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/miniUI </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-05-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> munsell </td>
   <td style="text-align:center;"> munsell </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/munsell </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/munsell </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nlme </td>
   <td style="text-align:center;"> nlme </td>
   <td style="text-align:center;"> 3.1.166 </td>
   <td style="text-align:center;"> 3.1-166 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/nlme </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/nlme </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-14 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> parallelly </td>
   <td style="text-align:center;"> parallelly </td>
   <td style="text-align:center;"> 1.38.0 </td>
   <td style="text-align:center;"> 1.38.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/parallelly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/parallelly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> patchwork </td>
   <td style="text-align:center;"> patchwork </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/patchwork </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/patchwork </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pbapply </td>
   <td style="text-align:center;"> pbapply </td>
   <td style="text-align:center;"> 1.7.2 </td>
   <td style="text-align:center;"> 1.7-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pbapply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pbapply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-06-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pillar </td>
   <td style="text-align:center;"> pillar </td>
   <td style="text-align:center;"> 1.9.0 </td>
   <td style="text-align:center;"> 1.9.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pillar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pillar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgbuild </td>
   <td style="text-align:center;"> pkgbuild </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgconfig </td>
   <td style="text-align:center;"> pkgconfig </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgload </td>
   <td style="text-align:center;"> pkgload </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgload </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgload </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plotly </td>
   <td style="text-align:center;"> plotly </td>
   <td style="text-align:center;"> 4.10.4 </td>
   <td style="text-align:center;"> 4.10.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plotly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plotly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plyr </td>
   <td style="text-align:center;"> plyr </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> png </td>
   <td style="text-align:center;"> png </td>
   <td style="text-align:center;"> 0.1.8 </td>
   <td style="text-align:center;"> 0.1-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/png </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/png </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polyclip </td>
   <td style="text-align:center;"> polyclip </td>
   <td style="text-align:center;"> 1.10.7 </td>
   <td style="text-align:center;"> 1.10-7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/polyclip </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/polyclip </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> profvis </td>
   <td style="text-align:center;"> profvis </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/profvis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/profvis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> progressr </td>
   <td style="text-align:center;"> progressr </td>
   <td style="text-align:center;"> 0.14.0 </td>
   <td style="text-align:center;"> 0.14.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/progressr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/progressr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> promises </td>
   <td style="text-align:center;"> promises </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/promises </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/promises </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> purrr </td>
   <td style="text-align:center;"> purrr </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/purrr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/purrr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R6 </td>
   <td style="text-align:center;"> R6 </td>
   <td style="text-align:center;"> 2.5.1 </td>
   <td style="text-align:center;"> 2.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R6 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RANN </td>
   <td style="text-align:center;"> RANN </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RANN </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RANN </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-25 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RColorBrewer </td>
   <td style="text-align:center;"> RColorBrewer </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> 1.1-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-04-03 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rcpp </td>
   <td style="text-align:center;"> Rcpp </td>
   <td style="text-align:center;"> 1.0.13 </td>
   <td style="text-align:center;"> 1.0.13 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rcpp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rcpp </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RcppAnnoy </td>
   <td style="text-align:center;"> RcppAnnoy </td>
   <td style="text-align:center;"> 0.0.22 </td>
   <td style="text-align:center;"> 0.0.22 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RcppHNSW </td>
   <td style="text-align:center;"> RcppHNSW </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppHNSW </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppHNSW </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> remotes </td>
   <td style="text-align:center;"> remotes </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/remotes </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/remotes </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reshape2 </td>
   <td style="text-align:center;"> reshape2 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reshape2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reshape2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-09 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reticulate </td>
   <td style="text-align:center;"> reticulate </td>
   <td style="text-align:center;"> 1.39.0 </td>
   <td style="text-align:center;"> 1.39.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reticulate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reticulate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rlang </td>
   <td style="text-align:center;"> rlang </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rlang </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rlang </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rmarkdown </td>
   <td style="text-align:center;"> rmarkdown </td>
   <td style="text-align:center;"> 2.28 </td>
   <td style="text-align:center;"> 2.28 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ROCR </td>
   <td style="text-align:center;"> ROCR </td>
   <td style="text-align:center;"> 1.0.11 </td>
   <td style="text-align:center;"> 1.0-11 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ROCR </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ROCR </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RSpectra </td>
   <td style="text-align:center;"> RSpectra </td>
   <td style="text-align:center;"> 0.16.2 </td>
   <td style="text-align:center;"> 0.16-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RSpectra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RSpectra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstudioapi </td>
   <td style="text-align:center;"> rstudioapi </td>
   <td style="text-align:center;"> 0.16.0 </td>
   <td style="text-align:center;"> 0.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rtsne </td>
   <td style="text-align:center;"> Rtsne </td>
   <td style="text-align:center;"> 0.17 </td>
   <td style="text-align:center;"> 0.17 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rtsne </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rtsne </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sass </td>
   <td style="text-align:center;"> sass </td>
   <td style="text-align:center;"> 0.4.9 </td>
   <td style="text-align:center;"> 0.4.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sass </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scales </td>
   <td style="text-align:center;"> scales </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scales </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scales </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scattermore </td>
   <td style="text-align:center;"> scattermore </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scattermore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scattermore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-06-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sctransform </td>
   <td style="text-align:center;"> sctransform </td>
   <td style="text-align:center;"> 0.4.1 </td>
   <td style="text-align:center;"> 0.4.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sctransform </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sctransform </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sessioninfo </td>
   <td style="text-align:center;"> sessioninfo </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Seurat </td>
   <td style="text-align:center;"> Seurat </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Seurat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Seurat </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeuratObject </td>
   <td style="text-align:center;"> SeuratObject </td>
   <td style="text-align:center;"> 5.0.2 </td>
   <td style="text-align:center;"> 5.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shiny </td>
   <td style="text-align:center;"> shiny </td>
   <td style="text-align:center;"> 1.9.1 </td>
   <td style="text-align:center;"> 1.9.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shiny </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shiny </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sp </td>
   <td style="text-align:center;"> sp </td>
   <td style="text-align:center;"> 2.1.4 </td>
   <td style="text-align:center;"> 2.1-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sp </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spam </td>
   <td style="text-align:center;"> spam </td>
   <td style="text-align:center;"> 2.10.0 </td>
   <td style="text-align:center;"> 2.10-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spam </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spam </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.data </td>
   <td style="text-align:center;"> spatstat.data </td>
   <td style="text-align:center;"> 3.1.2 </td>
   <td style="text-align:center;"> 3.1-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.explore </td>
   <td style="text-align:center;"> spatstat.explore </td>
   <td style="text-align:center;"> 3.3.2 </td>
   <td style="text-align:center;"> 3.3-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.geom </td>
   <td style="text-align:center;"> spatstat.geom </td>
   <td style="text-align:center;"> 3.3.3 </td>
   <td style="text-align:center;"> 3.3-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.random </td>
   <td style="text-align:center;"> spatstat.random </td>
   <td style="text-align:center;"> 3.3.2 </td>
   <td style="text-align:center;"> 3.3-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.sparse </td>
   <td style="text-align:center;"> spatstat.sparse </td>
   <td style="text-align:center;"> 3.1.0 </td>
   <td style="text-align:center;"> 3.1-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.univar </td>
   <td style="text-align:center;"> spatstat.univar </td>
   <td style="text-align:center;"> 3.0.1 </td>
   <td style="text-align:center;"> 3.0-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.univar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.univar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.utils </td>
   <td style="text-align:center;"> spatstat.utils </td>
   <td style="text-align:center;"> 3.1.0 </td>
   <td style="text-align:center;"> 3.1-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringi </td>
   <td style="text-align:center;"> stringi </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringr </td>
   <td style="text-align:center;"> stringr </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-14 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> survival </td>
   <td style="text-align:center;"> survival </td>
   <td style="text-align:center;"> 3.7.0 </td>
   <td style="text-align:center;"> 3.7-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/survival </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/survival </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> svglite </td>
   <td style="text-align:center;"> svglite </td>
   <td style="text-align:center;"> 2.1.3 </td>
   <td style="text-align:center;"> 2.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/svglite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/svglite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> systemfonts </td>
   <td style="text-align:center;"> systemfonts </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/systemfonts </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/systemfonts </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tensor </td>
   <td style="text-align:center;"> tensor </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tensor </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tensor </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2012-05-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tibble </td>
   <td style="text-align:center;"> tibble </td>
   <td style="text-align:center;"> 3.2.1 </td>
   <td style="text-align:center;"> 3.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tibble </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tibble </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyr </td>
   <td style="text-align:center;"> tidyr </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyselect </td>
   <td style="text-align:center;"> tidyselect </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyselect </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyselect </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-11 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> urlchecker </td>
   <td style="text-align:center;"> urlchecker </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/urlchecker </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/urlchecker </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> usethis </td>
   <td style="text-align:center;"> usethis </td>
   <td style="text-align:center;"> 3.0.0 </td>
   <td style="text-align:center;"> 3.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/usethis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/usethis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> utf8 </td>
   <td style="text-align:center;"> utf8 </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/utf8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/utf8 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> uwot </td>
   <td style="text-align:center;"> uwot </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/uwot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/uwot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vctrs </td>
   <td style="text-align:center;"> vctrs </td>
   <td style="text-align:center;"> 0.6.5 </td>
   <td style="text-align:center;"> 0.6.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vctrs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vctrs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> viridisLite </td>
   <td style="text-align:center;"> viridisLite </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/viridisLite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/viridisLite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xfun </td>
   <td style="text-align:center;"> xfun </td>
   <td style="text-align:center;"> 0.47 </td>
   <td style="text-align:center;"> 0.47 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xfun </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xfun </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xml2 </td>
   <td style="text-align:center;"> xml2 </td>
   <td style="text-align:center;"> 1.3.6 </td>
   <td style="text-align:center;"> 1.3.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xml2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xml2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xtable </td>
   <td style="text-align:center;"> xtable </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yaml </td>
   <td style="text-align:center;"> yaml </td>
   <td style="text-align:center;"> 2.3.10 </td>
   <td style="text-align:center;"> 2.3.10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yaml </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yaml </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zoo </td>
   <td style="text-align:center;"> zoo </td>
   <td style="text-align:center;"> 1.8.12 </td>
   <td style="text-align:center;"> 1.8-12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zoo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zoo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-04-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
</tbody>
</table>
