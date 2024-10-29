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
[mpepin\@stanford.edu](mailto:mpepin@stanford.edu){.email}\
**Affiliations**: Heidelberg University Hospital, Institute for Experimental Cardiology | Stanford University, Division of Cardiovascular Medicine \| Broad Institute of
Harvard and MIT\
**Location**: Stanford, CA

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
options(future.globals.maxSize = 16 * 1024^3)  # Set the future global memory limit to 16 GB
library(Seurat)
library(openxlsx)
library(dplyr)
library(DoubletFinder)
# Read sample info
Index <- read.xlsx("../1_Input/M020_Sample.Info_NvJ_mep.xlsx")
# Initialize an empty list to store Seurat objects
seurat_list <- list()

# Loop through all the samples in Index
for (i in 1:nrow(Index)) {
  data_dir <- paste0("../1_Input/snRNA/", Index$Treatment[i], "_", Index$Background[i], "_", Index$MouseID[i], "/filtered_feature_bc_matrix/")
  
  # Read the data
  snRNA_data <- Read10X(data.dir = data_dir)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = snRNA_data, project = "snRNA_NNT", min.cells = 3, min.features = 200)
  seurat_obj[["SampleID"]] <- Index$MouseID[i]
  seurat_obj[["Background"]] <- Index$Background[i]
  seurat_obj[["Treatment"]] <- Index$Treatment[i]
  
  # Add mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, assay = "RNA", pattern = "^mt-")
  
  # Subset data
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 5)
  
  # Normalization, variable features, and scaling
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 10000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- SCTransform(seurat_obj)
  
  # PCA
  seurat_obj <- RunPCA(seurat_obj, npcs = 50)
  
  # Perform clustering (necessary for DoubletFinder)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:40)  # Finding neighbors
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.2)  # Find clusters

  # DoubletFinder parameter sweep
  sweep_res <- paramSweep(seurat_obj, PCs = 1:40, sct = TRUE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  optimal_pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  
  # Estimate homotypic doublet proportion and expected number of doublets
  homotypic_prop <- modelHomotypic(seurat_obj$seurat_clusters)  # Adjust clusters if needed
  nExp_poi <- round(0.075 * nrow(seurat_obj@meta.data))  # Assuming 7.5% doublet rate
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  
  # DoubletFinder
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:40, pN = 0.25, pK = optimal_pK, nExp = nExp_poi_adj, reuse.pANN = FALSE, sct = TRUE)
  
  # Filter out doublets
 # Find the column name that starts with "DF.classifications"
df_class_col <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
# Add column that computes the doublet perentage
doublet_percentage <- seurat_obj@meta.data %>%
  group_by(SampleID) %>%
  summarise(PercentDoublets = mean((!!as.name(df_class_col)) == "Doublet") * 100)
# Subset based on the column dynamically
seurat_obj <- subset(seurat_obj, subset = (!!as.name(df_class_col)) == "Singlet")

# Step 4: Merge this doublet percentage information back into the Seurat metadata
seurat_obj@meta.data <- merge(seurat_obj@meta.data, doublet_percentage, by = "SampleID", all.x = TRUE)

  # Store the processed Seurat object in the list
  seurat_list[[i]] <- seurat_obj
}

# Ensure unique cell names and matching metadata rownames across Seurat objects
for (i in seq_along(seurat_list)) {
    seurat_list[[i]] <- RenameCells(object = seurat_list[[i]], add.cell.id = Index$MouseID[i])
    rownames(seurat_list[[i]]@meta.data) <- Cells(seurat_list[[i]])  # Sync metadata rownames with cell names
}

# Select features for integration
integ_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

# Merge Seurat objects
merged_seurat <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)], merge.data = TRUE)

# Set assay and variable features for merged object
DefaultAssay(merged_seurat) <- "SCT"
VariableFeatures(merged_seurat) <- integ_features

# Run PCA
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)

# Perform Harmony integration
library(harmony)
harmonized_seurat <- RunHarmony(merged_seurat, group.by.vars = "SampleID", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

# UMAP, neighbors, and clustering on Harmony results
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:40)
harmonized_seurat <- FindNeighbors(harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))

# Prepare for marker identification
harmonized_seurat <- PrepSCTFindMarkers(harmonized_seurat, assay = "SCT", verbose = TRUE)

# Save the final object
saveRDS(harmonized_seurat, file = "../1_Input/NNT_Integration_snRNAv2.rds")

# Record end time
end_time <- Sys.time()
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
Approximation and Projection (UMAP) was used for dimensional
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
snRNA.combined<-readRDS(file = "../1_Input/NNT_Integration_snRNAv2.rds")
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
seurat_obj <- readRDS(file = "../1_Input/NNT_Integration_snRNAv2.rds")
DefaultAssay(seurat_obj) <- "SCT"  # Specify the assay you want
hESCs <- as.SingleCellExperiment(seurat_obj, assay = "SCT")
# hESCs <- as.SingleCellExperiment(readRDS(file = "../1_Input/NNT_Integration_snRNAv2.rds"))
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
seurat_obj[["SingleR.labels"]] <- pred.hesSingleRpred.hesc$labels
Idents(seurat_obj) <- seurat_obj$SingleR.labels
UMAP_Treatment<-DimPlot(seurat_obj,  label = T, split.by = "Treatment") + NoLegend()
UMAP_Background<-DimPlot(seurat_obj,  label = T, split.by = "Background")
UMAP_Treatment + UMAP_Background
pdf(file = "../2_Output/UMAP_Clusters_unbiased.annotation.pdf", height = 5, width = 9)
UMAP_Treatment
UMAP_Background
dev.off()
saveRDS(snRNA.combined, , file = "../1_Input/snRNA_unbiased.Annnotation_snRNAv2.rds")
# Annotation performance with SingleR
## heatmap of annotation scores
pdf("../2_Output/annotation.heatmap.pdf")
plotScoreHeatmap(pred.hesSingleRpred.hesc)
dev.off()
## Delta Distribution
plotDeltaDistribution(pred.hesSingleRpred.hesc, ncol = 3)
```

# Gene Markers to Validate Cell Type Identification

We manually inspected cell-type-specific gene markers using density and ridge plots to validate cell identity that was assigned by the using the automated SingleR package (above). For each cell type, we plotted the combined expression of multiple marker genes. Additionally, we created individual density plots for representative genes within each cell type. To further analyze the distribution of marker genes, we used ridge plots to understand the distribution of expression within each cell type. The resulting Seurat object was subsequently saved for downstream analysis.


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
dev.off()

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

We created UMAP plots highlighting clustering by CellType and Treatment, followed by proportional bar plots that examined cell-type distributions across sample groups and treatment conditions. These proportional data were exported to CSV and Excel for further analysis. We also visualized PCA loadings and generated dimensional heatmaps, saving all plots as PDFs.


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


# Venn Diagram

We used differential expression analysis to compare gene expression of HFpEF (HFD+L-NAME) vs. Ctrl between two backgrounds (N vs. J). We then created matrices indicating the presence of significant DE genes across cell types for HFpEF and Ctrl, respectively, and visualized these patterns using "UpSet" plots (a high-order venn diagram). We performed Reactome and GO pathway enrichment analyses for each cell type, visualizing the results in a dot plot, showing enriched pathways for each cell type across conditions.


``` r
options(future.globals.maxSize = 8 * 1024^3)  # 2GB
library(Seurat)
library(dplyr)
library(openxlsx)
library(UpSetR)

# Load the Seurat object
TAA.combined <- readRDS(file = "../1_Input/NNT_Labelling_snRNA.rds")
# Define unique cell types
cell_states <- unique(TAA.combined$CellType)

# Initialize lists for storing DEGs per condition and cell type
deg_list_HFpEF <- list()
deg_list_Ctrl <- list()

# Loop through each cell type and calculate DEGs for HFpEF and Ctrl conditions separately
for (CELL in cell_states) {
  # Subset the Seurat object for HFpEF and Ctrl conditions
  TAA_HFpEF <- subset(TAA.combined, Treatment == "HFpEF" & CellType == CELL)
  TAA_Ctrl <- subset(TAA.combined, Treatment == "Ctrl" & CellType == CELL)
  
  # Recalculate SCT model to reset it
  TAA_HFpEF <- SCTransform(TAA_HFpEF, verbose = TRUE, assay = "RNA", return.only.var.genes = FALSE)
  TAA_Ctrl <- SCTransform(TAA_Ctrl, verbose = TRUE, assay = "RNA", return.only.var.genes = FALSE)

  # Prepare the SCT assay for FindMarkers
  TAA_HFpEF <- PrepSCTFindMarkers(TAA_HFpEF)
  TAA_Ctrl <- PrepSCTFindMarkers(TAA_Ctrl)
  
  # Find DEGs between "N" and "J" in HFpEF
  deg_HFpEF <- FindMarkers(TAA_HFpEF, ident.1 = "N", ident.2 = "J", group.by = "Background", logfc.threshold = 0, min.pct = 0.1)
  deg_HFpEF_filtered <- deg_HFpEF %>% filter(p_val_adj < 0.05)  # Filter significant DEGs
  deg_list_HFpEF[[CELL]] <- rownames(deg_HFpEF_filtered)  # Store significant DEGs
  
  # Find DEGs between "N" and "J" in Ctrl
  deg_Ctrl <- FindMarkers(TAA_Ctrl, ident.1 = "N", ident.2 = "J", group.by = "Background", logfc.threshold = 0, min.pct = 0.1)
  deg_Ctrl_filtered <- deg_Ctrl %>% filter(p_val_adj < 0.05)  # Filter significant DEGs
  deg_list_Ctrl[[CELL]] <- rownames(deg_Ctrl_filtered)  # Store significant DEGs
}

# Create a combined DEG matrix for UpSet plotting
deg_matrix_HFpEF <- data.frame(
  gene = unique(unlist(deg_list_HFpEF)),  # Unique genes for HFpEF condition
  stringsAsFactors = FALSE
)

deg_matrix_Ctrl <- data.frame(
  gene = unique(unlist(deg_list_Ctrl)),  # Unique genes for Ctrl condition
  stringsAsFactors = FALSE
)

# Add binary indicators for DEGs presence/absence in each cell type (HFpEF and Ctrl)
for (cell_type in names(deg_list_HFpEF)) {
  deg_matrix_HFpEF[[cell_type]] <- ifelse(deg_matrix_HFpEF$gene %in% deg_list_HFpEF[[cell_type]], 1, 0)
}

for (cell_type in names(deg_list_Ctrl)) {
  deg_matrix_Ctrl[[cell_type]] <- ifelse(deg_matrix_Ctrl$gene %in% deg_list_Ctrl[[cell_type]], 1, 0)
}

# Define the order of cell types
desired_order <- rev(c("Cardiomyocyte", "Fibroblast", "Mural_Cell", "EC", "Myeloid", "Lymphoid", "Mast_Cell", "Neural_Cell"))

# Plot UpSet for HFpEF condition
upset_HFpEF <- upset(deg_matrix_HFpEF, 
      sets = desired_order, 
      mb.ratio = c(0.4, 0.6),
      keep.order = TRUE,  
      main.bar.color = "#56B4E9",  
      sets.bar.color = "#D55E00",  
      order.by = "freq")  

# Plot UpSet for Ctrl condition
upset_Ctrl <- upset(deg_matrix_Ctrl, 
      sets = desired_order, 
      mb.ratio = c(0.4, 0.6),
      keep.order = TRUE,  
      main.bar.color = "#56B4E9",  
      sets.bar.color = "#D55E00",  
      order.by = "freq")

# Save the plots as PDF
pdf("../2_Output/UpSetR_HFpEF.pdf", height = 7, width = 7)
print(upset_HFpEF)
dev.off()

pdf("../2_Output/UpSetR_Ctrl.pdf", height = 7, width = 7)
print(upset_Ctrl)
dev.off()

# Load libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ReactomePA)
deg_list <- list()
for (CELL in cell_states) {
  markers <- openxlsx::read.xlsx(paste0("./HFpEF_", CELL, "_DEGs.xlsx"), rowNames = T) %>% top_n(n = 1000, wt = abs(log2FoldChange))
  markers <- markers[!is.na(rownames(markers)), ]
  deg_list[[CELL]] <- rownames(markers[markers$pvalue < 0.05, ])  # Extract DEGs with p-value < 0.05
}
# Prepare DEG list (Example DEG list structure)
common_genes <- Reduce(intersect, deg_list)
# Convert gene symbols to Entrez IDs
deg_list_entrez <- lapply(deg_list, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})
# Perform Reactome pathway enrichment analysis for each cell type
reactome_results <- lapply(deg_list_entrez, function(entrez_genes) {
  enrichPathway(gene = entrez_genes, organism = "mouse", qvalueCutoff = 0.05)
})
# Assign cell type names
names(reactome_results) <- names(deg_list)
# Create a compareCluster object for Reactome pathway enrichment
compare_cluster_results <- compareCluster(
  geneCluster = deg_list_entrez,
  fun = "enrichGO",
  OrgDb = "org.Mm.eg.db",
  # organism = "mm",
  pvalueCutoff = 0.05
)

compare_cluster_results@compareClusterResult$cell_type <- compare_cluster_results@compareClusterResult$Cluster
pdf("../2_Output/DotPlot_SPLIT.pdf", height = 4, width = 11)
dotplot(compare_cluster_results, showCategory = 3, x = "GeneRatio", split = "cell_type", label_format = 50) +
  facet_grid(. ~ cell_type) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10))
dev.off()
```

# Cell-Cell Interactome

This code initializes multiple analyses of cell-cell interactions, gene expression, and differential expression between cell types under different treatment conditions. Using LIANA for ligand-receptor interaction analysis of a seurat-generated single-cell dataset, we identify significant interactions, visualize these with heatmaps and chord diagrams, and plot the expression patterns of specific receptor genes across both Treatment (HFpEF vs Ctrl) and genetic background (N vs. J). Additionally, violin plots for expression patterns of multiple genes of interest are generated, with customized themes for better visual clarity. PDF files for each plot are saved to the specified directories.


``` r
library(Seurat)
library(tidyverse)
library(magrittr)
library(liana)
## Create Folder Structure
ifelse(!dir.exists(file.path(paste0("../2_Output/Regulation/"))), 
       dir.create(file.path(paste0("../2_Output/Regulation/"))), 
       FALSE)
# Import the seurat object
TAA.combined <- readRDS(file = "../1_Input/NNT_Labelling_snRNA.rds")
# Import the differentially-expressed genes based on HFpEF vs Ctrl in N mice
DEGs <- openxlsx::read.xlsx("../2_Output/Cardiomyocyte/Cardiomyocyte_Venn.Diagram.xlsx", sheet = "N_ONLY_p05")$GeneName
# Filter seurat object by the differentially-expressed genes
TAA.combined_DEG <- subset(TAA.combined, features = DEGs)
# Run LIANA for cell-type enrichment
liana_TAA <- liana_wrap(TAA.combined_DEG, resource = "MouseConsensus")
liana_TAA <- liana_aggregate(liana_TAA)
liana_TAA %>%
  liana_dotplot(source_groups = c("Cardiomyocyte"),
                target_groups = c("Cardiomyocyte", "Fibroblast", "EC","Myeloid"),
                ntop = 20)
liana_trunc <- liana_TAA %>%  filter(cellphonedb.pvalue <= 0.05) # note that these pvals are already corrected
heat_freq(liana_trunc) # heatmap
pdf("../2_Output/Regulation/LIANA_Heatmap.pdf")
heat_freq(liana_trunc) # heatmap
dev.off()
colors <- c(Cardiomyocyte = "coral2", 
          EC = "wheat", 
          Fibroblast = "steelblue4", 
          Mural_Cell = "deepskyblue3",
          Lymphoid = "azure4",
          Myeloid = "goldenrod2", 
          Mast_Cell = "tan2",
          Neural_Cell = "darkcyan")
# Load necessary libraries
library(circlize)
library(dplyr)
library(ComplexHeatmap)
#############################################################################################################
TAA.combined_N.HFpEF <- subset(TAA.combined, subset = Background == "N" & Treatment == "HFpEF") #%>% subset(., features = DEGs)
liana_TAA <- liana_wrap(TAA.combined_N.HFpEF, resource = "MouseConsensus")
liana_TAA <- liana_aggregate(liana_TAA)
liana_trunc <- liana_TAA %>%  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected
pdf("../2_Output/Regulation/CellCell_DotPlot_N.Ctrl_CM.targets.pdf", height = 7, width = 7)
liana_TAA %>%
  liana_dotplot(source_groups = c("Cardiomyocyte"),
                target_groups = c("Cardiomyocyte", "Fibroblast", "EC", "Mural_Cell", "Myeloid", "Lymphoid", "Mast_Cell", "Neural_Cell"),
                ntop = 20) + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 0), axis.title.x = element_text(size = 0))
dev.off()
#############################################
# Load necessary libraries
library(dplyr)
library(circlize)
# Filter interactions for the ligand Fgf13 from the filtered liana_trunc object
ligand_of_interest <- "Edn1"
fgf13_targets <- liana_trunc %>%
    filter(ligand.complex == ligand_of_interest) %>%
    select(ligand.complex, receptor.complex, source, target) %>%
    distinct()  %>%
    group_by(target) %>%
    mutate(receptor_position = dense_rank(paste(ligand.complex, receptor.complex, target, sep = "-"))) %>%
    ungroup()

# Define unique cell types
source_types <- unique(c(fgf13_targets$target,fgf13_targets$source))
source_types <- factor(source_types)
# Define colors for each cell type
grid.col <- c(
    Cardiomyocyte = "coral2", 
    EC = "wheat", 
    Fibroblast = "steelblue4", 
    Mural_Cell = "deepskyblue3",
    Lymphoid = "azure4",
    Myeloid = "goldenrod2", 
    Mast_Cell = "tan2",
    Neural_Cell = "darkcyan"
)

# Map each cell type to the maximum number of receptors within it
receptor_counts <- fgf13_targets %>% select(target,receptor.complex) %>% distinct() %>%
    group_by(target) %>%
    summarise(count = n())
source_counts <- fgf13_targets %>%
    filter(!(source %in% target)) %>%
    select(source) %>%
    distinct() %>%
    group_by(source) %>%
    mutate(count = 1) %>% distinct() # Add 1 to count the source
# Join the calculated counts back into the original fgf13_targets data frame
fgf13_targets <- dplyr::inner_join(fgf13_targets, receptor_counts)
# Define sector widths based on receptor counts per cell type
sector_widths <- setNames(c(receptor_counts$count, source_counts$count), c(receptor_counts$target, source_counts$source))

### Receptor Targets
pdf(paste0("../2_Output/Regulation/", ligand_of_interest, "_targets.pdf"), height = 5, width = 6)
# Initialize circos plot with main cell type sectors
circos.clear()
circos.par(gap.degree = 5)

# Initialize main sectors for each cell type with custom widths
circos.initialize(factors = source_types, xlim = cbind(rep(0, length(sector_widths)), sector_widths))
# Add an outer track with cell type names on the border
circos.trackPlotRegion(factors = source_types, ylim = c(0, 1), bg.border = NA, track.height = 0.05, 
                       panel.fun = function(x, y) {
    cell_type <- CELL_META$sector.index
    circos.text(CELL_META$xcenter, 1, cell_type, facing = "bending.outside", niceFacing = TRUE, 
                adj = c(0.5, 1), cex = 0.8, col = "black")
})
# Draw colored boxes for each cell type sector and receptor sub-sectors
circos.trackPlotRegion(factors = source_types, ylim = c(0, 1), panel.fun = function(x, y) {
    source_types <- CELL_META$sector.index
    cell_color <- grid.col[source_types]
    receptors <- fgf13_targets %>% filter(target == source_types) %>% pull(receptor.complex)
    receptor_positions <- fgf13_targets %>% filter(target == source_types) %>% pull(receptor_position)
    
    # Draw cell type background
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = cell_color, border = "black")
    
    # Add receptor names evenly spaced within each cell type sector
    for (j in seq_along(receptors)) {
        pos <- CELL_META$xlim[1] + receptor_positions[j] - 1
        circos.text(pos + 0.5, 0.5, receptors[j], facing = "bending.inside", niceFacing = TRUE, adj = c(0.5, 1), cex = 0.5)
    }
}, track.height = 0.15, bg.border = NA)
###
unique_sources <- unique(fgf13_targets$source)
for (source in unique_sources) {
    source_links <- fgf13_targets %>% filter(source == !!source)
    
    apply(source_links, 1, function(row) {
        target <- row["target"]
        receptor_position <- as.numeric(row["receptor_position"])  # Ensure numeric type for calculation
        
        # Calculate the exact position within the target sector
        target_pos <- receptor_position - 0.5
        
        # Draw the link with specified color and line width
        circos.link(sector.index1 = source, point1 = 0-0.2, 
                    sector.index2 = target, point2 = target_pos, 
                    col = grid.col[source], lwd = 2)
    })
}
# Add a title to the plot
mtext(paste0(ligand_of_interest, " Receptor Interactions"), side = 3, line = -1, outer = TRUE, cex = 1, font = 2)
dev.off()
circos.clear()
########################
## Fgf1 targets
# Create the UMAP plot with Fgf1 expression
pdf(paste0("../2_Output/Regulation/Ligand_", ligand_of_interest, "_umap.pdf"), height = 4, width = 4)
FeaturePlot(TAA.combined, features = ligand_of_interest, reduction = "umap", cols = c("darkgray","white", "firebrick4"), alpha = 0.7) +
  ggtitle(paste0(ligand_of_interest," Expression")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
# Create composive Fgf1 target score using the downstream genes identified via cell-cell interactome
genes_vector <- unique(fgf13_targets$receptor.complex)
gene_list <- genes_vector[genes_vector %in% rownames(TAA.combined)] 
gene_expr <- GetAssayData(object = TAA.combined, slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
TAA.combined$Fgf1_score <- scale(sum_expr)
# plot composite targets score
# Load necessary libraries
pdf(paste0("../2_Output/Regulation/Target_score_", ligand_of_interest, "_umap.pdf"), height = 4, width = 4)
FeaturePlot(TAA.combined, features = "Fgf1_score", reduction = "umap") +
  scale_color_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick4", 
    midpoint = 0,   # Center the color scale at 0
    limits = c(min(TAA.combined[["Fgf1_score"]]), max(TAA.combined[["Fgf1_score"]]))  # Adjust limits to match data range
  ) +
  ggtitle(paste0(ligand_of_interest," Target Score")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
library(Nebulosa)
pdf(paste0("../2_Output/Regulation/Target_score_", ligand_of_interest, "_density.pdf"), height = 4, width = 4)
plot_density(TAA.combined, gene_list, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
dev.off()

#######
gene_of_interest <- c("Scn5a", "Fgfr1","Fgfr2", "Egfr")
# Generate individual violin plots for each gene and store them in a list
vln_plots <- lapply(gene_of_interest, function(gene) {
    VlnPlot(
        TAA.combined,
        features = gene,
        group.by = "CellType",
        split.by = "Treatment",
        adjust = 3,
        pt.size = 0
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),  # Remove Y-axis label for individual plots
        legend.position = "none"         # Remove individual legends
    ) +
    labs(title = gene)  + # Title for each gene plot 
    xlab(NULL)  # Remove the x-axis label
})

# Combine the plots using patchwork
combined_plot <- wrap_plots(vln_plots, ncol = 4) +
    plot_layout(guides = "collect") &  # Collect legends into one
    plot_annotation(
        title = "Expression of Fgf1/13 and downstream Targets",
        theme = theme(
            plot.title = element_text(hjust = 0.5)
        )
    ) &
    theme(legend.position = "right")   # Position the legend at the bottom

# Adjust the layout to add a Y-axis label closer to the plots
combined_plot_with_ylabel <- wrap_elements(
    grid::textGrob("Expression Level", rot = 90, gp = gpar(cex = 1.2))
) + 
plot_spacer() +  # Minimal spacing
combined_plot + 
plot_layout(widths = c(0.02, 0.005, 1))  # Further reduced widths for tighter spacing

# Display the final combined plot
pdf("../2_Output/Regulation/Violin_Fgf13.targets.pdf", width = 12, height = 3.3)
combined_plot_with_ylabel
dev.off()
#############################################
# Aggregate the interactions to create a single weight for each source-target pair
liana_data_aggregated <- liana_trunc %>%
  mutate(weight = -log10(aggregate_rank)) %>%  # Transform rank for better visualization
  group_by(source, target) %>%
  summarise(weight = sum(weight, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.infinite(weight) & weight > 0)  # Remove infinite or zero weights
# Convert the tibble into a data frame for `circlize::chordDiagram()`
liana_data_aggregated <- as.data.frame(liana_data_aggregated)
# Save the plot to a PDF
pdf(file = "../2_Output/Regulation/CellCell_Chord_N.HFpEF.pdf", height = 6, width = 9)
# Load necessary libraries
library(circlize)
library(ComplexHeatmap)
# Define colors for each cell type (adjust colors as needed)
grid.col <- c(Cardiomyocyte = "coral2", 
          EC = "wheat", 
          Fibroblast = "steelblue4", 
          Mural_Cell = "deepskyblue3",
          Lymphoid = "azure4",
          Myeloid = "goldenrod2", 
          Mast_Cell = "tan2",
          Neural_Cell = "darkcyan")
# Create the chord diagram with the consolidated weights
circlize::chordDiagram(
  x = liana_data_aggregated,
  grid.col = grid.col,
  transparency = 0.5,          # Adjust transparency for better visibility
  directional = 1,             # Show directional arrows
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow", # Different heights and arrows to indicate direction
  annotationTrack = "grid",    # Add annotation tracks for labels
  preAllocateTracks = list(track.height = 0.1)  # Adjust space for labels
)
# Add cell type labels with more space (adjust for aesthetics if needed)
circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  # Adjust the y position to increase distance from the segment
  circos.text(
    CELL_META$xcenter, CELL_META$ylim[1] + 0.3, CELL_META$sector.index, 
    facing = "bending",  # Align labels with their segment
    niceFacing = TRUE,   # Ensure labels are readable
    adj = c(0.5, 0.5)    # Center the text on the sector
  )
}, bg.border = NA)
# Create a legend using ComplexHeatmap
legend <- ComplexHeatmap::Legend(
  at = names(grid.col), 
  title = "Cell Type", 
  legend_gp = gpar(fill = grid.col)
)
# Draw the legend
ComplexHeatmap::draw(legend, x = unit(1, "npc") - unit(5, "mm"), just = "right")
# Close the PDF device
dev.off()
################
## Dotplot of target expression
library(scCustomize)
pdf(file = "../2_Output/Regulation/Dot.plot_top5markers.pdf", height = 10, width = 7)
Clustered_DotPlot(seurat_object = TAA.combined, split.by = "CellType", group.by = "Treatment",features = gene_list, assay = "RNA")
dev.off()

# Load required packages
library(Seurat)
library(ggplot2)

# Define the list of receptor genes
receptor_genes <- c("Fgfr1", "Fgfr2", "Egfr", "Scn5a")  # Replace with your receptor genes

# Generate the dot plot
# Generate the dot plot with adjustments
DotPlot(TAA.combined, features = receptor_genes, group.by = "CellType", split.by = "Treatment") + 
    # scale_color_gradientn(colors = colorRampPalette(c("dodgerblue4", "white", "goldenrod1"))(20)) +  # Custom color scale
    facet_wrap(~ "Background", ncol = 1) +  # Facet by Background variable
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Receptor Gene Expression Across Cell Types", x = "Receptor Genes", y = "Cell Types")

#########################################################################
# N - Ctrl
TAA.combined_N.Ctrl <- subset(TAA.combined, subset = Background == "N" & Treatment == "Ctrl")
liana_TAA <- liana_wrap(TAA.combined_N.Ctrl, resource = "MouseConsensus")
liana_TAA <- liana_aggregate(liana_TAA)
liana_trunc <- liana_TAA %>%  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected
### Interaction DotPlot
pdf("../2_Output/Regulation/CellCell_DotPlot_N.Ctrl_CM.targets.pdf", height = 7, width = 7)
liana_TAA %>%
  liana_dotplot(source_groups = c("Cardiomyocyte"),
                target_groups = c("Cardiomyocyte", "Fibroblast", "EC", "Mural_Cell", "Myeloid", "Lymphoid", "Mast_Cell", "Neural_Cell"),
                ntop = 20) + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 0), axis.title.x = element_text(size = 0))
dev.off()
# Aggregate the interactions to create a single weight for each source-target pair
liana_data_aggregated <- liana_trunc %>%
  mutate(weight = -log10(aggregate_rank)) %>%  # Transform rank for better visualization
  group_by(source, target) %>%
  summarise(weight = sum(weight, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.infinite(weight) & weight > 0)  # Remove infinite or zero weights
# Convert the tibble into a data frame for `circlize::chordDiagram()`
liana_data_aggregated <- as.data.frame(liana_data_aggregated)
# Save the plot to a PDF
pdf(file = "../2_Output/CellCell_Chord_N.Ctrl.pdf", height = 6, width = 9)
# Load necessary libraries
library(circlize)
library(ComplexHeatmap)
# Define colors for each cell type (adjust colors as needed)
grid.col <- c(Cardiomyocyte = "coral2", 
          EC = "wheat", 
          Fibroblast = "steelblue4", 
          Mural_Cell = "deepskyblue3",
          Lymphoid = "azure4",
          Myeloid = "goldenrod2", 
          Mast_Cell = "tan2",
          Neural_Cell = "darkcyan")
# Create the chord diagram with the consolidated weights
circlize::chordDiagram(
  x = liana_data_aggregated,
  grid.col = grid.col,
  transparency = 0.5,          # Adjust transparency for better visibility
  directional = 1,             # Show directional arrows
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow", # Different heights and arrows to indicate direction
  annotationTrack = "grid",    # Add annotation tracks for labels
  preAllocateTracks = list(track.height = 0.1)  # Adjust space for labels
)
# Add cell type labels with more space (adjust for aesthetics if needed)
circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  # Adjust the y position to increase distance from the segment
  circos.text(
    CELL_META$xcenter, CELL_META$ylim[1] + 0.3, CELL_META$sector.index, 
    facing = "bending",  # Align labels with their segment
    niceFacing = TRUE,   # Ensure labels are readable
    adj = c(0.5, 0.5)    # Center the text on the sector
  )
}, bg.border = NA)
# Create a legend using ComplexHeatmap
legend <- ComplexHeatmap::Legend(
  at = names(grid.col), 
  title = "Cell Type", 
  legend_gp = gpar(fill = grid.col)
)
# Draw the legend
ComplexHeatmap::draw(legend, x = unit(1, "npc") - unit(5, "mm"), just = "right")
# Close the PDF device
dev.off()
#####################################################################
# J - Ctrl
TAA.combined_N.Ctrl <- subset(TAA.combined, subset = Background == "J" & Treatment == "Ctrl")
liana_TAA <- liana_wrap(TAA.combined_N.Ctrl, resource = "MouseConsensus")
liana_TAA <- liana_aggregate(liana_TAA)
liana_trunc <- liana_TAA %>%  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected
## Interaction DotPlot
pdf("../2_Output/CellCell_DotPlot_J.Ctrl_CM.targets.pdf", height = 7, width = 7)
liana_TAA %>%
  liana_dotplot(source_groups = c("Cardiomyocyte"),
                target_groups = c("Cardiomyocyte", "Fibroblast", "EC", "Mural_Cell", "Myeloid", "Lymphoid", "Mast_Cell", "Neural_Cell"),
                ntop = 20) + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1), plot.title = element_text(size = 0), axis.title.x = element_text(size = 0))
dev.off()
# Aggregate the interactions to create a single weight for each source-target pair
liana_data_aggregated <- liana_trunc %>%
  mutate(weight = -log10(aggregate_rank)) %>%  # Transform rank for better visualization
  group_by(source, target) %>%
  summarise(weight = sum(weight, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.infinite(weight) & weight > 0)  # Remove infinite or zero weights
# Convert the tibble into a data frame for `circlize::chordDiagram()`
liana_data_aggregated <- as.data.frame(liana_data_aggregated)
# Save the plot to a PDF
pdf(file = "../2_Output/CellCell_Chord_J.Ctrl.pdf", height = 6, width = 9)
# Load necessary libraries
library(circlize)
library(ComplexHeatmap)
# Define colors for each cell type (adjust colors as needed)
grid.col <- c(Cardiomyocyte = "coral2", 
          EC = "wheat", 
          Fibroblast = "steelblue4", 
          Mural_Cell = "deepskyblue3",
          Lymphoid = "azure4",
          Myeloid = "goldenrod2", 
          Mast_Cell = "tan2",
          Neural_Cell = "darkcyan")
# Create the chord diagram with the consolidated weights
circlize::chordDiagram(
  x = liana_data_aggregated,
  grid.col = grid.col,
  transparency = 0.5,          # Adjust transparency for better visibility
  directional = 1,             # Show directional arrows
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow", # Different heights and arrows to indicate direction
  annotationTrack = "grid",    # Add annotation tracks for labels
  preAllocateTracks = list(track.height = 0.1)  # Adjust space for labels
)
# Add cell type labels with more space (adjust for aesthetics if needed)
circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  # Adjust the y position to increase distance from the segment
  circos.text(
    CELL_META$xcenter, CELL_META$ylim[1] + 0.3, CELL_META$sector.index, 
    facing = "bending",  # Align labels with their segment
    niceFacing = TRUE,   # Ensure labels are readable
    adj = c(0.5, 0.5)    # Center the text on the sector
  )
}, bg.border = NA)
# Create a legend using ComplexHeatmap
legend <- ComplexHeatmap::Legend(
  at = names(grid.col), 
  title = "Cell Type", 
  legend_gp = gpar(fill = grid.col)
)
# Draw the legend
ComplexHeatmap::draw(legend, x = unit(1, "npc") - unit(5, "mm"), just = "right")
# Close the PDF device
dev.off()
###############################################################
# J - HFpEF
TAA.combined_N.Ctrl <- subset(TAA.combined, subset = Background == "J" & Treatment == "HFpEF")
liana_TAA <- liana_wrap(TAA.combined_N.Ctrl, resource = "MouseConsensus")
liana_TAA <- liana_aggregate(liana_TAA)
liana_trunc <- liana_TAA %>%  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected
# Interaction DotPLot
pdf("../2_Output/CellCell_DotPlot_J.HFpEF_CM.targets.pdf", height = 6, width = 7)
liana_TAA %>%
  liana_dotplot(source_groups = c("Cardiomyocyte"),
                target_groups = c("Cardiomyocyte", "Fibroblast", "EC", "Mural_Cell", "Myeloid", "Lymphoid", "Mast_Cell", "Neural_Cell"),
                ntop = 20) + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face = "bold"), plot.title = element_text(size = 0), axis.title.x = element_text(size = 0))
dev.off()
# Aggregate the interactions to create a single weight for each source-target pair
liana_data_aggregated <- liana_trunc %>%
  mutate(weight = -log10(aggregate_rank)) %>%  # Transform rank for better visualization
  group_by(source, target) %>%
  summarise(weight = sum(weight, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.infinite(weight) & weight > 0)  # Remove infinite or zero weights
# Convert the tibble into a data frame for `circlize::chordDiagram()`
liana_data_aggregated <- as.data.frame(liana_data_aggregated)
# Save the plot to a PDF
pdf(file = "../2_Output/CellCell_Chord_J.HFpEF.pdf", height = 6, width = 9)
# Load necessary libraries
library(circlize)
library(ComplexHeatmap)
# Define colors for each cell type (adjust colors as needed)
grid.col <- c(Cardiomyocyte = "coral2", 
          EC = "wheat", 
          Fibroblast = "steelblue4", 
          Mural_Cell = "deepskyblue3",
          Lymphoid = "azure4",
          Myeloid = "goldenrod2", 
          Mast_Cell = "tan2",
          Neural_Cell = "darkcyan")
# Create the chord diagram with the consolidated weights
circlize::chordDiagram(
  x = liana_data_aggregated,
  grid.col = grid.col,
  transparency = 0.5,          # Adjust transparency for better visibility
  directional = 1,             # Show directional arrows
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow", # Different heights and arrows to indicate direction
  annotationTrack = "grid",    # Add annotation tracks for labels
  preAllocateTracks = list(track.height = 0.1)  # Adjust space for labels
)
# Add cell type labels with more space (adjust for aesthetics if needed)
circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  # Adjust the y position to increase distance from the segment
  circos.text(
    CELL_META$xcenter, CELL_META$ylim[1] + 0.3, CELL_META$sector.index, 
    facing = "bending",  # Align labels with their segment
    niceFacing = TRUE,   # Ensure labels are readable
    adj = c(0.5, 0.5)    # Center the text on the sector
  )
}, bg.border = NA)
# Create a legend using ComplexHeatmap
legend <- ComplexHeatmap::Legend(
  at = names(grid.col), 
  title = "Cell Type", 
  legend_gp = gpar(fill = grid.col)
)
# Draw the legend
ComplexHeatmap::draw(legend, x = unit(1, "npc") - unit(5, "mm"), just = "right")
# Close the PDF device
dev.off()
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
library(devtools)
library(DT)

# Capture the session information
sinfo <- session_info()

# Create an interactive table of the packages data
datatable(
  sinfo$packages,  # Display session information about packages
  options = list(
    pageLength = 10,  # Show 10 entries per page
    autoWidth = TRUE,
    searching = TRUE,
    ordering = TRUE,
    dom = 'lfrtip'  # Length, filter, pagination controls
  ),
  rownames = FALSE
)
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-39b76fb009bac3deb2d8" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-39b76fb009bac3deb2d8">{"x":{"filter":"none","vertical":false,"data":[["bslib","cachem","cli","colorspace","devtools","digest","DT","ellipsis","evaluate","fastmap","fs","glue","htmltools","htmlwidgets","httpuv","jquerylib","jsonlite","kableExtra","knitr","later","lifecycle","magrittr","memoise","mime","miniUI","munsell","pkgbuild","pkgload","profvis","promises","purrr","R6","Rcpp","remotes","rlang","rmarkdown","rstudioapi","sass","scales","sessioninfo","shiny","stringi","stringr","svglite","systemfonts","urlchecker","usethis","vctrs","viridisLite","xfun","xml2","xtable","yaml"],["0.8.0","1.1.0","3.6.3","2.1.1","2.4.5","0.6.37","0.33","0.3.2","1.0.0","1.2.0","1.6.4","1.8.0","0.5.8.1","1.6.4","1.6.15","0.1.4","1.8.9","1.4.0","1.48","1.3.2","1.0.4","2.0.3","2.0.1","0.12","0.1.1.1","0.5.1","1.4.4","1.4.0","0.4.0","1.3.0","1.0.2","2.5.1","1.0.13","2.5.0","1.1.4","2.28","0.16.0","0.4.9","1.3.0","1.2.2","1.9.1","1.8.4","1.5.1","2.1.3","1.1.0","1.0.1","3.0.0","0.6.5","0.4.2","0.47","1.3.6","1.8.4","2.3.10"],["0.8.0","1.1.0","3.6.3","2.1-1","2.4.5","0.6.37","0.33","0.3.2","1.0.0","1.2.0","1.6.4","1.8.0","0.5.8.1","1.6.4","1.6.15","0.1.4","1.8.9","1.4.0","1.48","1.3.2","1.0.4","2.0.3","2.0.1","0.12","0.1.1.1","0.5.1","1.4.4","1.4.0","0.4.0","1.3.0","1.0.2","2.5.1","1.0.13","2.5.0","1.1.4","2.28","0.16.0","0.4.9","1.3.0","1.2.2","1.9.1","1.8.4","1.5.1","2.1.3","1.1.0","1.0.1","3.0.0","0.6.5","0.4.2","0.47","1.3.6","1.8-4","2.3.10"],["/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bslib","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cachem","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cli","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/colorspace","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/devtools","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/digest","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DT","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ellipsis","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/evaluate","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmap","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fs","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/glue","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmltools","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmlwidgets","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httpuv","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jquerylib","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jsonlite","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/kableExtra","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/knitr","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/later","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lifecycle","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/magrittr","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/memoise","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mime","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/miniUI","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/munsell","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgbuild","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgload","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/profvis","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/promises","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/purrr","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R6","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rcpp","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/remotes","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rlang","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rmarkdown","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rstudioapi","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sass","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scales","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sessioninfo","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shiny","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringi","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringr","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/svglite","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/systemfonts","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/urlchecker","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/usethis","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vctrs","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/viridisLite","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xfun","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xml2","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xtable","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yaml"],["/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bslib","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cachem","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cli","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/colorspace","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/devtools","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/digest","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DT","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ellipsis","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/evaluate","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmap","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fs","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/glue","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmltools","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmlwidgets","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httpuv","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jquerylib","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jsonlite","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/kableExtra","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/knitr","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/later","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lifecycle","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/magrittr","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/memoise","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mime","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/miniUI","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/munsell","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgbuild","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgload","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/profvis","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/promises","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/purrr","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R6","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rcpp","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/remotes","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rlang","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rmarkdown","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rstudioapi","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sass","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scales","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sessioninfo","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shiny","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringi","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringr","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/svglite","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/systemfonts","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/urlchecker","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/usethis","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vctrs","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/viridisLite","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xfun","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xml2","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xtable","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yaml"],[false,false,false,false,true,false,true,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,false,false],[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false],["2024-07-29","2024-05-16","2024-06-21","2024-07-26","2022-10-11","2024-08-19","2024-04-04","2021-04-29","2024-09-17","2024-05-15","2024-04-25","2024-09-30","2024-04-04","2023-12-06","2024-03-26","2021-04-26","2024-09-20","2024-01-24","2024-07-07","2023-12-06","2023-11-07","2022-03-30","2021-11-26","2021-09-28","2018-05-18","2024-04-01","2024-03-17","2024-06-28","2024-09-20","2024-04-05","2023-08-10","2021-08-19","2024-07-17","2024-03-17","2024-06-04","2024-08-17","2024-03-24","2024-03-15","2023-11-28","2021-12-06","2024-08-01","2024-05-06","2023-11-14","2023-12-08","2024-05-15","2021-11-30","2024-07-29","2023-12-01","2023-05-02","2024-08-17","2023-12-04","2019-04-21","2024-07-26"],["CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.1)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.1)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.1)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.1)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.1)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)","CRAN (R 4.4.0)"],[null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],["/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library","/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>package<\/th>\n      <th>ondiskversion<\/th>\n      <th>loadedversion<\/th>\n      <th>path<\/th>\n      <th>loadedpath<\/th>\n      <th>attached<\/th>\n      <th>is_base<\/th>\n      <th>date<\/th>\n      <th>source<\/th>\n      <th>md5ok<\/th>\n      <th>library<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":10,"autoWidth":true,"searching":true,"ordering":true,"dom":"lfrtip","columnDefs":[{"name":"package","targets":0},{"name":"ondiskversion","targets":1},{"name":"loadedversion","targets":2},{"name":"path","targets":3},{"name":"loadedpath","targets":4},{"name":"attached","targets":5},{"name":"is_base","targets":6},{"name":"date","targets":7},{"name":"source","targets":8},{"name":"md5ok","targets":9},{"name":"library","targets":10}],"order":[],"orderClasses":false},"selection":{"mode":"multiple","selected":null,"target":"row","selectable":null}},"evals":[],"jsHooks":[]}</script>
```
