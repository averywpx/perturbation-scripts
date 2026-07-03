install.packages(c("shiny", "plotly", "miniUI"), repos="https://cloud.r-project.org")

install.packages("BiocManager")
install.packages("Seurat")
install.packages("tidyverse")
install.packages("devtools")

install.packages("remotes")
install.packages("hdf5r", repos="https://cloud.r-project.org")
remotes::install_github("mojaveazure/seurat-disk")

install.packages("usethis")
install.packages("gitcreds")

usethis::create_github_token()
gitcreds::gitcreds_set()

remotes::install_github('cellgeni/sceasy')

install.packages('mixtools')
install.packages('devtools')
devtools::install_github('immunogenomics/presto')


library(shiny)
library(BiocManager)
library(dplyr)
library(Seurat)
library(patchwork)

library(Seurat)
library(SeuratData)
library(SeuratDisk)

library(sceasy)
library(mixtools)
library(scales)
library(ggplot2)
library(reshape2)

# Setup custom theme for plotting.
custom_theme <- theme(
  plot.title = element_text(size=16, hjust = 0.5), 
  legend.key.size = unit(0.7, "cm"), 
  legend.text = element_text(size = 14))



sceasy::convertFormat(
  "/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_sub_ntc.h5ad", 
  from="anndata", 
  to="seurat", 
  outFile="/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_sub_nt.rds")

# Load the newly created object
x_atlas_obj <- readRDS("/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_sub_nt.rds")

x_atlas_obj <- NormalizeData(x_atlas_obj)
x_atlas_obj <- FindVariableFeatures(x_atlas_obj)
x_atlas_obj <- ScaleData(x_atlas_obj)
x_atlas_obj <- RunPCA(x_atlas_obj)

x_atlas_obj <- CalcPerturbSig(
  object = x_atlas_obj,
  assay = "RNA",
  slot = "data",
  gd.class = "gene_target",
  nt.cell.class = "Non-Targeting",
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = "PRTB"
)

x_atlas_obj <- RunMixscape(
  object = x_atlas_obj,
  assay = "PRTB",
  slot = "data",
  labels = "gene_target",
  nt.class.name = "Non-Targeting",
  new.class.name = "mixscape_class",
  de.assay = "RNA",
  prtb.type = "KO"
)


x_atlas_obj <- RunMixscape(
  object = x_atlas_obj,
  assay = "RNA",
  slot = "scale.data",
  labels = "gene_target",
  nt.class.name = "Non-Targeting",
  new.class.name = "mixscape_class",
  de.assay = "RNA",
  prtb.type = "KO"
)

table(x_atlas_obj$gene_target, x_atlas_obj$mixscape_class.global)
prop.table(table(x_atlas_obj$gene_target, x_atlas_obj$mixscape_class.global), margin = 1)

tab <- prop.table(table(x_atlas_obj$gene_target, x_atlas_obj$mixscape_class.global), margin = 1)
tab[, "KO"]

#========================================================================================================

# Read Replogle data
# raw hepg2 cell line
sceasy::convertFormat(
  "/data/storage/Replogle/raw_data/GSE264667_hepg2_raw_singlecell_01_sub.h5ad", 
  from="anndata", 
  to="seurat", 
  outFile="/data/storage/Replogle/raw_data/GSE264667_hepg2_raw_singlecell_01_sub.rds")

# Load the newly created object
raw_hepg2 <- readRDS("/data/storage/Replogle/raw_data/GSE264667_hepg2_raw_singlecell_01_sub.rds")

raw_hepg2 <- NormalizeData(raw_hepg2)
raw_hepg2 <- FindVariableFeatures(raw_hepg2)
raw_hepg2 <- ScaleData(raw_hepg2)
raw_hepg2 <- RunPCA(raw_hepg2)
raw_hepg2 <- RunUMAP(object = raw_hepg2, dims = 1:15)

rp1 <- DimPlot(
  object = raw_hepg2, 
  group.by = 'gene', 
  label = F, 
  pt.size = 0.7, 
  reduction = "umap", cols = "Dark2", repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("genes") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  custom_theme

rp1

rp2 <- DimPlot(
  object = raw_hepg2, 
  group.by = 'perturbed', 
  label = F, 
  pt.size = 0.7, 
  reduction = "umap", cols = "Dark2", repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("genes") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  custom_theme

rp2

rp3 <- DimPlot(
  object = raw_hepg2,
  group.by = 'perturbed',
  reduction = 'umap', 
  split.by = "perturbed", 
  ncol = 1, 
  pt.size = 0.7, 
  cols = c("grey39","goldenrod3")) +
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

rp3

raw_hepg2 <- CalcPerturbSig(
  object = raw_hepg2,
  assay = "RNA",
  slot = "data",
  gd.class = "gene",
  nt.cell.class = "non-targeting",
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = "PRTB"
)

raw_hepg2 <- RunMixscape(
  object = raw_hepg2,
  assay = "PRTB",
  slot = "data",
  labels = "gene",
  nt.class.name = "non-targeting",
  new.class.name = "mixscape_class",
  de.assay = "RNA",
  prtb.type = "KO"
)


table(raw_hepg2$gene, raw_hepg2$mixscape_class.global)
prop.table(table(raw_hepg2$gene, raw_hepg2$mixscape_class.global), margin = 1)

rtab <- prop.table(table(raw_hepg2$gene, raw_hepg2$mixscape_class.global), margin = 1)
rtab[, "KO"]

rp4 <- DimPlot(
  object = raw_hepg2, 
  group.by = 'gene', 
  label = F, 
  pt.size = 1.5, 
  reduction = "prtbumap", cols = "Dark2", repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("genes") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  custom_theme

rp4


# Remove non-perturbed cells and run LDA to reduce the dimensionality of the data.
Idents(raw_hepg2) <- "mixscape_class.global"
rawsub <- subset(raw_hepg2, idents = c("KO", "non-targeting"))

# Run LDA.
rawsub <- MixscapeLDA(
  object = rawsub, 
  assay = "RNA", 
  pc.assay = "PRTB", 
  labels = "gene", 
  nt.label = "non-targeting", 
  npcs = 10, 
  logfc.threshold = 0.25, 
  verbose = F)

# Use LDA results to run UMAP and visualize cells on 2-D. 
# Here, we note that the number of the dimensions to be used is equal to the number of 
# labels minus one (to account for NT cells).
rawsub <- RunUMAP(
  object = rawsub,
  dims = 1:7,
  reduction = 'lda',
  reduction.key = 'ldaumap',
  reduction.name = 'ldaumap')

# Visualize UMAP clustering results.
Idents(rawsub) <- "mixscape_class"
rawsub$mixscape_class <- as.factor(rawsub$mixscape_class)

## Set colors for each perturbation.
#col = setNames(object = hue_pal()(12),nm = levels(sub$mixscape_class))
#names(col) <- c(names(col)[1:7], "non-targeting", names(col)[9:12])
#col[8] <- "grey39"

genes <- levels(factor(rawsub$gene))
col <- setNames(hue_pal()(length(genes)), genes)

rmp2 <- DimPlot(object = rawsub, 
              reduction = "ldaumap", 
              repel = T, 
              label.size = 5, 
              label = T,
              group.by = "gene",
              cols = col) + 
  scale_color_manual(values = col, drop = FALSE) +
  theme(legend.position = "right")
#NoLegend()

rmp2

# Run PCA to reduce the dimensionality of the data.
raw_hepg2 <- RunPCA(object = raw_hepg2, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
raw_hepg2 <- RunUMAP(
  object = raw_hepg2, 
  dims = 1:15, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')

rmp3 <- DimPlot(
  object = raw_hepg2,
  reduction = "prtbumap",
  group.by = "mixscape_class.global",
  pt.size = 1.5,
  cols = c("KO" = "goldenrod3",
           "NP" = "grey39",
           "non-targeting" = "steelblue")
)

rmp3

# =============================================
# processed hepg2 cell line
sceasy::convertFormat(
  "/data/storage/Replogle/processed_data/replogle_hepg2_sub.h5ad", 
  from="anndata", 
  to="seurat", 
  outFile="/data/storage/Replogle/processed_data/replogle_hepg2_sub.rds")

# Load the newly created object
processed_hepg2 <- readRDS("/data/storage/Replogle/processed_data/replogle_hepg2_sub.rds")

processed_hepg2 <- NormalizeData(processed_hepg2)
processed_hepg2 <- FindVariableFeatures(processed_hepg2)
processed_hepg2 <- ScaleData(processed_hepg2)
processed_hepg2 <- RunPCA(processed_hepg2)
processed_hepg2 <- RunUMAP(object = processed_hepg2, dims = 1:15)

pp1 <- DimPlot(
  object = processed_hepg2, 
  group.by = 'gene', 
  label = F, 
  pt.size = 1.5, 
  reduction = "umap", cols = "Dark2", repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("genes") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  custom_theme

pp1

p5 <- DimPlot(
  object = processed_hepg2, 
  group.by = 'gene', 
  label = F, 
  pt.size = 1.0, 
  reduction = "umap", cols = "Dark2", repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("genes") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  custom_theme

p5

pp2 <- DimPlot(
  object = processed_hepg2, 
  group.by = 'perturbed', 
  label = F, 
  pt.size = 1.5, 
  reduction = "umap", cols = "Dark2", repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("genes") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  custom_theme

pp2

pp3 <- DimPlot(
  object = processed_hepg2, 
  group.by = 'perturbed', 
  pt.size = 0.7, 
  reduction = "umap", 
  split.by = "perturbed", 
  ncol = 1, 
  cols = c("grey39","goldenrod3")) + 
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme
pp3

processed_hepg2 <- CalcPerturbSig(
  object = processed_hepg2,
  assay = "RNA",
  slot = "data",
  gd.class = "gene",
  nt.cell.class = "non-targeting",
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = "PRTB"
)

# Prepare PRTB assay for dimensionality reduction: 
# Normalize data, find variable features and center data.
DefaultAssay(object = processed_hepg2) <- 'PRTB'

# Use variable features from RNA assay.
VariableFeatures(object = processed_hepg2) <- VariableFeatures(object = processed_hepg2[["RNA"]])
processed_hepg2 <- ScaleData(object = processed_hepg2, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
processed_hepg2 <- RunPCA(object = processed_hepg2, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
processed_hepg2 <- RunUMAP(
  object = processed_hepg2, 
  dims = 1:15, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')

q3 <- DimPlot(
  object = processed_hepg2,
  group.by = 'perturbed',
  reduction = 'prtbumap', 
  split.by = "perturbed", 
  ncol = 1, 
  pt.size = 0.7, 
  cols = c("grey39","goldenrod3")) +
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

q3

processed_hepg2 <- RunMixscape(
  object = processed_hepg2,
  assay = "PRTB",
  slot = "data",
  labels = "gene",
  nt.class.name = "non-targeting",
  new.class.name = "mixscape_class",
  de.assay = "RNA",
  prtb.type = "KO"
)

table(processed_hepg2$gene, processed_hepg2$mixscape_class.global)
prop.table(table(processed_hepg2$gene, processed_hepg2$mixscape_class.global), margin = 1)

ptab <- prop.table(table(processed_hepg2$gene, processed_hepg2$mixscape_class.global), margin = 1)
ptab[, "KO"]

# Run UMAP to visualize clustering in 2-D.
processed_hepg2 <- RunUMAP(
  object = processed_hepg2, 
  dims = 1:15, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')

pp4 <- DimPlot(
  object = processed_hepg2, 
  group.by = 'gene', 
  label = F, 
  pt.size = 1.5, 
  reduction = "prtbumap", cols = "Dark2", repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("genes") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  custom_theme

pp4

q4 <- DimPlot(
  object = processed_hepg2,
  group.by = 'perturbed',
  reduction = 'prtbumap', 
  split.by = "perturbed", 
  ncol = 1, 
  pt.size = 0.7, 
  cols = c("grey39","goldenrod3")) +
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

q4

q5 <- DimPlot(
  object = processed_hepg2,
  group.by = 'gene',
  reduction = 'prtbumap', 
  split.by = "gene", 
  ncol = 1, 
  pt.size = 0.7, 
  cols = c("grey39","goldenrod3")) +
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

q5

q5 <- DimPlot(
  object = processed_hepg2, 
  group.by = 'gene', 
  label = F, 
  pt.size = 1.5, 
  reduction = "prtbumap", cols = "Dark2", repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("genes") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  custom_theme

q5

q6 <- DimPlot(
  object = processed_hepg2,
  reduction = "prtbumap",
  group.by = "mixscape_class.global",
  pt.size = 1.5,
  cols = c("KO" = "goldenrod3",
           "NP" = "grey39",
           "non-targeting" = "steelblue")
)

q6


# Convert("/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_st14_nt.h5ad", dest="h5seurat", overwirte=TRUE)

# obj <- LoadH5Seurate("/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_st14_nt.h5seurat")

# ============== UMAP ======================
# Remove non-perturbed cells and run LDA to reduce the dimensionality of the data.
Idents(processed_hepg2) <- "mixscape_class.global"
sub <- subset(processed_hepg2, idents = c("KO", "non-targeting"))

# Run LDA.
sub <- MixscapeLDA(
  object = sub, 
  assay = "RNA", 
  pc.assay = "PRTB", 
  labels = "gene", 
  nt.label = "non-targeting", 
  npcs = 10, 
  logfc.threshold = 0.25, 
  verbose = F)

# Use LDA results to run UMAP and visualize cells on 2-D. 
# Here, we note that the number of the dimensions to be used is equal to the number of 
# labels minus one (to account for NT cells).
sub <- RunUMAP(
  object = sub,
  dims = 1:7,
  reduction = 'lda',
  reduction.key = 'ldaumap',
  reduction.name = 'ldaumap')

# Visualize UMAP clustering results.
Idents(sub) <- "mixscape_class"
sub$mixscape_class <- as.factor(sub$mixscape_class)

## Set colors for each perturbation.
#col = setNames(object = hue_pal()(12),nm = levels(sub$mixscape_class))
#names(col) <- c(names(col)[1:7], "non-targeting", names(col)[9:12])
#col[8] <- "grey39"

genes <- levels(factor(sub$gene))
col <- setNames(hue_pal()(length(genes)), genes)

pmp3 <- DimPlot(object = sub, 
             reduction = "ldaumap", 
             repel = T, 
             label.size = 5, 
             label = T,
             group.by = "gene",
             cols = col) + 
  scale_color_manual(values = col, drop = FALSE) +
  theme(legend.position = "right")
  #NoLegend()

pmp3

p2 <- p+ 
  scale_color_manual(values=col, drop=FALSE) + 
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme
p2


# ============== Save results ============

ko_obj <- subset(
  processed_hepg2,
  subset = grepl(" KO$", mixscape_class)
)

np_obj <- subset(
  processed_hepg2,
  subset = grepl(" NP$", mixscape_class)
)

table(ko_obj$mixscape_class)
table(np_obj$mixscape_class)

writeLines(Cells(ko_obj), "/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/processed_hepg2_ko_cells.txt")
writeLines(Cells(np_obj), "/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/processed_hepg2_np_cells.txt")

# ==============

ko_obj <- subset(
  raw_hepg2,
  subset = grepl(" KO$", mixscape_class)
)

np_obj <- subset(
  raw_hepg2,
  subset = grepl(" NP$", mixscape_class)
)

table(ko_obj$mixscape_class)
table(np_obj$mixscape_class)

writeLines(Cells(ko_obj), "/data/storage/Replogle/raw_data/GSE264667_hepg2_raw_singlecell_01_sub_ko_cells.rds")
writeLines(Cells(np_obj), "/data/storage/Replogle/raw_data/GSE264667_hepg2_raw_singlecell_01_sub_np_cells.rds")















