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

library(SeuratDisk)

library(sceasy)
library(mixtools)


sceasy::convertFormat(
  "/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_st14_nt.h5ad", 
  from="anndata", 
  to="seurat", 
  outFile="/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_st14_nt.rds")

# Load the newly created object
obj <- readRDS("/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_st14_nt.rds")


# Read Replogle data
# raw hepg2 cell line
sceasy::convertFormat(
  "/data/storage/Replogle/raw_data/GSE264667_hepg2_raw_singlecell_01.h5ad", 
  from="anndata", 
  to="seurat", 
  outFile="/data/storage/Replogle/raw_data/GSE264667_hepg2_raw_singlecell_01.rds")

# Load the newly created object
raw_hepg2 <- readRDS("/data/storage/Replogle/raw_data/GSE264667_hepg2_raw_singlecell_01.rds")


# processed hepg2 cell line
sceasy::convertFormat(
  "/data/storage/Replogle/processed_data/replogle_hepg2.h5ad", 
  from="anndata", 
  to="seurat", 
  outFile="/data/storage/Replogle/processed_data/replogle_hepg2.rds")

# Load the newly created object
processed_hepg2 <- readRDS("/data/storage/Replogle/processed_data/replogle_hepg2.rds")



obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)

obj <- RunMixscape(
  object = obj,
  assay = "RNA",
  slot = "scale.data",
  labels = "gene_target",
  nt.class.name = "Non-Targeting",
  new.class.name = "mixscape_class",
  de.assay = "RNA",
  prtb.type = "KO"
)

table(obj$gene_target, obj$mixscape_class.global)
prop.table(table(obj$gene_target, obj$mixscape_class.global), margin = 1)

tab <- prop.table(table(obj$gene_target, obj$mixscape_class.global), margin = 1)
tab[, "KO"]

# Convert("/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_st14_nt.h5ad", dest="h5seurat", overwirte=TRUE)

# obj <- LoadH5Seurate("/data/storage/X-Atlas/HCT116/h5ad/intermediate_data/xatlas_st14_nt.h5seurat")















