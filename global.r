library(shiny)
library(rhdf5)
library(plotrix)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(plotly)
library(tidyr)
library(shiny)
library(rhdf5)
library(Seurat)
library(shinydashboard)
library(ggthemes)
library(shinycssloaders)


source("functions.R")
gene_de.file <- "data/AllDEtable.txt"
expression.file <- ifelse(file.exists("data/AgingData.h5"),"data/AgingData.h5", "data/scaledData.h5")
cell_info.file <- "data/metadata.txt"
tSNE.file <- "data/tSNE_coord.txt"
genes.file <- "data/genes.txt"
solubility.file <- "data/Solubility.RData"
protein_bulk.file <- "data/BulkProtein.RData"
markers.file <- "data/AllMarkersCelltypes.txt"


load(solubility.file)
load(protein_bulk.file)
markers_table <- fread(markers.file)
gene_de_table <- read.delim(gene_de.file, check.names = F)
gene_de_table <- as.data.table(gene_de_table, keep.rownames = T)
colnames(gene_de_table) <- gsub("\\s", "__", colnames(gene_de_table))
cell_info <- as.data.table(read.table(cell_info.file), keep.rownames = T)
tsne_coord <- read.table(tSNE.file, sep = "\t")
genes <- scan(what = character(), genes.file, sep = "\n")
cell_types <- unique(cell_info$celltype) %>% as.character()
cell_types <- cell_types[which(!is.na(cell_types))]

print(paste("tSNE file:",tSNE.file))
print(paste("Genes file:",genes.file))
print(paste("Gene DE file:",gene_de.file))
print(paste("Cell info file:",cell_info.file))
print(paste("Expression file:",expression.file))

