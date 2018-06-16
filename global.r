library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(ggthemes)
library(rhdf5)
library(data.table)
library(plotrix)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(plotly)
library(tidyr)
library(shiny)
library(Seurat)
library(parallel)
library(future)
library(shinyjs)


source("functions.R")
gene_de.file <- "data/AllDEtable.txt"
expression.file <- ifelse(file.exists("data/AgingData.h5"),"data/AgingData.h5", "data/scaledData.h5")
cell_info.file <- "data/metadata.txt"
tSNE.file <- "data/tSNE_coord.txt"
genes.file <- "data/genes.txt"
solubility.file <- "data/Solubility.RData"
protein_bulk.file <- "data/BulkProtein.RData"
markers.file <- "data/AllMarkersCelltypes.txt"
enrichment.file <- "data/enrichment_table.csv"

# Load solubility file
load(solubility.file)
# Load protein bulk file
# Load markers table
load(protein_bulk.file)
markers_table <- fread(markers.file)
# Load gene_de_table
gene_de_table <- read.delim(gene_de.file, check.names = F)
gene_de_table <- as.data.table(gene_de_table, keep.rownames = T)
colnames(gene_de_table) <- gsub("\\s", "__", colnames(gene_de_table))
# Load cell_info file
cell_info <- as.data.table(read.table(cell_info.file), keep.rownames = T)
# Extract cell types from cell_info
cell_types <- unique(cell_info$celltype) %>% as.character()
cell_types <- cell_types[which(!is.na(cell_types))]
bad <-
  c("red_blood_cells",
    "Gamma-Delta_T_cells",
    "low_quality_cells")
cell_types <- cell_types[-which(cell_types %in% bad)]
cell_types <- cell_types[order(cell_types)] 
# Load tsne file
tsne_coord <- read.table(tSNE.file, sep = "\t")
# Load gene file
genes <- scan(what = character(), genes.file, sep = "\n")
# Read encrichment table
enrichment_table <- fread(enrichment.file)
enrichment_table[ , Column := tstrsplit(Column, "\\s", keep=1)]
setnames(enrichment_table, c("Column"), c("Cell type"))
enrichment_table <- enrichment_table[, c("Cell type", "Type", "Name", "Score", "Benj. Hoch. FDR"), with=F]
enrichment_types <- c("All", (enrichment_table$Type %>% unique))

print(paste("tSNE file:",tSNE.file))
print(paste("Genes file:",genes.file))
print(paste("Gene DE file:",gene_de.file))
print(paste("Cell info file:",cell_info.file))
print(paste("Expression file:",expression.file))