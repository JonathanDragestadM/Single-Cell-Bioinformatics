library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(fgsea)
library(Signac)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)

# here()

# Load Genome Annotations
addArchRGenome("hg38")

# Set Threads to be used
# threads <- 16
addArchRThreads()

metadata <- read.table("../project/data/hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
metadata <- metadata |> filter(ATAC_Sample %in% c('A001-C-007-D', 'A015-C-001-D', 'CRC-1-8810-D', 'CRC-2-15564-D', 'CRC-3-11773-D',
                                                  'B001-A-401-D', 'B001-A-406-D', 'B001-A-501-D', 'B004-A-008-D', 'B004-A-204-D'))
pre_defined_doublets <- TRUE
# cells_to_include <- read.table("./cells_in_initial_clustering.txt")

inputFiles1 <- c(
  "A001-C-007-D" = "../project/data/GSM6058771_A001-C-007-D_20200817_fragments.tsv.gz",
  "A015-C-001-D" = "../project/data/GSM6058824_A015-C-001-D_20200817_fragments.tsv.gz",
  "B001-A-401-D" = "../project/data/GSM6058852_B001-A-401-D_fragments.tsv.gz",
  "B001-A-406-D" = "../project/data/GSM6058853_B001-A-406-D_fragments.tsv.gz"
)

inputFiles2 <- c(
  "B001-A-501-D" = "../project/data/GSM6058854_B001-A-501-D_fragments.tsv.gz",
  "B004-A-008-D" = "../project/data/GSM6058857_B004-A-008-D_20200817_fragments.tsv.gz",
  "B004-A-204-D" = "../project/data/GSM6058858_B004-A-204-D_20200702_fragments.tsv.gz"
)

inputFiles3 <- c(
  "CRC-1-8810-D" = "../project/data/GSM6058842_CRC-1-8810-D_20200917_fragments.tsv.gz",
  "CRC-2-15564-D" = "../projectdata/GSM6058843_CRC-2-15564-D_20200917_fragments.tsv.gz",
  "CRC-3-11773-D" = "../project/data/GSM6058844_CRC-3-11773-D_20200917_fragments.tsv.gz"
)

ArrowFiles1 <- createArrowFiles(
  inputFiles = inputFiles1,
  sampleNames = names(inputFiles1), minFrags = 1500,
)

ArrowFiles2 <- createArrowFiles(
  inputFiles = inputFiles2,
  sampleNames = names(inputFiles2), minFrags = 1500,
)

ArrowFiles3 <- createArrowFiles(
  inputFiles = inputFiles3,
  sampleNames = names(inputFiles3), minFrags = 1500,
)

doubScores <- addDoubletScores(ArrowFiles1, k = 10, knnMethod = "UMAP", LSIMethod = 1)
doubScores <- addDoubletScores(ArrowFiles2, k = 10, knnMethod = "UMAP", LSIMethod = 1)
doubScores <- addDoubletScores(ArrowFiles3, k = 10, knnMethod = "UMAP", LSIMethod = 1)

ArrowFiles <- c("A001-C-007-D.arrow",
                "A015-C-001-D.arrow",
                "CRC-1-8810-D.arrow",
                "CRC-2-15564-D.arrow",
                "CRC-3-11773-D.arrow",
                "B001-A-401-D.arrow",
                "B001-A-406-D.arrow",
                "B001-A-501-D.arrow",
                "B004-A-008-D.arrow",
                "B004-A-204-D.arrow")

# #Create ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "scATACseq_data"
)

for (j in 2:dim(metadata)[2]){
  # initialize list
  cellsNamesToAdd <- c()
  annotationToAdd <- c()
  for (i in 1:dim(metadata)[1]){
    idxSample <- BiocGenerics::which(getCellColData(proj, "Sample") %in% metadata[i,"ATAC_Sample"])
    cellsSample <- proj$cellNames[idxSample[["Sample"]]]
    cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
    annotationToAdd <- append(annotationToAdd, rep(metadata[i,j], length(cellsSample)))
  }

  proj <- addCellColData(ArchRProj = proj, data = paste0(annotationToAdd), cells = paste0(cellsNamesToAdd), name = colnames(metadata)[j], force = TRUE)
}
saveArchRProject(ArchRProj = proj, outputDirectory = "ATACseq", load = FALSE, overwrite = FALSE)

p2 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "nFrags",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p2,p3, name = "QC-Sample-Statistics-TSS-Sample-Violin.pdf", ArchRProj = proj, addDOC = FALSE, width = 8, height = 6)

# For reproducibility, can use previously defined set of cells after doublet removal
proj <- filterDoublets(proj, filterRatio = 1.2)

p2 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "nFrags",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p2,p3, name = "QC-Sample-Statistics-TSS-Sample-Violin-Doublets-Filtered.pdf", ArchRProj = proj, addDOC = FALSE, width = 18, height = 14)

pal <- paletteDiscrete(unique(getCellColData(proj)$Sample))
for (i in names(pal)){
  pal[i] <- "#D51F26"
}
p2 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  name = "TSSEnrichment",
  colorBy = "cellColData",
  plotAs = "violin",
  alpha = 1,
  addBoxPlot = FALSE, pal = pal)
p2 <- p2+geom_boxplot(outlier.shape = NA, alpha = 1)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))

plotPDF(p2, name = "Temp-QC-Sample-Statistics-TSS-Sample-Violin-Doublets-Filtered.pdf", ArchRProj = proj, addDOC = FALSE, width = 30, height = 20)


####step 5
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000, #nCells(proj),
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  sampleCellsPre = 50000,
  sampleCellsFinal = 50000, force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine", force = TRUE
)

# Identify Clusters from Iterative LSI
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1.7, sampleCells = 50000
)

cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
plotPDF(p, name = "Sample-Cluster-Confusion-Matrix", width = 6, height = 6,  ArchRProj = proj, addDOC = FALSE)


# Plot the UMAP Embedding with Metadata Overlayed such as Experimental Sample and Clusters.
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DiseaseState", embedding = "UMAP")
p8 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p9 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
p10 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Location", embedding = "UMAP")
p11 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Donor", embedding = "UMAP")
plotPDF(p1,p2,p3,p8,p9,p10,p11, name = "UMAP-Samples-Clusters", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)

p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DiseaseState", pal = c("#89288F", "#D7CEC7", "#D51F26", "#D7CEC7"), embedding = "UMAP")
plotPDF(p5, name = "UMAP-Disease_State_Only2", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)

p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DiseaseState", pal = c("#D51F26", "#D7CEC7", "#89288F", "#D7CEC7"), embedding = "UMAP")
plotPDF(p5, name = "UMAP-Disease_State_Only3", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj, outputDirectory = "all_samples_HuBMAP_HTAN", load = FALSE, overwrite = FALSE)




