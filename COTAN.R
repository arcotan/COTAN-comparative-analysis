pid <- Sys.getpid()
cat("Process ID:", pid, "\n")

library(COTAN)
library(Seurat)
library(tibble)
library(ggplot2)
library(zeallot)
library(rjson)
source("./utils.R")
options(parallelly.fork.enable = TRUE)

numCores = 4
defaultGDIThreshold = 1.4
numTopMarkers = 500
minClusterSize = 0

datasetName = 'PBMC4' # modify this
datasetPath = paste("./data/", datasetName, sep='')
inDir = paste(datasetPath, '/filtered/', sep='')
outDir = paste(datasetPath, '/COTAN/', sep='')
outDirDefault = paste(outDir, 'default/', sep='')
outDirCelltypist = paste(outDir, 'celltypist/', sep='')
outDirAntibody = paste(outDir, 'antibody/', sep='')

if (!dir.exists(outDirDefault)) {
  dir.create(outDirDefault, recursive = TRUE)
}
if (!dir.exists(outDirCelltypist)) {
  dir.create(outDirCelltypist, recursive = TRUE)
}
if (!dir.exists(outDirAntibody)) {
  dir.create(outDirAntibody, recursive = TRUE)
}

numClusterThreshold = 1/10
numClusterCelltypist = fromJSON(file=paste(datasetPath, '/celltypist/nclusters.json', sep=''))$nclusters
minNumClusterCelltypist = round(numClusterCelltypist - numClusterThreshold*numClusterCelltypist)
maxNumClusterCelltypist = round(numClusterCelltypist + numClusterThreshold*numClusterCelltypist)
numClusterAntibody = fromJSON(file=paste(datasetPath, '/antibody_annotation/nclusters_postproc.json', sep=''))$nclusters
minNumClusterAntibody = round(numClusterAntibody - numClusterThreshold*numClusterAntibody)
maxNumClusterAntibody = round(numClusterAntibody + numClusterThreshold*numClusterAntibody)

setLoggingLevel(2)
setLoggingFile(file.path(outDir, paste(datasetName, ".log", sep='')))


# Data loading


PBMC <- readRDS(file = file.path(inDir, paste0(datasetName, ".cotan.RDS")))


# Calculate genes’ COEX


PBMC <- proceedToCoex(
  PBMC,
  calcCoex = TRUE,
  cores = numCores,
  saveObj = TRUE,
  outDir = outDir
)



gdiData <- calculateGDI(PBMC)
genesToLabel <- head(
  rownames(gdiData[order(gdiData[["GDI"]], decreasing = TRUE), ]),
  n = 50L
)
sort(genesToLabel)



gdiData[genesToLabel[[50L]], "GDI"]



gdiPlot <- GDIPlot(
  PBMC,
  GDIIn = gdiData,
  GDIThreshold = defaultGDIThreshold,
  genes = list("Top 10 GDI genes" = genesToLabel[1L:10L])
)
plot(gdiPlot)


# Clustering


c(splitClusters, splitCoexDF) %<-%
  cellsUniformClustering(
    PBMC,
    GDIThreshold = defaultGDIThreshold,
    cores = numCores,
    saveObj = TRUE,
    outDir = outDir
  )

PBMC <- addClusterization(
  PBMC,
  clName = "split",
  clusters = splitClusters,
  coexDF = splitCoexDF,
  override = TRUE
)



splitClusters <- getClusterizationData(PBMC, clName = "split")[[1]]
table(splitClusters)


# Merge uniform clusters

c(mergedClusters, mergedCoexDF) %<-%
  mergeUniformCellsClusters(
    PBMC,
    clusters = splitClusters,
    GDIThreshold = defaultGDIThreshold,
    cores = numCores,
    saveObj = TRUE,
    outDir = outDir
  )

PBMC <- addClusterization(
  PBMC,
  clName = "merged",
  clusters = mergedClusters,
  coexDF = mergedCoexDF,
  override = TRUE
)

mergedClusters <- getClusterizationData(PBMC, clName = "merged")[[1]]
clusters_sizes = table(mergedClusters)
cat(paste('Number of clusters found: ', nrow(clusters_sizes), sep=''))

# Save the clustering result

clustersDf = as.data.frame(mergedClusters)
colnames(clustersDf)[colnames(clustersDf) == deparse(substitute(mergedClusters))] <- "cluster"
clustersDf$cluster = match(clustersDf$cluster, levels(mergedClusters))
clustersDf$cell = rownames(clustersDf)
rownames(clustersDf) = c(1:nrow(clustersDf))
write_clustering(outDirDefault, clustersDf, "cell", "cluster")

# Differential expression

coexDF <- DEAOnClusters(PBMC, clusters = mergedClusters)
PBMC <- addClusterizationCoex(
  PBMC,
  clName = "merged",
  coexDF = coexDF
)

saveRDS(PBMC, file = file.path(outDirDefault, paste0(datasetName, ".cotan.RDS")))

# Save markers

write_cotan_markes = function(outDir, clusters, coex, numTopMarkers, decreasing) {
    markers = data.frame()

    for (cluster_id in (levels(clusters))) {
        cluster_markers = data.frame()
        pv = coex[, cluster_id]
        names(pv) = rownames(coex)
        sorted_pv = sort(pv, decreasing = decreasing)
        cluster_markers = data.frame(gene = names(sorted_pv)[1:numTopMarkers],
                                    cluster = match(cluster_id, levels(clusters)),
                                    rank = 1:numTopMarkers)
        markers = rbind(markers, cluster_markers)
    }
    colnames(markers) = c("gene","cluster","rank")

    write.csv(markers, paste(outDir, "markers.csv", sep=''), row.names = FALSE)
}

numCells <- getNumCells(PBMC)
pvalsDefault <- pValueFromDEA(coexDF, numCells)
write_cotan_markes(outDirDefault, mergedClusters, pvalsDefault, numTopMarkers, FALSE)

# Merge clusters according to celltypist

# get ids of clusters bigger than minClusterSize cells
mapping = read.csv(paste(datasetPath, '/celltypist/celltypist_mapping.csv', sep=''))
counts = read.csv(paste(datasetPath, '/celltypist/celltypist_annotation_counts.csv', sep=''))
mappingCounts = merge(mapping, counts, by.x = "go", by.y = "cluster.ids")
mappingCounts = subset(mappingCounts, count > minClusterSize)
clusterIdsToKeep = mappingCounts$id

# get barcodes of cells in clusters bigger than minClusterSize cells
celltypistLabels = read.csv(paste(datasetPath, '/celltypist/celltypist_labels.csv', sep=''))
celltypistLabels = subset(celltypistLabels, cluster.ids %in% clusterIdsToKeep)
barcodesToKeep = celltypistLabels$cell
barcodesToKeep = substring(barcodesToKeep, 1, nchar(barcodesToKeep) - 2)

# keep only cells in clusters bigger than minClusterSize cells
barcodesToDrop = getCells(PBMC)[!getCells(PBMC) %in% barcodesToKeep]
if (length(barcodesToDrop) != 0) {
  PBMCCelltypist <- dropGenesCells(PBMC, cells = barcodesToDrop)
} else {
  PBMCCelltypist <- PBMC
}

# binary search on GDIThreshold to match the number of clusters found by celltypist
cat(paste("binary search on GDIThreshold to match the number of clusters found by celltypist\n", sep=''))
cat(paste("minNumClusterCelltypist:", minNumClusterCelltypist, ' maxNumClusterCelltypist:', maxNumClusterCelltypist, '\n', sep=''))

maxGDIThreshold = 3
minGDIThreshold = 1.35

repeat{
  GDIThreshold = (maxGDIThreshold + minGDIThreshold) / 2
  cat(paste("Trying GDIThreshold ", GDIThreshold, sep=''))
  c(celltypistClusters, celltypistCoexDF) %<-%
  mergeUniformCellsClusters(
    PBMCCelltypist,
    clusters = splitClusters,
    GDIThreshold = GDIThreshold,
    cores = numCores,
    saveObj = TRUE,
    outDir = outDirCelltypist
  )
  PBMCCelltypist <- addClusterization(
    PBMCCelltypist,
    clName = "merged_celltypist",
    clusters = celltypistClusters,
    coexDF = celltypistCoexDF,
    override = TRUE
  )
  celltypistClusters <- getClusterizationData(PBMCCelltypist, clName = "merged_celltypist")[[1]]
  numClusters = nrow(table(celltypistClusters))
  cat(paste("Got ", numClusters, ' clusters\n', sep=''))
  cat(paste("GDI: ", GDIThreshold,'\n' sep=''))
  if (numClusters >= minNumClusterCelltypist & numClusters <= maxNumClusterCelltypist){
    cat(paste("exit the cycle\n", sep=''))
    write(toJSON(list(GDI_threshold=GDIThreshold)), file=paste(outDirCelltypist, 'GDIthreshold.json',sep=''))
    break
  }
  else if (numClusters < minNumClusterCelltypist){
    cat(paste("reducing the GDI \n", sep=''))
    maxGDIThreshold = GDIThreshold
  }
  else{
    cat(paste("increasing the GDI \n", sep=''))
    minGDIThreshold = GDIThreshold
  }
}

# Save the clustering result

celltypistClustersDf = as.data.frame(celltypistClusters)
colnames(celltypistClustersDf)[colnames(celltypistClustersDf) == deparse(substitute(celltypistClusters))] <- "cluster"
celltypistClustersDf$cluster = match(celltypistClustersDf$cluster, levels(celltypistClusters))
celltypistClustersDf$cell = rownames(celltypistClustersDf)
rownames(celltypistClustersDf) = c(1:nrow(celltypistClustersDf))
write_clustering(outDirCelltypist, celltypistClustersDf, "cell", "cluster")

# Differential expression

celltypistCoexDF <- DEAOnClusters(PBMCCelltypist, clusters = celltypistClusters)
PBMCCelltypist <- addClusterizationCoex(
  PBMCCelltypist,
  clName = "merged_celltypist",
  coexDF = celltypistCoexDF
)

saveRDS(PBMC, file = file.path(outDirCelltypist, paste0(datasetName, ".cotan.RDS")))

# Save markers

celltypistNumCells <- getNumCells(PBMCCelltypist)
celltypistPvals <- pValueFromDEA(celltypistCoexDF, celltypistNumCells)
write_cotan_markes(outDirCelltypist, celltypistClusters, celltypistPvals, numTopMarkers, FALSE)

# Merge clusters according to protein surface

# get barcodes of cells labelled using antibody data
antibodyLabels = read.csv(paste(datasetPath, '/antibody_annotation/antibody_labels_postproc.csv', sep=''))
barcodesToKeep = antibodyLabels$cell

# barcodesToKeep = substring(barcodesToKeep, 1, nchar(barcodesToKeep) - 2) # TODO: ?

# keep only cells labelled with antibody data
barcodesToDrop = getCells(PBMC)[!getCells(PBMC) %in% barcodesToKeep]
if (length(barcodesToDrop) != 0) {
  PBMCAntibody <- dropGenesCells(PBMC, cells = barcodesToDrop)
} else {
  PBMCAntibody <- PBMC
}


# binary search on GDIThreshold to match the number of clusters found by antibody
cat(paste("binary search on GDIThreshold to match the number of clusters found by celltypist\n", sep=''))
cat(paste("minNumClusterAntibody:", minNumClusterAntibody, ' maxNumClusterAntibody:', maxNumClusterAntibody, '\n', sep=''))

maxGDIThreshold = 3
minGDIThreshold = 1.35
repeat{
  GDIThreshold = (maxGDIThreshold + minGDIThreshold) / 2
  cat(paste("Trying GDIThreshold ", GDIThreshold, sep=''))
  c(abtibodyClusters, antibodyCoexDF) %<-%
  mergeUniformCellsClusters(
    PBMCAntibody,
    clusters = splitClusters,
    GDIThreshold = GDIThreshold,
    cores = numCores,
    saveObj = TRUE,
    outDir = outDirAntibody
  )
  PBMCAntibody <- addClusterization(
    PBMCAntibody,
    clName = "merged_antibody",
    clusters = antibodyClusters,
    coexDF = antibodyCoexDF,
    override = TRUE
  )
  antibodyClusters <- getClusterizationData(PBMCAntibody, clName = "merged_antibody")[[1]]
  numClusters = nrow(table(antibodyClusters))
  cat(paste("Got ", numClusters, ' clusters\n', sep=''))
  cat(paste("Got ", numClusters, ' clusters\n', sep=''))
  if (numClusters >= minNumClusterAntibody & numClusters <= maxNumClusterAntibody){
    write(toJSON(list(GDI_threshold=GDIThreshold)), file=paste(outDirAntibody, 'GDIthreshold.json',sep=''))
    cat(paste("exit the cycle\n", sep=''))
    break
  }
  else if (numClusters < minNumClusterAntibody){
    cat(paste("reducing the GDI \n", sep=''))
    maxGDIThreshold = GDIThreshold
  }
  else{
    cat(paste("increasing the GDI \n", sep=''))
    minGDIThreshold = GDIThreshold
  }
}

# Save the clustering result

antibodyClustersDf = as.data.frame(antibodyClusters)
colnames(antibodyClustersDf)[colnames(antibodyClustersDf) == deparse(substitute(antibodyClusters))] <- "cluster"
antibodyClustersDf$cluster = match(antibodyClustersDf$cluster, levels(antibodyClusters))
antibodyClustersDf$cell = rownames(antibodyClustersDf)
rownames(antibodyClustersDf) = c(1:nrow(antibodyClustersDf))
write_clustering(outDirAntibody, antibodyClustersDf, "cell", "cluster")

# Differential expression

antibodyCoexDF <- DEAOnClusters(PBMCAntibody, clusters = antibodyClusters)
PBMCAntibody <- addClusterizationCoex(
  PBMCAntibody,
  clName = "merged_antibody",
  coexDF = antibodyCoexDF
)

saveRDS(PBMC, file = file.path(outDirAntibody, paste0(datasetName, ".cotan.RDS")))

# Save markers

numCells <- getNumCells(PBMCAntibody)
antibodyPvals <- pValueFromDEA(antibodyCoexDF, numCells)
write_cotan_markes(OUT_PROTEIN_DIR, antibodyClusters, antibodyPvals, numTopMarkers, FALSE)

sessionInfo()