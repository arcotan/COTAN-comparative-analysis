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

datasetName = 'E175'
datasetPath = paste("./data/", datasetName, sep='')

outDir = paste(datasetPath, '/COTAN/', sep='')
outDirDefault = paste(outDir, 'default/', sep='')
outDirCelltypist = paste(outDir, 'celltypist/', sep='')

if (!dir.exists(outDirDefault)) {
  dir.create(outDirDefault, recursive = TRUE)
}
if (!dir.exists(outDirCelltypist)) {
  dir.create(outDirCelltypist, recursive = TRUE)
}

numClusterThreshold = 1/10
numClusterCelltypist = fromJSON(file=paste(datasetPath, '/nclusters.json', sep=''))$nclusters
minNumClusterCelltypist = round(numClusterCelltypist - numClusterThreshold*numClusterCelltypist)
maxNumClusterCelltypist = round(numClusterCelltypist + numClusterThreshold*numClusterCelltypist)

setLoggingLevel(2)
setLoggingFile(file.path(outDir, paste(datasetName, ".log", sep='')))

# Data loading

matrix <- read.csv(paste(datasetPath, 'CorticalCells_GSM2861514_E175_cleaned.csv', sep='/'), header = TRUE, row.names = 1)
sampleCond <- datasetName
e175 <- COTAN(raw = matrix)
e175 <- initializeMetaDataset(
  e175,
  GEO = datasetName,
  sequencingMethod = "",
  sampleCond = sampleCond
)

# Calculate genes’ COEX

e175 <- proceedToCoex(
  e175,
  calcCoex = TRUE,
  cores = numCores,
  saveObj = TRUE,
  outDir = outDir
)

gdiData <- calculateGDI(e175)
genesToLabel <- head(
  rownames(gdiData[order(gdiData[["GDI"]], decreasing = TRUE), ]),
  n = 50L
)
sort(genesToLabel)

gdiData[genesToLabel[[50L]], "GDI"]

gdiPlot <- GDIPlot(
  e175,
  GDIIn = gdiData,
  GDIThreshold = defaultGDIThreshold,
  genes = list("Top 10 GDI genes" = genesToLabel[1L:10L])
)
plot(gdiPlot)

# Clustering

c(splitClusters, splitCoexDF) %<-%
  cellsUniformClustering(
    e175,
    GDIThreshold = defaultGDIThreshold,
    cores = numCores,
    saveObj = TRUE,
    outDir = outDir
  )

e175 <- addClusterization(
  e175,
  clName = "split",
  clusters = splitClusters,
  coexDF = splitCoexDF,
  override = TRUE
)

splitClusters <- getClusterizationData(e175, clName = "split")[[1]]
table(splitClusters)

# Save the clustering result

clustersDf = as.data.frame(splitClusters)
colnames(clustersDf)[colnames(clustersDf) == deparse(substitute(splitClusters))] <- "cluster"
clustersDf$cluster = match(clustersDf$cluster, levels(splitClusters))
clustersDf$cell = rownames(clustersDf)
rownames(clustersDf) = c(1:nrow(clustersDf))
write_clustering(outDirDefault, clustersDf, "cell", "cluster")

# Merge uniform clusters

c(mergedClusters, mergedCoexDF) %<-%
  mergeUniformCellsClusters(
    e175,
    clusters = splitClusters,
    GDIThreshold = 2,#defaultGDIThreshold,
    cores = numCores,
    saveObj = TRUE,
    outDir = outDir
  )

e175 <- addClusterization(
  e175,
  clName = "merged",
  clusters = mergedClusters,
  coexDF = mergedCoexDF,
  override = TRUE
)

mergedClusters <- getClusterizationData(e175, clName = "merged")[[1]]
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

coexDF <- DEAOnClusters(e175, clusters = mergedClusters)
e175 <- addClusterizationCoex(
  e175,
  clName = "merged",
  coexDF = coexDF
)

saveRDS(e175, file = file.path(outDirDefault, paste0(datasetName, ".cotan.RDS")))

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

numCells <- getNumCells(e175)
pvalsDefault <- pValueFromDEA(coexDF, numCells)
write_cotan_markes(outDirDefault, mergedClusters, pvalsDefault, numTopMarkers, FALSE)

# Merge clusters according to celltypist

# binary search on GDIThreshold to match the number of clusters found by celltypist
maxGDIThreshold = 3
minGDIThreshold = 0
repeat{
  GDIThreshold = (maxGDIThreshold + minGDIThreshold) / 2
  cat(paste("Trying GDIThreshold ", GDIThreshold, '\n', sep=''))
  c(celltypistClusters, celltypistCoexDF) %<-%
  mergeUniformCellsClusters(
    e175,
    clusters = splitClusters,
    GDIThreshold = GDIThreshold,
    cores = numCores,
    saveObj = TRUE,
    outDir = outDirCelltypist
  )
  e175 <- addClusterization(
    e175,
    clName = "merged_celltypist",
    clusters = celltypistClusters,
    coexDF = celltypistCoexDF,
    override = TRUE
  )
  celltypistClusters <- getClusterizationData(e175, clName = "merged_celltypist")[[1]]
  numClusters = nrow(table(celltypistClusters))
  cat(paste("Got ", numClusters, ' clusters', sep=''))
  if (numClusters >= minNumClusterCelltypist & numClusters <= maxNumClusterCelltypist){
    cat(paste("GDI: ", GDIThreshold, sep=''))
    write(toJSON(list(GDI_threshold=GDIThreshold)), file=paste(outDirCelltypist, 'GDIthreshold.json',sep=''))
    break
  }
  else if (numClusters < minNumClusterCelltypist){
    # increase GDIThreshold
    maxGDIThreshold = GDIThreshold
  }
  else{
    # reduce GDIThreshold
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

celltypistCoexDF <- DEAOnClusters(e175, clusters = celltypistClusters)
e175 <- addClusterizationCoex(
  e175,
  clName = "merged_celltypist",
  coexDF = celltypistCoexDF
)

saveRDS(e175, file = file.path(outDirCelltypist, paste0(datasetName, ".cotan.RDS")))

# Save markers

celltypistNumCells <- getNumCells(e175)
celltypistPvals <- pValueFromDEA(celltypistCoexDF, celltypistNumCells)
write_cotan_markes(outDirCelltypist, celltypistClusters, celltypistPvals, numTopMarkers, FALSE)


sessionInfo()
