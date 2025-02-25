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
numClusterCelltypist = 26
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

# Merge uniform clusters

c(mergedClusters, mergedCoexDF) %<-%
  mergeUniformCellsClusters(
    e175,
    clusters = splitClusters,
    GDIThreshold = defaultGDIThreshold,
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

# # Differential expression

# coexDF <- DEAOnClusters(e175, clusters = mergedClusters)
# e175 <- addClusterizationCoex(
#   e175,
#   clName = "merged",
#   coexDF = coexDF
# )

# saveRDS(e175, file = file.path(outDirDefault, paste0(datasetName, ".cotan.RDS")))

# # Save markers

# write_cotan_markes = function(outDir, clusters, coex, numTopMarkers, decreasing) {
#     markers = data.frame()

#     for (cluster_id in (levels(clusters))) {
#         cluster_markers = data.frame()
#         pv = coex[, cluster_id]
#         names(pv) = rownames(coex)
#         sorted_pv = sort(pv, decreasing = decreasing)
#         cluster_markers = data.frame(gene = names(sorted_pv)[1:numTopMarkers],
#                                     cluster = match(cluster_id, levels(clusters)),
#                                     rank = 1:numTopMarkers)
#         markers = rbind(markers, cluster_markers)
#     }
#     colnames(markers) = c("gene","cluster","rank")

#     write.csv(markers, paste(outDir, "markers.csv", sep=''), row.names = FALSE)
# }

# numCells <- getNumCells(e175)
# pvalsDefault <- pValueFromDEA(coexDF, numCells)
# write_cotan_markes(outDirDefault, mergedClusters, pvalsDefault, numTopMarkers, FALSE)

# Merge clusters according to celltypist

# get ids of clusters bigger than minClusterSize cells
# mapping = read.csv(paste(datasetPath, '/celltypist/celltypist_mapping.csv', sep=''))
# counts = read.csv(paste(datasetPath, '/celltypist/celltypist_annotation_counts.csv', sep=''))
# mappingCounts = merge(mapping, counts, by.x = "go", by.y = "cluster.ids")
# mappingCounts = subset(mappingCounts, count > minClusterSize)
# clusterIdsToKeep = mappingCounts$id

# # get barcodes of cells in clusters bigger than minClusterSize cells
# celltypistLabels = read.csv(paste(datasetPath, '/celltypist/celltypist_labels.csv', sep=''))
# celltypistLabels = subset(celltypistLabels, cluster.ids %in% clusterIdsToKeep)
# barcodesToKeep = celltypistLabels$cell
# barcodesToKeep = substring(barcodesToKeep, 1, nchar(barcodesToKeep) - 2)

# # keep only cells in clusters bigger than minClusterSize cells
# barcodesToDrop = getCells(e175)[!getCells(e175) %in% barcodesToKeep]
# if (length(barcodesToDrop) != 0) {
#   e175Celltypist <- dropGenesCells(e175, cells = barcodesToDrop)
# } else {
#   e175Celltypist <- e175
# }
# splitClustersCT = splitClusters[!names(splitClusters) %in% barcodesToDrop]

# e175Celltypist <- proceedToCoex(
#   e175Celltypist,
#   calcCoex = FALSE,
#   cores = numCores,
#   saveObj = TRUE,
#   outDir = outDir
# )



# binary search on GDIThreshold to match the number of clusters found by celltypist
cat(paste("\nbinary search on GDIThreshold to match the number of clusters found by celltypist\n", sep=''))
cat(paste("minNumClusterCelltypist:", minNumClusterCelltypist, ' maxNumClusterCelltypist:', maxNumClusterCelltypist, '\n', sep=''))

maxGDIThreshold = 3
minGDIThreshold = 1.35

repeat{
  GDIThreshold = (maxGDIThreshold + minGDIThreshold) / 2
  cat(paste("Trying GDIThreshold ", "\n", GDIThreshold, sep=''))
  c(celltypistClusters, celltypistCoexDF) %<-%
  mergeUniformCellsClusters(
    e175,
    clusters = splitClusters,
    GDIThreshold = GDIThreshold,
    cores = numCores,
    saveObj = TRUE,
    outDir = outDirCelltypist
  )
  e175Celltypist <- addClusterization(
    e175Celltypist,
    clName = "merged_celltypist",
    clusters = celltypistClusters,
    coexDF = celltypistCoexDF,
    override = TRUE
  )
  celltypistClusters <- getClusterizationData(e175Celltypist, clName = "merged_celltypist")[[1]]
  numClusters = nrow(table(celltypistClusters))
  cat(paste("Got ", numClusters, ' clusters\n', sep=''))
  cat(paste("GDI: ", GDIThreshold,'\n', sep=''))
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

# # Differential expression

# celltypistCoexDF <- DEAOnClusters(e175Celltypist, clusters = celltypistClusters)
# e175Celltypist <- addClusterizationCoex(
#   e175Celltypist,
#   clName = "merged_celltypist",
#   coexDF = celltypistCoexDF
# )

# saveRDS(e175, file = file.path(outDirCelltypist, paste0(datasetName, ".cotan.RDS")))

# # Save markers

# celltypistNumCells <- getNumCells(e175Celltypist)
# celltypistPvals <- pValueFromDEA(celltypistCoexDF, celltypistNumCells)
# write_cotan_markes(outDirCelltypist, celltypistClusters, celltypistPvals, numTopMarkers, FALSE)


sessionInfo()