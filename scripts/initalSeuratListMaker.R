source("scripts/generic.R")
source("scripts/adtProcessing.R")


.processFromCellranger <- function(
  dataDir,
  sampleID,
  runN,
  plotDir,
  minCells = 3,
  minFeatures = 200
) {
  
  seuTmp <- Read10X(data.dir = dataDir, strip.suffix = TRUE)
  seuTmp <- CreateSeuratObject(
    counts = seuTmp,
    min.cells = minCells,
    min.features = minFeatures)
  
  seuTmp <- RenameCells(seuTmp, add.cell.id = sampleID)
  seuTmp <- AddMetaData(seuTmp, runN, col.name = "runN")
  seuTmp <- AddMetaData(seuTmp, sampleID, col.name = "hashedSampledID")
  
  seuTmp[["percent.mt"]] <- PercentageFeatureSet(seuTmp, pattern = "^MT-")
  ggsave(
    VlnPlot(seuTmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),
    filename = paste0(plotDir, "/", sampleID, "_seuOrigVlnPlts.png"),
    height = 5,
    width = 8)
  
  return(seuTmp)
}

.processHTOFromKallisto <- function(
  smpl,
  seuTmp,
  adt,
  htoOutDir
) {
  
  hto <- createMatFromKallisto(
    sample = smpl,
    dir = htoOutDir,
    dataMode = "hto",
    sampleSuffix = "_hto_output_count"
  )
  
  colnames(hto) <- paste0(smpl, "_", colnames(hto))
  
  seuTmp <- subset(seuTmp, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 12.5)
  allowableCbcs <- Reduce(intersect, list(colnames(seuTmp), colnames(hto)))
  seuTmp <- subset(seuTmp, cells = allowableCbcs)
  
  seuTmp <- NormalizeData(seuTmp)
  
  seuTmp[["hto"]] <- CreateAssayObject(counts = hto[, allowableCbcs])
  
  seuTmp <- NormalizeData(seuTmp, assay = "hto", normalization.method = "CLR")
  seuTmp <- MULTIseqDemux(seuTmp, assay = "hto", autoThresh = TRUE)
  
  return(seuTmp)
}


initialSeuratListMaker <- function(
  version,
  sampleIds,
  sampleRuns,
  rnaOutDir,
  adtOutDir,
  htoOutDir,
  plotDir = "outs/qc"
) {
  
  seuList <- lapply(seq_along(sampleIds), function(i) {
    x <- sampleIds[i]
    runN <- sampleRuns[i]
    
    print(paste0("working on: ", x))

    seuFileRdsFn <- paste0("rds/seuOrig", version, "/", x, "_seuOrig.rds")
    seuDemuxRdsFn <- paste0("rds/seuDemux", version, "/", x, "_seuDemux.rds")
    
    if (!file.exists(seuDemuxRdsFn)) {
      seu <- rdsCheckAndRun(
        fn = seuFileRdsFn,
        f = .processFromCellranger,
        dataDir = paste0(rnaOutDir, "/run_count_", x, "__RNA/outs/filtered_feature_bc_matrix"),
        sampleID = x,
        runN = runN,
        plotDir = plotDir)
      
      # adt creation
      seu <- rdsCheckAndRun(
        fn = seuDemuxRdsFn,
        f = .processHTOFromKallisto,
        smpl = x,
        seuTmp = seu,
        htoOutDir = htoOutDir)
      
      return(seuDemuxRdsFn)
    } else {
      # seu <- readRDS(seuDemuxRdsFn)
      return(seuDemuxRdsFn)
    }
  })
  
  return(seuList)  
}

  