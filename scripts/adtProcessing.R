#' Read in kallisto list
#'
#' @description
#' Saves and returns a RDS containing a list of count matrix objects by reading files from kallisto bustools output.
#'
#' @param sample A character that corresponds to the sample prefix
#' of the kallisto bustools output. For example, "A01_1" will look for
#' "A01_1_count_out", which uses the default "_count_out" suffix.
#'
#' @param dir A character denoting where the ADT output directions are stored.
#' 
#' @param dataMode A character: either "adt" or "hto"
#'
#' @param emptyDropsLower Numeric value that is passed to emptyDrops() lower
#' parameter to setting the lower bound for automatic calling of empty droplets.
#' Default is 100, which is emptyDrops's default value.
#'
#' @param outRDSPath A character denoting where to save the RDS containing the RDS.
#' Default is paste0("rds/", sample, "_post_emptydrops")
#'
#' @param sampleSuffix A character value indicating the suffix of the directory containing the
#' ADT data. Default: _count_out.
#'
#' @return A list of count matrices.


createMatFromKallisto <- function(
  sample,
  dir,
  dataMode,
  emptyDropsLower = 100,
  outRDSPath = paste0("rds/", sample, "_post_emptydrops"),
  sampleSuffix = "_count_out",
  kallistoPrefix = "output")
{
  
  stopifnot(dataMode %in% c("adt", "hto"))
  
  bus_output <- read.table(file = glue("{dir}/{sample}{sampleSuffix}/{kallistoPrefix}.mtx"), header = FALSE, skip = 4)
  feats <- read.table(file = glue("{dir}/{sample}{sampleSuffix}/{kallistoPrefix}.genes.txt"), header = FALSE)
  cbcs <- read.table(file = glue("{dir}/{sample}{sampleSuffix}/{kallistoPrefix}.barcodes.txt"), header = FALSE)
  
  count_matrix <- Matrix::sparseMatrix(
    i = bus_output$V2,
    j = bus_output$V1,
    x = bus_output$V3)
  
  rownames(count_matrix) <- feats$V1
  colnames(count_matrix) <- cbcs$V1
  
  if (dataMode == "adt" && !file.exists(outRDSPath)) {
    empty_calc <- emptyDrops(count_matrix, lower = emptyDropsLower)
    adt_filt <- count_matrix[, which(empty_calc$FDR < 0.01)]
    
    saveRDS(adt_filt, file = outRDSPath)
    return(adt_filt)
    
  } else if (dataMode == "adt") {
    adt <- readRDS(outRDSPath)
    
    return(adt)
  } else {
    return(count_matrix)
  }
}