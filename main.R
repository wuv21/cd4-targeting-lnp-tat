# uncomment if running in RProject
# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

library("dplyr")
library("tidyr")
library("stringr")
library("ggplot2")
library("Seurat")
library("glue")
library("patchwork")
library("ComplexHeatmap")

# general scripts
source("scripts/generic.R")
source("scripts/adtProcessing.R")
source("scripts/initalSeuratListMaker.R")

# start of main code

set.seed(42)

rnaOut <- "data/rna"
htoOut <- "data/hto"

unhashedSampleMeta <- read.csv("unhashed_sample_meta.csv")


# load merged dataset if it exists
seuMergedFn <- "rds/seuMerged.rds"
harmonizedFn <- "rds/seuMergedAndHarmonized.rds"
umapRdsFn <- "rds/seuMergedAndHarmonized_withUMAP.rds"

if (file.exists(seuMergedFn) & file.exists(harmonizedFn)) {
  message("harmonized and merged file already exists, so skipping initial merging")
  
} else if (file.exists(seuMergedFn)) {
  message("loading previously merged rds file")
  seuMerged <- readRDS(seuMergedFn)
  
} else {
  seuListFns <- rdsCheckAndRun(
    fn = "rds/allSeuList_v1.rds",
    f = initialSeuratListMaker,
    version = "Merge",
    sampleIds = unhashedSampleMeta$tube_ID,
    sampleRuns = unhashedSampleMeta$chip,
    rnaOutDir = rnaOut,
    htoOutDir = htoOut,
    plotDir = "outs/qc")
  
  seuList <- lapply(seuListFns, function(x) {
    return(readRDS(x))
  })
  
  names(seuList) <- unhashedSampleMeta$tube_ID
  
  sampleMeta <- read.csv("sample_meta.csv")
  sampleMetaDict <- sampleMeta %>%
    group_split(pool_type)
  
  names(sampleMetaDict) <- sapply(sampleMetaDict, function(x) {return(as.character(x[1, "pool_type"]))})
  poolNumDict <- c("Nil pool", "vehicle pool", "luc pool", "tat pool")
  
  for (i in seq(1, length(seuList))) {
    xName <- names(seuList)[i]
    
    message(paste0("working on ", xName))
    
    seuList[[i]] <- subset(seuList[[i]], subset = MULTI_ID != "Doublet" & MULTI_ID != "Negative")
    
    # add metadata
    poolNum <- as.numeric(substr(as.character(seuList[[i]]$hashedSampledID[1]), 1, 1))
    poolType <- poolNumDict[poolNum]
    
    seuDf <- data.frame(
      HTO_antibody = as.character(seuList[[i]]$MULTI_ID)
    ) %>%
      left_join(sampleMetaDict[[poolType]]) %>%
      select(-HTO_antibody)
    
    rownames(seuDf) <- colnames(seuList[[i]])
    seuList[[i]] <- AddMetaData(seuList[[i]], seuDf)
  }
  
  message('merging datasets')
  # filter out specific donors due to poor viability during thaw and initial qc
  seuList <- lapply(seuList, function(x) {
    message("finding variable features.")
    seuList[[i]] <- FindVariableFeatures(seuList[[i]])
    
    return(x)
  })
  
  seuMerged <- merge(seuList[[1]], seuList[2:length(seuList)], merge.data = TRUE)
  
  seuMerged <- JoinLayers(seuMerged)
  
  VariableFeatures(seuMerged, assay = "RNA") <- SelectIntegrationFeatures(
    object.list = seuList,
    assay = rep("RNA", length(seuList)))
  
  # remove TAT86R5M4 and LUC2 from being factored into PCA calculations
  varFeatures <- VariableFeatures(seuMerged, assay = "RNA")
  VariableFeatures(seuMerged, assay = "RNA") <- varFeatures[!(varFeatures %in% c("TAT86R5M4", "LUC2"))]
  
  seuMerged <- ScaleData(seuMerged, assay = "RNA") %>%
    RunPCA(assay = "RNA")
  
  saveRDS(seuMerged, seuMergedFn)
  
  rm(seuList)
  gc()
}

################################################################################
# harmony batch effect correction against run number and donor ID.
################################################################################
if (file.exists(harmonizedFn) & file.exists(umapRdsFn)) {
  message("next rds (clustered + umap) already found so skipping import of harmonized rds")
  
} else if (file.exists(harmonizedFn)) {
  message("loading previously merged and harmonized rds")
  seuMerged <- readRDS(harmonizedFn)
} else {
  # ensure that all metadata is character
  # for some reason, harmony doesn't like integer metadata
  harmonyGroupingVars <- c("runN", "donor")
  
  for (v in harmonyGroupingVars) {
    tmpHolder <- as.character(seuMerged[[v]][, 1])
    seuMerged <- AddMetaData(seuMerged, tmpHolder, col.name = v)
    
    stopifnot(sum(tmpHolder == seuMerged[[v]][, 1]) == length(tmpHolder))
  }
  
  message('harmonizing RNA')
  seuMerged <- harmony::RunHarmony(
    object = seuMerged,
    reduction.save = "harmonyRNA",
    reduction = "pca",
    assay.use = "RNA",
    project.dim = FALSE,
    dims.use = 1:50,
    verbose = TRUE,
    group.by.vars = harmonyGroupingVars)
  
  saveRDS(seuMerged, harmonizedFn)
}

################################################################################
# cluster finding by RNA modality only
################################################################################
if (!file.exists(umapRdsFn)) {
  seuMerged <- FindNeighbors(seuMerged, assay = "RNA", reduction = "harmonyRNA", dims = 1:30)
  seuMerged <- FindClusters(
    seuMerged,
    algorithm = 2,
    # method = "igraph",
    verbose = TRUE)
  
  saveRDS(seuMerged, "rds/seuMerged_withClusters")
  
  seuMerged <- RunUMAP(seuMerged,
    dims = 1:30,
    reduction.name = "rna.umap",
    reduction.key = "rnaUMAP_",
    n.neighbors = 50,
    min.dist = 0.2,
    reduction = "harmonyRNA")
  
  saveRDS(seuMerged, umapRdsFn)
  
} else {
  seuMerged <- readRDS(umapRdsFn)
}

p1 <- DimPlot(seuMerged, group.by = "pool_type")
p2 <- DimPlot(seuMerged, group.by = "runN")
p3 <- DimPlot(seuMerged, group.by = "donor")
p4 <- DimPlot(seuMerged, group.by = "RNA_snn_res.0.8", label = TRUE)
# 
# # standard phenotypic
FeaturePlot(seuMerged,
  features = c("CD3E", "CD8A", "IL27RA", "GZMB")) &
  theme(legend.position = "blank",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank())

RidgePlot(seuMerged,
  features = c("CD3E", "CD8A", "IL27RA", "FCGR3A"),
  group.by = "RNA_snn_res.0.8", ncol = 2)

# transfer annotations
satijaReference <- readRDS("ref/pbmc_multimodal_2023.rds")
anchor <- FindTransferAnchors(
  reference = satijaReference,
  query = seuMerged,
  reference.reduction = "spca",
  normalization.method = "SCT",
  dims = 1:50
)

seuMerged <- MapQuery(
  anchorset = anchor,
  query = seuMerged,
  reference = satijaReference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2"
  ),
  reduction.model = "wnn.umap"
)


save.image("rds/working_state_v1.rdata")

#### start here
load("rds/working_state_v1.rdata")

# rename pools
seuMerged$pool_type <- case_when(
  seuMerged$pool_type == "luc pool" ~ "Ib-Luc",
  seuMerged$pool_type == "tat pool" ~ "Ib-Tat",
  seuMerged$pool_type == "vehicle pool" ~ "Ib-empty",
  seuMerged$pool_type == "Nil pool" ~ "N/A"
)

# fig 1
BASELINEWIDTH <- 12/14

umapTheme <- theme(
  legend.position = "bottom",
  legend.justification = "center",
  legend.spacing.x = unit(3, 'pt'),
  legend.text = element_text(size = BASEPTFONTSIZE - 2),
  legend.title = element_blank(),
  # axis.title = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(linewidth = BASELINEWIDTH, color = "#000000"),
  plot.margin = margin(t = 0, b = 0, l = 3, r = 3),
  plot.title = element_text(size = BASEPTFONTSIZE - 2, color = "#000000", face = "plain", margin = margin(b = -2)),
  legend.key.size = unit(BASEPTFONTSIZE - 4, 'points'),
  legend.key.spacing.y = unit(BASEPTFONTSIZE - 6, "points"),
  legend.background = element_rect(fill = "transparent", colour = NA))

fig1A <- DimPlot(seuMerged, group.by = "pool_type") + umapTheme +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Treatment") +
  guides(color = guide_legend(nrow = 2))

fig1B <- DimPlot(seuMerged, group.by = "donor") + umapTheme +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Donor") +
  guides(color = guide_legend(nrow = 2))

fig1C <- DimPlot(seuMerged, group.by = "predicted.celltype.l1") + umapTheme +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Cell type") +
  guides(color = guide_legend(nrow = 4))

# fig1 <- fig1A + fig1B + fig1C + 
#   plot_layout(nrow = 1, axis_titles = "collect") &
#   theme(
#     legend.position = "bottom",
#     legend.margin = margin(t = -10),
#     axis.title = element_blank())


# fig 2
Idents(seuMerged) <- "pool_type"
tatVsLuc <- FindMarkers(
  subset(seuMerged, subset = predicted.celltype.l1 == "CD4 T" & predicted.celltype.l1.score >= 0.9),
  ident.1 = "Ib-Tat",
  ident.2 = "Ib-Luc",
  min.pct = 0.3
)

tmpForJilian <- FindMarkers(
  subset(seuMerged, subset = predicted.celltype.l1 == "CD4 T" & predicted.celltype.l1.score >= 0.9),
  ident.1 = "Ib-Luc",
  ident.2 = "N/A",
  min.pct = 0.3
)

tmpForJilianMono <- FindMarkers(
  subset(seuMerged, subset = predicted.celltype.l1 == "Mono" & predicted.celltype.l1.score >= 0.9),
  ident.1 = "Ib-Luc",
  ident.2 = "N/A",
  min.pct = 0.3
)

tatVsLucDf <- tatVsLuc %>%
  mutate(gene = rownames(.)) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(abs(avg_log2FC) > 0.25) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "Up in Ib-Tat", "Up in Ib-Luc")) %>%
  dplyr::group_by(direction) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  dplyr::mutate(gene = factor(gene, levels = gene))

seuCD4 <- subset(seuMerged, subset = predicted.celltype.l1 == "CD4 T" & predicted.celltype.l1.score >= 0.9)

tatVsLucHeatdf <- AverageExpression(
  subset(seuCD4, subset = pool_type %in% c("Ib-Tat", "Ib-Luc")),
  features = tatVsLucDf$gene,
  assay = "RNA",
  group.by = c("donor", "pool_type"))

tatVsLucHeatMat <- tatVsLucHeatdf$RNA
tatVsLucHeatMat <- t(scale(t(as.matrix(tatVsLucHeatMat))))

fig2 <- Heatmap(
  matrix = tatVsLucHeatMat,
  name = " ",
  row_title_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4),
  show_column_names = FALSE,
  row_dend_width = unit(0.2, "cm"),
  column_dend_height = unit(0.2, "cm"),
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    legend_width = unit(2, "cm"),
    legend_height = unit(0.2, "cm"),
    labels_gp = gpar(fontsize = 6)),
  bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
    `Donor` = str_split_fixed(colnames(tatVsLucHeatMat), "_", n = 2)[, 1],
    `Treatment` = str_split_fixed(colnames(tatVsLucHeatMat), "_", n = 2)[, 2],
    annotation_legend_param = list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6)
    ),
    annotation_name_gp = gpar(fontsize = 6),
    simple_anno_size = unit(0.3, "cm")
  ))



activationGoi <- c("JUNB", "FOS", "LCK", "MYC", "RELA", "NFKB1")

fig3 <- VlnPlot(seuMerged,
  features = activationGoi,
  group.by = "pool_type", ncol = 2, pt.size = 0.1, alpha = 0.01) +
  plot_layout(axis_titles = "collect", axes = "collect_x") &
  theme(
    plot.margin = margin(t = 2, b = 2, l = 5, r = 5),
    axis.line = element_line(linewidth = BASELINEWIDTH, color = "#000000"),
    axis.line.x = element_line(color = "#000000", linewidth = BASELINEWIDTH),
    plot.title = element_text(color = "#000000", size = BASEPTFONTSIZE - 2, face = "plain", margin = margin(b = -2)),
    axis.text = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    axis.title = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    axis.title.x = element_blank()
)


FindMarkers(
  seuCD4,
  features = activationGoi,
  ident.1 = "Ib-Tat",
  ident.2 = "Ib-Luc",
  min.pct = 0,
  logfc.threshold = 0
)

supFig <- VlnPlot(seuMerged,
  features = c("CCR5", "HLA-DRA", "CD69", "IL2RA", "PDCD1"),
  group.by = "pool_type", ncol = 2, pt.size = 0.1, alpha = 0.01) +
  plot_layout(axis_titles = "collect", axes = "collect") &
  theme(
    plot.margin = margin(t = 2, b = 2, l = 5, r = 5),
    axis.line = element_line(linewidth = BASELINEWIDTH, color = "#000000"),
    plot.title = element_text(color = "#000000", size = BASEPTFONTSIZE - 2, face = "plain"),
    axis.text = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    axis.title = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    axis.title.x = element_blank()
  )

# layout <- c(
#   area(1, 1, 2, 3), # umap
#   area(1, 4, 3, 9), # tracks
#   area(3, 1, 4, 3), # vln plot
#   area(5, 1, 6, 3), # deg volcano 72 hrs
#   area(4, 4, 6, 9) # dot plot
# )
# 
# top_half_p <- list(
#   "A" = figA,
#   "B" = fig2TrackFinal,
#   "C" = fig2C,
#   "D" = fig2d
# )
# saveRDS(top_half_p, "outs/rds/fig3_tophalf.rds")
# 
# final_p <- patchwork::free(figA, side = "trb") + 
#   patchwork::free(fig2TrackFinal, side = "lr") +
#   patchwork::free(fig2C, side = "trb") +
#   patchwork::free(fig2d, side = "trb") +
#   figDotPlotApop +
#   plot_layout(design = layout) + 
#   plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E"))) &
#   theme(plot.tag = element_text(size = 10, color = "#000000", face = "plain"))
# 
# 
# savePlot(plot = final_p, fn = "fig2", devices = c("png", "pdf"), gwidth = 8, gheight = 5)

top_half <- readRDS("/data/tat_multiome/tat_multiome_early_analysis/outs/rds/fig3_tophalf.rds")

layout <- c(
  area(1, 1, 2, 6), # umap
  area(1, 7, 3, 18), # tracks
  area(3, 1, 4, 6), # vln plot
  area(5, 1, 6, 6), # deg volcano 72 hrs
  area(4, 7, 6, 18), # dot plot
  area(7, 1, 9, 9), # umap
  # area(7, 4, 9, 6), # umap
  # area(7, 7, 9, 8), # umap
  area(7, 10, 13, 18), # deg heatmap
  area(10, 1, 13, 9) # goi
)

final_p <- patchwork::free(top_half[["A"]] + theme(axis), side = "trb") + 
  patchwork::free(top_half[["B"]], side = "lr") +
  patchwork::free(top_half[["C"]], side = "trb") +
  patchwork::free(top_half[["D"]] + theme(legend.position = "blank"), side = "trb") +
  patchwork::free(top_half[["E"]], side = "tlrb") +
  patchwork::free(fig1A + fig1B + fig1C & theme(legend.box.margin = margin(t = -10, b = 5)), side = "trb") +
  wrap_elements(
    panel = grid.grabExpr(ComplexHeatmap::draw(
      fig2,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legend = TRUE,
      padding = unit(c(0, 2, 0, 2), "mm")))
  ) +
  patchwork::free(fig3, side = "tblr") +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "", "", "G", "H", "", "", "", "", ""))) &
  theme(plot.tag = element_text(size = 10, color = "#000000", face = "plain"))

savePlot(plot = final_p, fn = "fig3_complete", devices = c("png", "pdf"), gwidth = 8, gheight = 10.2)


# layout <- c(
#   area(1, 1, 2, 4), # umap
#   area(1, 5, 8, 9), # deg heatmap
#   area(3, 1, 8, 4) # goi
# )

# final_p <- free(fig1, side = "b") + 
#   wrap_elements(
#     panel = grid.grabExpr(ComplexHeatmap::draw(
#       fig2,
#       heatmap_legend_side = "bottom",
#       annotation_legend_side = "bottom",
#       merge_legend = TRUE))
#   ) +
#   free(fig3, side = "tblr") +
#   plot_layout(design = layout) + 
#   plot_annotation(tag_levels = list(c("E", "F", "G", "H", "I"))) &
#   theme(plot.tag = element_text(size = 10, color = "#000000", face = "plain"))

