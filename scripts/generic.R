BASEPTFONTSIZE <- 8
BASEFONTSIZE <- BASEPTFONTSIZE / ggplot2:::.pt

DISEASESTATUSCOLORS <- c(
  "ND" = "#1b9e77",
  "AAb+" = "#7570b3",
  "T1D" = "#d95f02"
)


rdsCheckAndRun <- function(fn, f, ...) {
  if (!file.exists(fn)) {
    print(paste0("File (", fn, ") does not exist. Generating now."))
    
    res <- do.call(f, list(...))
    
    print(paste0("Saving file (", fn, ")..."))
    saveRDS(res, fn)
    
  } else {
    print(paste0("File (", fn, ") exists. Loading now."))
    res <- readRDS(fn)
  }
  
  return(res)
}

savePlot <- function(
  plot,
  fn,
  parentOutDir = "outs",
  devices = c("rds", "png"),
  gheight,
  gwidth,
  rdsPlot = NULL,
  scale = 1,
  customSavePlot = NULL
) {
  
  if (!is.vector(devices)) {
    devices <- c(devices)
  }
  
  for (d in devices) {
    gfn <- glue("{parentOutDir}/{d}/{fn}.{d}")
    
    if (d == "rds" & !is.null(rdsPlot)) {
      saveRDS(rdsPlot, gfn)
    } else if (d == "rds") {
      saveRDS(plot, gfn)
    } else if (!is.null(customSavePlot)) {
      ggsave(gfn, plot = customSavePlot, dpi = "retina", device = d, width = gwidth, height = gheight, scale = scale)
    } else {
      ggsave(gfn, plot = plot, dpi = "retina", device = d, width = gwidth, height = gheight, scale = scale)  
    }
  }
}

### GENERIC THEMES

textSizeOnlyTheme <- theme(
  plot.title = element_text(size = BASEPTFONTSIZE),
  plot.subtitle = element_text(size = BASEPTFONTSIZE),
  text = element_text(size = BASEPTFONTSIZE),
  legend.text = element_text(size = BASEPTFONTSIZE),
  axis.text = element_text(size = BASEPTFONTSIZE),
  axis.title = element_text(size = BASEPTFONTSIZE)
)
