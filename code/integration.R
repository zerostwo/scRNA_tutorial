#### Information ----
# Title   :   Integrate multiple sample data 
# File    :   integration.R
# Author  :   Songqi Duan
# Contact :   songqi.duan@outlook.com
# License :   Copyright (C) by Songqi Duan
# Created :   2021/06/28 14:05:14
# Updated :   none

#### 加载包 ----
library(Seurat)
library(patchwork)

#### 导入数据 ----
data_dir <- "./data/ifnb/control/"
expr1 <- Read10X(data.dir = data_dir)
dim(expr1)
seurat_obj1 <- CreateSeuratObject(
  counts = expr1,
  project = "control",
  min.cells = 3,
  min.features = 200
)

data_dir <- "./data/ifnb/stim/"
expr2 <- Read10X(data.dir = data_dir)
dim(expr2)
seurat_obj2 <- CreateSeuratObject(
  counts = expr2,
  project = "stim",
  min.cells = 3,
  min.features = 200
)

merged.list <- list(seurat_obj1, seurat_obj2)
rm(seurat_obj1, seurat_obj2, expr1, expr2)
# normalize and identify variable features for each dataset independently
merged.list <- lapply(
  X = merged.list,
  FUN = function(x) {
    x <- NormalizeData(x)
    x <-
      FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  }
)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = merged.list)

immune.anchors <-
  FindIntegrationAnchors(object.list = merged.list, anchor.features = features)

#### Perform integration ----
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

#### Perform an integrated analysis ----
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, verbose = FALSE)
ElbowPlot(immune.combined, ndims = 50)
ndims <- "30"
immune.combined <-
  RunUMAP(immune.combined, reduction = "pca", dims = 1:ndims)
immune.combined <-
  FindNeighbors(immune.combined, reduction = "pca", dims = 1:ndims)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <-
  DimPlot(immune.combined,
          reduction = "umap", group.by = "orig.ident")
p2 <-
  DimPlot(
    immune.combined,
    reduction = "umap",
    label = TRUE,
    repel = TRUE
  )
p1 + p2

DimPlot(immune.combined,
        reduction = "umap", split.by = "orig.ident")
saveRDS(immune.combined, file = "./output/seurat/ifnb.rds")
#### Identify conserved cell type markers ----
# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(immune.combined) <- "RNA"
nk.markers <-
  FindConservedMarkers(
    immune.combined,
    ident.1 = 6,
    grouping.var = "orig.ident",
    verbose = FALSE
  )
head(nk.markers)

FeaturePlot(
  immune.combined,
  features = c(
    "CD3D",
    "SELL",
    "CREM",
    "CD8A",
    "GNLY",
    "CD79A",
    "FCGR3A",
    "CCL2",
    "PPBP"
  ),
  min.cutoff = "q9"
)

immune.combined <-
  RenameIdents(
    immune.combined,
    `0` = "CD14 Mono",
    `1` = "CD4 Naive T",
    `2` = "CD4 Memory T",
    `3` = "CD16 Mono",
    `4` = "B",
    `5` = "CD8 T",
    `6` = "NK",
    `7` = "T activated",
    `8` = "DC",
    `9` = "B Activated",
    `10` = "Mk",
    `11` = "pDC",
    `12` = "Eryth",
    `13` = "Mono/Mk Doublets",
    `14` = "HSPC"
  )
DimPlot(immune.combined, label = TRUE)

Idents(immune.combined) <-
  factor(
    Idents(immune.combined),
    levels = c(
      "HSPC",
      "Mono/Mk Doublets",
      "pDC",
      "Eryth",
      "Mk",
      "DC",
      "CD14 Mono",
      "CD16 Mono",
      "B Activated",
      "B",
      "CD8 T",
      "NK",
      "T activated",
      "CD4 Naive T",
      "CD4 Memory T"
    )
  )
markers.to.plot <-
  c(
    "CD3D",
    "CREM",
    "HSPH1",
    "SELL",
    "GIMAP5",
    "CACYBP",
    "GNLY",
    "NKG7",
    "CCL5",
    "CD8A",
    "MS4A1",
    "CD79A",
    "MIR155HG",
    "NME1",
    "FCGR3A",
    "VMO1",
    "CCL2",
    "S100A9",
    "HLA-DQA1",
    "GPR183",
    "PPBP",
    "GNG11",
    "HBA2",
    "HBB",
    "TSPAN13",
    "IL3RA",
    "IGJ",
    "PRSS57"
  )
DotPlot(
  immune.combined,
  features = markers.to.plot,
  cols = c("blue", "red"),
  dot.scale = 8,
  split.by = "stim"
) +
  RotatedAxis()

#### Identify differential expressed genes across conditions ----
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <-
  as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <-
  as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ISG15",
                   "LY6E",
                   "IFI6",
                   "ISG20",
                   "MX1",
                   "IFIT2",
                   "IFIT1",
                   "CXCL10",
                   "CCL8")
p1 <-
  ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1,
                  points = genes.to.label,
                  repel = TRUE)
p2 <-
  ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2,
                  points = genes.to.label,
                  repel = TRUE)
p1 + p2

immune.combined$celltype.stim <-
  paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
b.interferon.response <-
  FindMarkers(
    immune.combined,
    ident.1 = "B_STIM",
    ident.2 = "B_CTRL",
    verbose = FALSE
  )
head(b.interferon.response, n = 15)

FeaturePlot(
  immune.combined,
  features = c("CD3D", "GNLY", "IFI6"),
  split.by = "stim",
  max.cutoff = 3,
  cols = c("grey", "red")
)

plots <-
  VlnPlot(
    immune.combined,
    features = c("LYZ", "ISG15", "CXCL10"),
    split.by = "stim",
    group.by = "celltype",
    pt.size = 0,
    combine = FALSE
  )
wrap_plots(plots = plots, ncol = 1)
