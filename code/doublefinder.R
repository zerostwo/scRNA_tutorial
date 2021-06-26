#### 加载包 ----
library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)
rm(list = ls())

#### 导入Seurat对象 ----
seurat_obj <- readRDS("./output/seurat/seurat_obj.rds")
dims <- 10
#### 寻找最优pK值 ----
sweep.res.list <- paramSweep_v3(seurat_obj,
                                PCs = 1:dims,
                                sct = F,
                                num.cores = 10)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <-
  bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = 0.039                     # 5000细胞对应的doublets rate是3.9%
homotypic.prop <-
  modelHomotypic(seurat_obj$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate * ncol(seurat_obj))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

## 使用确定好的参数鉴定doublets
seurat_obj <-
  doubletFinder_v3(
    seurat_obj,
    PCs = 1:dims,
    pN = 0.25,
    pK = pK_bcmvn,
    nExp = nExp_poi.adj,
    reuse.pANN = F,
    sct = F
  )

## 结果展示，分类结果在seurat_obj@meta.data中
DimPlot(
  seurat_obj,
  reduction = "tsne",
  group.by = colnames(seurat_obj@meta.data)[grep("^DF", colnames(seurat_obj@meta.data))],
  order = T
)
table(seurat_obj@meta.data[, colnames(seurat_obj@meta.data)[grep("^DF", colnames(seurat_obj@meta.data))]]) ##查看双细胞及单细胞数量
seurat_obj <-
  seurat_obj[, seurat_obj@meta.data[, colnames(seurat_obj@meta.data)[grep("^DF", colnames(seurat_obj@meta.data))]] == "Singlet"]
