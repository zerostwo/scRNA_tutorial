#### Information ----
# Title   :   Seurat standard analysis process
# File    :   seurat_basic.R
# Author  :   Songqi Duan
# Contact :   songqi.duan@outlook.com
# License :   Copyright (C) by Songqi Duan
# Created :   2021/06/21 17:30:11
# Updated :   none

#### 导入包 ----
library(Seurat)
library(dplyr)
library(ggplot2)

#### 数据导入 ----
data_dir <- "./data/pbmc3k/" # 更换为自己的文件
list.files(data_dir)
# [1] "barcodes.tsv" "genes.tsv"    "matrix.mtx"，显示这三个文件才对
expr <- Read10X(data.dir = data_dir)
dim(expr) # 查看表达矩阵维度，行是基因，列是细胞

#### 创建Seurat对象 ----
seurat_obj <- CreateSeuratObject(
  counts = expr,
  project = "pbmc3k",
  # 样本名字，可根据实际情况更改
  min.cells = 3,
  # 如果一个基因在3个细胞以下表达，则过滤掉，一般默认
  min.features = 200 # 如果一个细胞表达200个基因以下，则过滤掉，一般默认
)
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# 这个警告信息是将基因名字里含有下划线(_)的替换成了减号(-)

#### 质控和选择细胞 ----
seurat_obj[["percent.mt"]] <-
  PercentageFeatureSet(seurat_obj, pattern = "^MT-") # 人源是MT,鼠源是mt，注意更换
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
seurat_obj <- subset(seurat_obj, subset =
                       nFeature_RNA > 200 &
                       nFeature_RNA < 2500 &
                       percent.mt < 5)
#### 标准化数据 ----
# 消除测序深度的影响
seurat_obj <-
  NormalizeData(seurat_obj,
                normalization.method = "LogNormalize",
                scale.factor = 10000)
# boxplot(as.matrix(seurat_obj@assays$RNA@counts)[,1:50])
# boxplot(as.matrix(seurat_obj@assays$RNA@data)[,1:50])

#### 高变基因的选择 ----
# 下游分析选择2000个高变基因进行分析，不但提升了分析速度而且对后续结果没影响
seurat_obj <- FindVariableFeatures(seurat_obj,
                                   selection.method = "vst",
                                   nfeatures = 2000)

#### 归一化数据 ----
# 归一化可以将每个细胞内基因的表达值调整到0的上下左右
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
# 如果细胞数量过多，可以使用高变基因进行归一化，这样速度更快
# seurat_obj <- ScaleData(seurat_obj)
# 如果想消去线粒体效应，可以使用以下指令
# seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "percent.mt")
# boxplot(as.matrix(seurat_obj@assays$RNA@scale.data)[,1:50])

#### PCA分析 ----
seurat_obj <- RunPCA(seurat_obj,
                     features = VariableFeatures(object = seurat_obj))
DimPlot(seurat_obj, reduction = "pca")

#### 确定数据集的“主成分个数” ----
ElbowPlot(seurat_obj, ndims = 50)
dims <- "10"

#### 细胞聚类 ----
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims)
seurat_obj <-
  FindClusters(seurat_obj, resolution = 0.5) # 解析度越大分群越多，反之亦然
table(seurat_obj@active.ident)

#### 非线性降维 ----
seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims)
DimPlot(seurat_obj, reduction = "umap", label = T) + NoLegend()

seurat_obj <- RunTSNE(seurat_obj, dims = 1:dims)
DimPlot(seurat_obj, reduction = "tsne", label = T) + NoLegend()

# DimPlot(seurat_obj, reduction = "pca", label = T) + NoLegend()

#### 细胞周期回归 ----
seurat_obj <- CellCycleScoring(seurat_obj,
                               g2m.features = cc.genes$g2m.genes,
                               s.features = cc.genes$s.genes)
DimPlot(seurat_obj,
        reduction = "tsne",
        label = T,
        group.by = "Phase") +
  NoLegend()

#### 差异分析 ----
# 单个Cluster的差异基因分析
cluster.markers <- FindMarkers(
  seurat_obj,
  ident.1 = 1,
  only.pos = TRUE,
  # 如果是TRUE，则只显示上调的差异基因
  min.pct = 0.25,
  # 如果是做GSEA，这里改成0
  logfc.threshold = 0.25 # 如果是做GSEA，这里改成0
)

# 所有Cluster的差异基因分析
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  # 如果是TRUE，则只显示上调的差异基因
  min.pct = 0.25,
  # 如果是做GSEA，这里改成0
  logfc.threshold = 0.25 # 如果是做GSEA，这里改成0
)
top.markers <-
  markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#### marker基因可视化 -----
# 散点图
FeaturePlot(seurat_obj,
            features = c("PPBP", "MS4A1", "CD3E", "LYZ"),
            order = T)
# 小提琴图
VlnPlot(seurat_obj,
        features = c("PPBP", "MS4A1", "CD3E", "LYZ"),
        ncol = 2)

# 气泡图
DotPlot(seurat_obj, features = c("PPBP", "MS4A1", "CD3E", "LYZ"))

# 山脊图
RidgePlot(seurat_obj, c("PPBP", "MS4A1", "CD3E", "LYZ"), ncol = 2)

# 热图
DoHeatmap(seurat_obj, features = top.markers$gene)

#### 注释 ----
new.cluster.ids <- c(
  "T cells",
  # 00
  "Myeloid cells",
  # 01
  "T cells",
  # 02
  "B cells",
  # 03
  "T cells",
  # 04
  "Myeloid cells",
  # 05
  "T cells",
  # 06
  "Myeloid cells",
  # 07
  "Platelet" # 08
)
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

DimPlot(seurat_obj, reduction = "umap", label = T) + NoLegend()
DimPlot(seurat_obj, reduction = "tsne", label = T) + NoLegend()

#### 保存数据 ----
saveRDS(seurat_obj, file = "./output/seurat/seurat_obj.rds")

#### 亚群重聚类 ----
sub_seurat_obj <- subset(seurat_obj, idents = c("Myeloid cells"))
sub_seurat_obj <- FindVariableFeatures(sub_seurat_obj,
                                       selection.method = "vst",
                                       nfeatures = 2000)
all.genes <- rownames(sub_seurat_obj)
sub_seurat_obj <- ScaleData(sub_seurat_obj, features = all.genes)
sub_seurat_obj <- RunPCA(sub_seurat_obj,
                         features = VariableFeatures(sub_seurat_obj))
ElbowPlot(sub_seurat_obj, ndims = 50)
ndims <- 10
sub_seurat_obj <- FindNeighbors(sub_seurat_obj, dims = 1:ndims)
sub_seurat_obj <- FindClusters(sub_seurat_obj, resolution = 0.2)
sub_seurat_obj <- RunUMAP(sub_seurat_obj, dims = 1:ndims)
sub_seurat_obj <- RunTSNE(sub_seurat_obj, dims = 1:ndims)
DimPlot(sub_seurat_obj, reduction = "umap")

FeaturePlot(sub_seurat_obj,
            features = c("FCER1A", "FCGR3A", "CD14"))

new.cluster.ids <- c("CD14+ Mono", # 00
                     "FCGR3A+ Mono", # 01
                     "DC" # 02)
                     names(new.cluster.ids) <- levels(sub_seurat_obj)
                     sub_seurat_obj <- RenameIdents(sub_seurat_obj, new.cluster.ids)
                     
                     DimPlot(sub_seurat_obj, reduction = "umap", label = T) + NoLegend()
                     saveRDS(sub_seurat_obj, file = "./output/seurat/myeloid.rds")
                     