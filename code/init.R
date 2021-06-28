#### Information ----
# Title   :   Initialize the project file 
# File    :   init.R
# Author  :   Songqi Duan
# Contact :   songqi.duan@outlook.com
# License :   Copyright (C) by Songqi Duan
# Created :   2021/06/26 17:22:34
# Updated :   none

# 创建保存Seurat对象文件夹
dir.create("./output/seurat")

# 创建保存细胞通讯文件夹
dir.create("./output/cellphonedb")

# 创建保存差异基因文件夹
dir.create("./output/deg")

# 创建保存功能富集文件夹
dir.create("./output/function")
dir.create("./output/function/go")
dir.create("./output/function/kegg")

