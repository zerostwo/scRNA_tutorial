#### Information ----
# Title   :   Birthday color table
# File    :   color.R
# Author  :   Songqi Duan
# Contact :   songqi.duan@outlook.com
# License :   Copyright (C) by Songqi Duan
# Created :   2021/06/26 18:08:10
# Updated :   none

color <- c(
  # 朋友
  "#83ccd2",
  # Fei Xie, 0813
  "#f19ca7",
  # Renshu Li, 0627
  "#0098d4",
  # Liujia Yuan, 0802
  "#fbbaa8",
  # Xinyu Tao, 1211
  "#bf783e",
  # Shiteng Yang, 0914
  # 博二
  "#b4766b",
  # Lan Zhang, 1118
  # 博一
  "#0f2350",
  # Tingting Tang, 0125
  "#665a1a",
  # Yuxian You, 0608
  # 研三
  "#4c5e74",
  # Yiwen Li, 0217
  "#caa980",
  # Yansu Yan, 1202
  "#f6ad49",
  # Yuan Yuan, 0818
  # 研二
  "#473a17",
  # Weiming Huang, 1115
  "#3f61a1",
  # Xiaoyu Duan, 0810
  "#009854",
  # Linjian Li, 0502
  # 研一
  "#bb5520",
  # Yu Ke, 1015
  "#fdb86d",
  # Hangyan Dan, 0925
  "#d0af4c",
  # Xuejiao Li, 1023
  "#ddeab2",
  # Feng Zhou, 0605
  "#b16268",
  # Yixi Liu, 1215
  "#288c66",
  # Yaowen Tan, 0120
  # 大四
  "#264939",
  # Songqi Duan, 0524
  # 大三
  "#6b76ae",
  # Mingyue Zhang, 0717
  "#b370a6",
  # Jiaqi Ding, 0715
  "#409ecc",
  # Chengyu Bai, 0803
  # 大二
  "#f7b894",
  # Zepeng Gu, 0615
  "#94adda",
  # Wanyun Han, 0716
  "#bea2ca" # Changyi Zhang, 0601
  
)

scales::show_col(color, borders = NA, labels = F)
scales::show_col(sample(color,
                        length(color)), borders = NA, labels = F)
