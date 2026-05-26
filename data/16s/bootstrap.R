# 加载必要的R包
required_packages <- c("vegan", "ape", "ggplot2", "dplyr", "tidyr", "reshape2",
                       "pairwiseAdonis", "microbiome", "microbiomeMarker",
                       "ANCOMBC", "ALDEx2", "igraph", "ggpubr", "GUniFrac", "VennDiagram",
                       "proxy", "ggraph", "tidygraph", "ggrepel", "corrplot", "RColorBrewer");

# 检查并安装缺失的包
for(pkg in required_packages) {
    library(pkg, character.only = TRUE);
}

