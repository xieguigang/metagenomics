# 加载必要的R包
required_packages <- c("vegan", "ape", "ggplot2", "dplyr", "tidyr", "reshape2",
                       "pairwiseAdonis", "microbiome", "microbiomeMarker",
                       "ANCOMBC", "ALDEx2", "igraph", "ggpubr", "GUniFrac", "VennDiagram",
                       "proxy", "ggraph", "tidygraph", "ggrepel", "corrplot", "RColorBrewer");

# 检查并安装缺失的包
for(pkg in required_packages) {
    library(pkg, character.only = TRUE);
}

pkg_dir = "G:\\metagenomics\\data\\16s";

import = function(h) {
    if (tolower(tools::file_ext(h)) == "r") {
        # is script
        source(file.path(pkg_dir, h));
    } else {
        # is folder
        for(file in list.files(file.path(pkg_dir,h), pattern = "\\.[rR]$", full.names = TRUE, recursive = FALSE)) {
            source(file);
        }
    }
}

import("./process_analysis.R");
import("./export_and_visualize_results.R");
import("./generate_summary_report.R");

import("./tools/");
import("./analyze/");
import("./export/");