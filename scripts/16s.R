# ==============================================================================
# 16S Microbiome Comprehensive Analysis Pipeline
# Date: 2025-11-02
# Input: otu_table.txt (SampleID, Group, OTUs...)
# Optional Input: clinical_data.xlsx, PICRUSt2 outputs
# Output: A series of PDF/PNG figures and XLSX result tables.
# ==============================================================================


# --- 0. Environment Setup, Library Installation, and Data Loading ---

# 0.1 Clean environment
rm(list = ls())

# 0.2 Install and load required packages
# Create a list of packages to be installed/loaded
packages <- c(
  "tidyverse", "openxlsx", "phyloseq", "microbiome", "vegan", "picante",
  "DESeq2", "randomForest", "pROC", "VennDiagram", "igraph", "ggraph",
  "Hmisc", "corrplot", "pheatmap", "ape", "ggtree", "ggpubr", "RColorBrewer"
)

# Check if packages are installed, install if not, then load them
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  library(pkg, character.only = TRUE)
}

# Bioconductor packages
bioc_packages <- c("phyloseq", "microbiome", "DESeq2")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    }
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}


# 0.3 Set up output directories
main_dir <- "Microbiome_Analysis_Results"
dir.create(main_dir, showWarnings = FALSE)
sub_dirs <- c("01_Community_Structure", "02_Diversity", "03_Differential_Analysis", 
              "04_ROC_Analysis", "05_Network_Analysis", "06_Clinical_Correlation",
              "07_Function_Prediction")
sapply(sub_dirs, function(dir) dir.create(file.path(main_dir, dir), showWarnings = FALSE))

# 0.4 Load and prepare data
# Load OTU table
otu_data <- read.csv("GROUP2/species.csv", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1);
otu_data$SampleID <- NULL
meta_data <- as.data.frame(otu_data$class)
colnames(meta_data) <- "Group"
rownames(meta_data) <- rownames(otu_data)
otu_data$class <- NULL

# Create OTU matrix for phyloseq
otu_matrix <- as.matrix(otu_data)

# Create taxonomy table from OTU names (ASSUMPTION: names are like k__...;p__...;...)
parse_taxonomy <- function(tax_str) {
  tax_vec <- unlist(strsplit(tax_str, ";"))
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_df <- data.frame(matrix(NA, nrow = 1, ncol = length(tax_levels)))
  colnames(tax_df) <- tax_levels
  for (i in 1:length(tax_vec)) {
    if (grepl("__", tax_vec[i])) {
      tax_df[1, i] <- gsub("^[kpcofgs]__", "", tax_vec[i])
    }
  }
  return(tax_df)
}
tax_list <- apply(data.frame(OTU = colnames(otu_matrix)), 1, function(x) parse_taxonomy(x))
tax_table <- do.call(rbind, tax_list)
rownames(tax_table) <- colnames(otu_matrix)
tax_table <- as.matrix(tax_table)

# Create phyloseq object
PS <- phyloseq(otu_table(otu_matrix, taxa_are_rows = FALSE),
               sample_data(meta_data),
               tax_table(tax_table))

# Filter out low-prevalence taxa (e.g., present in less than 10% of samples)
PS <- filter_taxa(PS, function(x) sum(x > 0) > (0.1 * length(x)), TRUE)

# --- 1. Community Structure Analysis ---

# 1.1 Community Composition (Stacked Bar Plot)
plot_composition <- function(ps, level) {
  ps_glom <- tax_glom(ps, taxrank = level)
  ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
  df <- psmelt(ps_rel)
  df$Abundance <- df$Abundance * 100
  
  p <- ggplot(df, aes(x = Sample, y = Abundance, fill = get(level))) +
    geom_bar(stat = "identity", width = 0.7) +
    facet_wrap(~Group, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(x = "Samples", y = "Relative Abundance (%)", fill = level)
  
  # Save plots
  ggsave(filename = file.path(main_dir, "01_Community_Structure", paste0("Composition_", level, ".pdf")), plot = p, width = 12, height = 8)
  ggsave(filename = file.path(main_dir, "01_Community_Structure", paste0("Composition_", level, ".png")), plot = p, width = 12, height = 8, dpi = 300)
  
  # Save data
  write.xlsx(df, file = file.path(main_dir, "01_Community_Structure", paste0("Composition_", level, "_data.xlsx")))
}

for (lvl in c("Phylum", "Class", "Order", "Family", "Genus")) {
  if (lvl %in% rank_names(PS)) {
    plot_composition(PS, lvl)
  }
}

# 1.2 Venn Diagram (Shared/Unique OTUs)
group_list <- split(rownames(otu_table(PS)), sample_data(PS)$Group)
otu_lists <- lapply(group_list, function(g) {
  taxa_names(PS)[which(rowSums(otu_table(PS)[g, ]) > 0)]
})
venn.plot <- venn.diagram(
  x = otu_lists,
  category.names = names(otu_lists),
  filename = NULL,
  output = TRUE,
  col = "transparent",
  fill = RColorBrewer::brewer.pal(length(otu_lists), "Set3"),
  alpha = 0.50,
  cex = 1.5,
  cat.cex = 1.2,
  cat.col = "black",
  margin = 0.2
)
pdf(file.path(main_dir, "01_Community_Structure", "Venn_Diagram.pdf"), width = 8, height = 8)
grid.draw(venn.plot)
dev.off()
png(file.path(main_dir, "01_Community_Structure", "Venn_Diagram.png"), width = 2400, height = 2400, res = 300)
grid.draw(venn.plot)
dev.off()
write.xlsx(otu_lists, file.path(main_dir, "01_Community_Structure", "Venn_OTU_lists.xlsx"))

# 1.3 Sample-Species Heatmap (Top abundant species)
ps_rel <- transform_sample_counts(PS, function(x) x / sum(x))
# Select top N abundant genera (or other level)
top_taxa <- names(sort(taxa_sums(ps_rel), TRUE))[1:30]
ps_top <- prune_taxa(top_taxa, ps_rel)

# Prepare data for pheatmap
df_heatmap <- psmelt(ps_top)

# --- 关键修正点 1: 简化数据框以避免 dcast 的问题 ---
# 只保留构建矩阵必需的列，防止其他分类学信息被 dcast 误用为 id.vars
df_heatmap_simplified <- df_heatmap[, c("Sample", "Genus", "Abundance")]

# 将长格式数据转换为宽格式矩阵
df_heatmap_cast <- dcast(df_heatmap_simplified, Sample ~ Genus, value.var = "Abundance", fill = 0)
rownames(df_heatmap_cast) <- df_heatmap_cast$Sample
df_heatmap_cast$Sample <- NULL

# --- 关键修正点 2: 确保矩阵是数值型 ---
# as.matrix 在 data.frame 包含非数值列时会将整个矩阵转为字符型
# 现在由于 df_heatmap_cast 只包含数值列，可以安全转换
mat_heatmap <- as.matrix(df_heatmap_cast)

# --- 关键修正点 3: 准备注释数据并修正参数 ---
# 创建一个干净的注释数据框
annotation_df <- as.data.frame(sample_data(PS))
# 确保行名与热图矩阵的行名（样本名）匹配
rownames(annotation_df) <- sample_names(PS)
# 只保留分组信息列，防止其他列干扰
annotation_df <- annotation_df[, "Group", drop = FALSE]

# 绘制热图
# 注意：因为样本是行，所以应该使用 annotation_row 而不是 annotation_col
pheatmap(mat_heatmap,
         annotation_row = annotation_df,
         show_rownames = FALSE,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         filename = file.path(main_dir, "01_Community_Structure", "Sample_Species_Heatmap.pdf"),
         width = 10, height = 12)

# PNG版本
pheatmap(mat_heatmap,
         annotation_row = annotation_df,
         show_rownames = FALSE,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         filename = file.path(main_dir, "01_Community_Structure", "Sample_Species_Heatmap.png"),
         width = 10, height = 12)

# 保存用于绘图的数据
write.xlsx(df_heatmap, file = file.path(main_dir, "01_Community_Structure", "Heatmap_data.xlsx"))



# --- 2. Diversity Analysis ---

# 2.1 Alpha Diversity
alpha_div <- estimate_richness(PS, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
alpha_div$Group <- sample_data(PS)$Group
alpha_div$SampleID <- rownames(alpha_div)

# Good's Coverage
goods_coverage <- function(x) { 1 - sum(x == 1) / sum(x) }
coverage <- apply(otu_table(PS), 1, goods_coverage)
alpha_div$GoodsCoverage <- coverage[match(rownames(alpha_div), names(coverage))]

# Save alpha diversity table
write.xlsx(alpha_div, file = file.path(main_dir, "02_Diversity", "Alpha_Diversity_Table.xlsx"))

# Plot alpha diversity
plot_alpha <- function(measure) {
  p <- ggplot(alpha_div, aes(x = Group, y = get(measure), fill = Group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 1.5) +
    theme_bw() +
    stat_compare_means(method = "kruskal.test", label = "p.signif") +
    labs(x = "Group", y = measure, title = paste("Alpha Diversity:", measure))
  
  ggsave(file.path(main_dir, "02_Diversity", paste0("Alpha_", measure, ".pdf")), plot = p, width = 8, height = 6)
  ggsave(file.path(main_dir, "02_Diversity", paste0("Alpha_", measure, ".png")), plot = p, width = 8, height = 6, dpi = 300)
}
for (m in c("Observed", "Chao1", "Shannon", "Simpson", "GoodsCoverage")) {
  plot_alpha(m)
}

# 2.2 Beta Diversity
# Calculate distance matrix
dist_bray <- distance(PS, method = "bray")
dist_unifrac <- distance(PS, method = "unifrac", weighted = TRUE)

# Ordination functions
plot_ordination_ps <- function(ordination_method, dist_method, title) {
  ord <- ordinate(PS, method = ordination_method, distance = dist_method)
  p <- plot_ordination(PS, ord, color = "Group", shape = "Group") +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = title)
  
  ggsave(file.path(main_dir, "02_Diversity", paste0(ordination_method, "_", dist_method, ".pdf")), plot = p, width = 8, height = 6)
  ggsave(file.path(main_dir, "02_Diversity", paste0(ordination_method, "_dist_method", ".png")), plot = p, width = 8, height = 6, dpi = 300)
}

# PCA
plot_ordination_ps("PCA", "bray", "PCA (Bray-Curtis)")
# PCoA
plot_ordination_ps("PCoA", "bray", "PCoA (Bray-Curtis)")
plot_ordination_ps("PCoA", "unifrac", "PCoA (Weighted UniFrac)")
# NMDS
plot_ordination_ps("NMDS", "bray", "NMDS (Bray-Curtis)")

# Statistical tests for beta diversity
group_factor <- sample_data(PS)$Group
anosim_res <- anosim(dist_bray, group_factor)
adonis_res <- adonis2(dist_bray ~ group_factor)

beta_stats <- rbind(
  ANOSIM_R = c(R = anosim_res$statistic, `P-value` = anosim_res$signif),
  Adonis_R2 = c(R = adonis_res$R2[1], `P-value` = adonis_res$`Pr(>F)`[1])
)
write.xlsx(beta_stats, file = file.path(main_dir, "02_Diversity", "Beta_Diversity_Statistics.xlsx"))

# UPGMA Clustering
hc <- hclust(dist_bray, method = "average")
plot(hc, main = "UPGMA Clustering (Bray-Curtis)", xlab = "", sub = "")
dev.copy(pdf, file.path(main_dir, "02_Diversity", "UPGMA_Tree.pdf"), width = 8, height = 6)
dev.off()
dev.copy(png, file.path(main_dir, "02_Diversity", "UPGMA_Tree.png"), width = 2400, height = 1800, res = 300)
dev.off()


# --- 3. Differential Abundance Analysis ---

# 3.1 LEfSe Analysis (using lem package)
# Note: lem requires specific formatting. We provide the phyloseq object.
# It's good to aggregate to a certain taxonomic level, e.g., Genus.
PS_genus <- tax_glom(PS, "Genus")
lem_res <- lem(PS_genus, group = "Group", wilcoxon_cutoff = 0.05)
# lem_res is a list, we are interested in the LDA results
lda_res <- lem_res$lda_res
lda_res <- lda_res[order(lda_res$LDA, decreasing = TRUE), ]
write.xlsx(lda_res, file = file.path(main_dir, "03_Differential_Analysis", "LEfSe_Results.xlsx"))

# Plot LEfSe results (Cladogram is complex, a bar plot of LDA scores is common)
p_lda <- ggplot(lda_res[1:20, ], aes(x = reorder(Taxon, LDA), y = LDA, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  labs(x = "Taxon", y = "LDA Score", title = "LEfSe LDA Scores (Top 20)")
ggsave(file.path(main_dir, "03_Differential_Analysis", "LEfSe_LDA_Plot.pdf"), plot = p_lda, width = 10, height = 8)
ggsave(file.path(main_dir, "03_Differential_Analysis", "LEfSe_LDA_Plot.png"), plot = p_lda, width = 10, height = 8, dpi = 300)


# 3.2 Group-wise significance analysis (using DESeq2)
# DESeq2 requires raw counts, not normalized. We use the original PS object.
diagdds <- phyloseq_to_deseq2(PS, ~ Group)
# Calculate geometric means for normalization
gm_mean = function(x) {
  if (all(x == 0)) {
    return(0)
  }
  exp(mean(log(x[x > 0])))
}
geo_means <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geo_means)
diagdds <- DESeq(diagdds, test = "Wald", fitType = "local")

# Results for pairwise comparisons
groups <- unique(sample_data(PS)$Group)
pairwise_comb <- combn(groups, 2, simplify = FALSE)
res_list <- list()

for (p in pairwise_comb) {
  contrast <- c("Group", p[1], p[2])
  res <- results(diagdds, contrast = contrast, alpha = 0.05)
  res <- res[order(res$padj), ]
  sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
  res_list[[paste(p[1], "_vs_", p[2], sep = "")]] <- as.data.frame(sig_res)
}
write.xlsx(res_list, file = file.path(main_dir, "03_Differential_Analysis", "DESeq2_Results.xlsx"))

# 3.3 Random Forest Analysis
# Convert to matrix suitable for randomForest
otu_mat <- t(otu_table(PS))
group_vec <- as.factor(sample_data(PS)$Group)

set.seed(123)
rf_model <- randomForest(x = otu_mat, y = group_vec, importance = TRUE, ntree = 500)
rf_importance <- importance(rf_model)
rf_importance_df <- data.frame(Taxon = rownames(rf_importance), MeanDecreaseGini = rf_importance[, "MeanDecreaseGini"])
rf_importance_df <- rf_importance_df[order(rf_importance_df$MeanDecreaseGini, decreasing = TRUE), ]
write.xlsx(rf_importance_df, file = file.path(main_dir, "03_Differential_Analysis", "RandomForest_Importance.xlsx"))

p_rf <- ggplot(rf_importance_df[1:20, ], aes(x = reorder(Taxon, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(x = "Taxon (OTU)", y = "Mean Decrease in Gini", title = "Random Forest Importance (Top 20)")
ggsave(file.path(main_dir, "03_Differential_Analysis", "RandomForest_Plot.pdf"), plot = p_rf, width = 10, height = 8)
ggsave(file.path(main_dir, "03_Differential_Analysis", "RandomForest_Plot.png"), plot = p_rf, width = 10, height = 8, dpi = 300)


# --- 4. Disease Diagnosis Model (ROC Curve Analysis) ---
# Note: ROC is for binary classification. We will perform pairwise comparisons.
# We use a top biomarker from LEfSe or Random Forest.
biomarker <- rf_importance_df$Taxon[1]

roc_list <- list()
for (p in pairwise_comb) {
  group1_samples <- rownames(sample_data(PS))[sample_data(PS)$Group == p[1]]
  group2_samples <- rownames(sample_data(PS))[sample_data(PS)$Group == p[2]]
  
  df_roc <- data.frame(
    Abundance = otu_mat[c(group1_samples, group2_samples), biomarker],
    Group = factor(c(rep(p[1], length(group1_samples)), rep(p[2], length(group2_samples))))
  )
  
  roc_obj <- roc(df_roc$Group, df_roc$Abundance, levels = rev(levels(df_roc$Group)), direction = "<")
  roc_list[[paste(p[1], "_vs_", p[2], sep = "")]] <- list(auc = auc(roc_obj), ci = ci.auc(roc_obj))
  
  p_roc <- ggroc(roc_obj) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype = "dashed") +
    theme_bw() +
    labs(title = paste("ROC Curve for", biomarker, ":", p[1], "vs", p[2]),
         subtitle = paste("AUC =", round(auc(roc_obj), 3)))
  
  ggsave(filename = file.path(main_dir, "04_ROC_Analysis", paste0("ROC_", p[1], "_vs_", p[2], ".pdf")), plot = p_roc, width = 6, height = 5)
  ggsave(filename = file.path(main_dir, "04_ROC_Analysis", paste0("ROC_", p[1], "_vs_", p[2], ".png")), plot = p_roc, width = 6, height = 5, dpi = 300)
}
write.xlsx(roc_list, file = file.path(main_dir, "04_ROC_Analysis", "ROC_AUC_Results.xlsx"))


# --- 5. Species Correlation Network Analysis ---
# We'll build a network based on Spearman correlation.
# Filter to a specific group or use all samples. Here we use all.
cor_mat <- rcorr(as.matrix(otu_mat), type = "spearman")
r_values <- cor_mat$r
p_values <- cor_mat$P

# Set thresholds
r_threshold <- 0.6
p_threshold <- 0.05

# Create edge list
edges <- which(p_values < p_threshold & abs(r_values) > r_threshold, arr.ind = TRUE)
edges <- edges[edges[, 1] < edges[, 2], ] # Remove duplicates
edge_list <- data.frame(
  from = rownames(r_values)[edges[, 1]],
  to = colnames(r_values)[edges[, 2]],
  weight = r_values[edges]
)
write.xlsx(edge_list, file = file.path(main_dir, "05_Network_Analysis", "Network_Edges.xlsx"))

# Create igraph object
net <- graph_from_data_frame(edge_list, directed = FALSE)
net <- simplify(net) # Remove loops and multiple edges

# Add node attributes (e.g., Phylum)
node_info <- data.frame(Taxon = V(net)$name)
node_info$Phylum <- tax_table(PS)[node_info$Taxon, "Phylum"]
V(net)$Phylum <- node_info$Phylum

# Plot network
ggraph(net, layout = "fr") +
  geom_edge_link(aes(width = abs(weight), color = weight > 0), alpha = 0.5) +
  scale_edge_color_manual(values = c("red", "blue"), guide = "none") +
  geom_node_point(aes(color = Phylum), size = 5) +
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Spearman Correlation Network")

ggsave(file.path(main_dir, "05_Network_Analysis", "Correlation_Network.pdf"), width = 10, height = 10)
ggsave(file.path(main_dir, "05_Network_Analysis", "Correlation_Network.png"), width = 10, height = 10, dpi = 300)


# --- 6. Correlation with Clinical Indicators ---
# IMPORTANT: This section requires a 'clinical_data.xlsx' file.
if (file.exists("clinical_data.xlsx")) {
  clinical_data <- read.xlsx("clinical_data.xlsx")
  # Ensure SampleID is the first column and is character
  clinical_data$SampleID <- as.character(clinical_data[[1]])
  rownames(clinical_data) <- clinical_data$SampleID
  clinical_data$SampleID <- NULL
  
  # Merge with phyloseq sample data
  sample_data(PS) <- merge(sample_data(PS), clinical_data, by = "row.names")
  colnames(sample_data(PS))[1] <- "SampleID"
  rownames(sample_data(PS)) <- sample_data(PS)$SampleID
  sample_data(PS)$SampleID <- NULL
  
  # 6.1 Spearman correlation between species (e.g., Genus) and clinical factors
  ps_genus_rel <- transform_sample_counts(tax_glom(PS, "Genus"), function(x) x / sum(x))
  genus_mat <- t(otu_table(ps_genus_rel))
  clinical_mat <- as.matrix(sample_data(PS)[, !(colnames(sample_data(PS)) %in% "Group")])
  
  # Ensure samples match
  common_samples <- intersect(rownames(genus_mat), rownames(clinical_mat))
  genus_mat <- genus_mat[common_samples, ]
  clinical_mat <- clinical_mat[common_samples, ]
  
  cor_res <- rcorr(genus_mat, clinical_mat, type = "spearman")
  
  # Plot correlation heatmap
  corrplot(cor_res$r, type = "upper", tl.col = "black", tl.srt = 45, 
           p.mat = cor_res$P, sig.level = 0.05, insig = "blank",
           title = "Spearman Correlation: Genus vs Clinical", mar = c(0,0,1,0))
  dev.copy(pdf, file.path(main_dir, "06_Clinical_Correlation", "Spearman_Correlation_Heatmap.pdf"), width = 10, height = 10)
  dev.off()
  dev.copy(png, file.path(main_dir, "06_Clinical_Correlation", "Spearman_Correlation_Heatmap.png"), width = 3000, height = 3000, res = 300)
  dev.off()
  
  # Save correlation results
  cor_r_df <- as.data.frame(cor_res$r)
  cor_p_df <- as.data.frame(cor_res$P)
  write.xlsx(list(`Correlation_r` = cor_r_df, `P_value` = cor_p_df), 
             file = file.path(main_dir, "06_Clinical_Correlation", "Spearman_Correlation_Results.xlsx"))
  
  # 6.2 Environmental Factor Analysis (RDA/CCA)
  # Use Hellinger transformation for species data
  genus_hel <- decostand(genus_mat, method = "hellinger")
  
  # Perform RDA
  rda_res <- rda(genus_hel ~ ., data = as.data.frame(clinical_mat))
  
  # Plot RDA
  plot(rda_res, type = "text")
  dev.copy(pdf, file.path(main_dir, "06_Clinical_Correlation", "RDA_Plot.pdf"), width = 8, height = 8)
  dev.off()
  dev.copy(png, file.path(main_dir, "06_Clinical_Correlation", "RDA_Plot.png"), width = 2400, height = 2400, res = 300)
  dev.off()
  
  # Save RDA results
  rda_summary <- summary(rda_res)
  write.xlsx(list(`Summary` = capture.output(rda_summary), `Scores` = scores(rda_res)), 
             file = file.path(main_dir, "06_Clinical_Correlation", "RDA_Results.xlsx"))
  
} else {
  message("Skipping Clinical Correlation Analysis: 'clinical_data.xlsx' not found.")
}


# --- 7. Functional Prediction (PICRUSt2) ---
# IMPORTANT: This section assumes PICRUSt2 has been run externally and outputs are available.
picrust_dir <- "picrust2_output"
if (dir.exists(picrust_dir)) {
  
  # Function to analyze PICRUSt2 output
  analyze_picrust <- function(filename, analysis_name) {
    filepath <- file.path(picrust_dir, filename)
    if (!file.exists(filepath)) {
      message(paste("File not found:", filename, "Skipping", analysis_name))
      return(NULL)
    }
    
    func_data <- read.delim(filepath, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
    # Ensure sample names match
    func_data <- func_data[, intersect(colnames(func_data), rownames(sample_data(PS)))]
    sample_data(PS)$Group <- factor(sample_data(PS)$Group) # ensure factor
    
    # Melt for plotting
    df_melt <- psmelt(phyloseq(otu_table(as.matrix(func_data), taxa_are_rows = TRUE), sample_data(PS)))
    
    # Plot top 20 abundant pathways
    top_pathways <- df_melt %>% group_by(OTU) %>% summarise(mean_ab = mean(Abundance)) %>% top_n(20, mean_ab) %>% pull(OTU)
    df_top <- df_melt %>% filter(OTU %in% top_pathways)
    
    p <- ggplot(df_top, aes(x = Group, y = Abundance, fill = Group)) +
      geom_boxplot() +
      facet_wrap(~OTU, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(filename = file.path(main_dir, "07_Function_Prediction", paste0(analysis_name, "_Top20_Boxplot.pdf")), plot = p, width = 12, height = 10)
    ggsave(filename = file.path(main_dir, "07_Function_Prediction", paste0(analysis_name, "_Top20_Boxplot.png")), plot = p, width = 12, height = 10, dpi = 300)
    
    # Differential analysis using DESeq2
    diagdds_func <- phyloseq_to_deseq2(phyloseq(otu_table(as.matrix(func_data), taxa_are_rows = TRUE), sample_data(PS)), ~ Group)
    geo_means_func <- apply(counts(diagdds_func), 1, gm_mean)
    diagdds_func <- estimateSizeFactors(diagdds_func, geoMeans = geo_means_func)
    diagdds_func <- DESeq(diagdds_func, test = "Wald", fitType = "local")
    
    res_list_func <- list()
    for (p in pairwise_comb) {
      contrast <- c("Group", p[1], p[2])
      res <- results(diagdds_func, contrast = contrast, alpha = 0.05)
      res <- res[order(res$padj), ]
      sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
      if (nrow(sig_res) > 0) {
        res_list_func[[paste(analysis_name, p[1], "_vs_", p[2], sep = "_")]] <- as.data.frame(sig_res)
      }
    }
    
    if (length(res_list_func) > 0) {
      write.xlsx(res_list_func, file = file.path(main_dir, "07_Function_Prediction", paste0(analysis_name, "_DESeq2_Results.xlsx")))
    }
  }
  
  analyze_picrust("pathabundance.tsv", "Pathways")
  analyze_picrust("ko_abundance.tsv", "KO")
  analyze_picrust("cog_abundance.tsv", "COG")
  analyze_picrust("ec_abundance.tsv", "EC")
  analyze_picrust("pfam_abundance.tsv", "PFAM")
  analyze_picrust("tigrfam_abundance.tsv", "TIGRFAM")
  
} else {
  message("Skipping Functional Prediction: 'picrust2_output' directory not found.")
}


# --- 8. Phenotype Prediction ---
# NOTE: Phenotype prediction (e.g., BugBase) is typically done via web servers.
# It cannot be directly performed in an R script.
# Please upload your OTU table to a service like BugBase (https://bugbase.cs.umn.edu/)
# to get phenotype predictions (e.g., Gram-positive, Gram-negative, Potential Pathogens, etc.).
# The output from BugBase can then be analyzed similarly to the functional prediction data above.
message("Phenotype Prediction (e.g., BugBase) requires external web tools and is not included in this script.")


# --- End of Script ---
message("All analyses completed. Results are saved in the 'Microbiome_Analysis_Results' directory.")

