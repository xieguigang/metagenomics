# 加载必要的包
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(microbiome)
library(DESeq2)
library(pheatmap)
library(reshape2)
library(openxlsx)
library(multcomp)
library(ggpubr)
library(RColorBrewer)

# 创建输出目录
dir.create("16S_Results", showWarnings = FALSE)
dir.create("16S_Results/Plots", showWarnings = FALSE)
dir.create("16S_Results/Tables", showWarnings = FALSE)

# 1. 数据导入与预处理
# --------------------
# 读取OTU表格
otu_data <- read.csv("16s.csv", row.names = 1, check.names = FALSE)

# 提取样本信息
sample_data <- otu_data[, 1, drop = FALSE]
colnames(sample_data) <- "Group"
sample_data$Sample <- rownames(sample_data)

# 提取OTU矩阵
otu_matrix <- as.matrix(otu_data[, -1])
rownames(otu_matrix) <- rownames(otu_data)

# 解析物种注释信息
taxa_names <- colnames(otu_matrix)
tax_table <- data.frame(
  OTU = taxa_names,
  Kingdom = sapply(strsplit(taxa_names, ";"), `[`, 1),
  Phylum = sapply(strsplit(taxa_names, ";"), `[`, 2),
  Class = sapply(strsplit(taxa_names, ";"), `[`, 3),
  Order = sapply(strsplit(taxa_names, ";"), `[`, 4),
  Family = sapply(strsplit(taxa_names, ";"), `[`, 5),
  Genus = sapply(strsplit(taxa_names, ";"), `[`, 6),
  Species = sapply(strsplit(taxa_names, ";"), `[`, 7),
  stringsAsFactors = FALSE
)

# 清理分类信息
clean_tax <- function(x) {
  gsub("^[a-z]__", "", x)
}
tax_table[, -1] <- lapply(tax_table[, -1], clean_tax)
rownames(tax_table) <- taxa_names

# 创建phyloseq对象
OTU <- otu_table(otu_matrix, taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(tax_table))
SAM <- sample_data(sample_data)
physeq <- phyloseq(OTU, TAX, SAM)

# 过滤低丰度OTU (在至少20%样本中相对丰度>0.01%)
physeq <- filter_taxa(physeq, function(x) sum(x > 0) > length(x)*0.2, TRUE)

# 2. Alpha多样性分析
# --------------------
# 计算多样性指数
alpha_div <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$Group <- sample_data(physeq)$Group
alpha_div$Sample <- sample_names(physeq)

# 保存表格
write.xlsx(alpha_div, "16S_Results/Tables/Alpha_Diversity.xlsx")

# 绘制多样性指数箱线图
plot_alpha <- function(metric) {
  ggplot(alpha_div, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2) +
    stat_compare_means(method = "wilcox.test") +
    labs(title = paste("Alpha Diversity:", metric),
         y = metric,
         x = "Group") +
    theme_classic() +
    scale_fill_brewer(palette = "Set2")
}

metrics <- c("Observed", "Shannon", "Simpson")
lapply(metrics, function(m) {
  p <- plot_alpha(m)
  ggsave(paste0("16S_Results/Plots/Alpha_", m, ".png"), p, width = 8, height = 6, dpi = 300)
  ggsave(paste0("16S_Results/Plots/Alpha_", m, ".pdf"), p, width = 8, height = 6)
})

# 3. Beta多样性分析
# --------------------
# 转换为相对丰度
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# 计算Bray-Curtis距离
dist_bc <- phyloseq::distance(physeq_rel, method = "bray")

# PCoA分析
ord <- ordinate(physeq_rel, method = "PCoA", distance = dist_bc)

# 提取PCoA坐标
pcoa_df <- data.frame(sample_data(physeq),
                      scores(ord, display = "sites"))

# 绘制PCoA图
p_pcoa <- ggplot(pcoa_df, aes(Axis.1, Axis.2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  labs(title = "PCoA (Bray-Curtis)",
       x = paste0("PCoA1 (", round(ord$values$Relative_eig[1]*100, 1), "%)"),
       y = paste0("PCoA2 (", round(ord$values$Relative_eig[2]*100, 1), "%)")) +
  theme_classic() +
  scale_color_brewer(palette = "Set1")

ggsave("16S_Results/Plots/PCoA_Bray.png", p_pcoa, width = 8, height = 6, dpi = 300)
ggsave("16S_Results/Plots/PCoA_Bray.pdf", p_pcoa, width = 8, height = 6)

# PERMANOVA检验
adonis_res <- adonis2(dist_bc ~ Group, data = as(sample_data(physeq), "data.frame"))
write.xlsx(as.data.frame(adonis_res), "16S_Results/Tables/PERMANOVA.xlsx")

# 4. 物种组成分析
# --------------------
# 门水平组成
physeq_phylum <- tax_glom(physeq_rel, "Phylum")
phylum_rel <- transform_sample_counts(physeq_phylum, function(x) x / sum(x))
phylum_df <- psmelt(phylum_rel) %>%
  group_by(Sample, Group, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 绘制堆叠柱状图
p_phylum <- ggplot(phylum_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~Group, scales = "free_x") +
  labs(title = "Phylum Level Composition",
       y = "Relative Abundance") +
  theme_classic() +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("16S_Results/Plots/Phylum_Composition.png", p_phylum, width = 12, height = 6, dpi = 300)
ggsave("16S_Results/Plots/Phylum_Composition.pdf", p_phylum, width = 12, height = 6)

# 属水平热图
physeq_genus <- tax_glom(physeq_rel, "Genus")
genus_rel <- transform_sample_counts(physeq_genus, function(x) x / sum(x))
genus_mat <- otu_table(genus_rel)

# 选择前20个丰度最高的属
top_genus <- names(sort(colSums(genus_mat), decreasing = TRUE))[1:20]
genus_mat <- genus_mat[, top_genus]

# 绘制热图
pheatmap(genus_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = data.frame(Group = sample_data(physeq)$Group),
         show_rownames = FALSE,
         fontsize_col = 8,
         filename = "16S_Results/Plots/Genus_Heatmap.png",
         width = 10,
         height = 8)

# 5. 差异物种分析
# --------------------
# DESeq2差异分析
diagdds <- phyloseq_to_deseq2(physeq, ~ Group)
diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")

# 提取结果
res <- results(diagdds, contrast = c("Group", "OB", "CON"))
res_sig <- res[which(res$padj < 0.05), ]

# 保存结果
res_df <- as.data.frame(res)
write.xlsx(res_df, "16S_Results/Tables/DESeq2_Results.xlsx")

# 绘制火山图
volcano_df <- data.frame(
  log2FoldChange = res$log2FoldChange,
  padj = res$padj,
  sig = ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Yes", "No")
)

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Differential Abundance (OB vs CON)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_classic()

ggsave("16S_Results/Plots/Volcano_DESeq2.png", p_volcano, width = 8, height = 6, dpi = 300)
ggsave("16S_Results/Plots/Volcano_DESeq2.pdf", p_volcano, width = 8, height = 6)

# 6. 系统发育树（可选）
# --------------------
# 如果有系统发育树文件，可添加以下代码：
# tree <- read_tree("phylogenetic_tree.nwk")
# physeq_tree <- merge_phyloseq(physeq, tree)
# plot_tree(physeq_tree, color = "Group", label.tips = "Genus")

# 7. 生成分析报告摘要
# --------------------
report <- data.frame(
  Analysis = c("Sample Count", "OTU Count", "Alpha Diversity (Shannon)", "Beta Diversity (PERMANOVA R2)"),
  Result = c(
    nrow(sample_data),
    ntaxa(physeq),
    paste(round(mean(alpha_div$Shannon), 2), "±", round(sd(alpha_div$Shannon), 2)),
    paste(round(adonis_res$R2[1], 3), "p =", round(adonis_res$`Pr(>F)`[1], 3))
  )
)

write.xlsx(report, "16S_Results/Tables/Analysis_Summary.xlsx")

# 结束提示
cat("\n16S分析完成！结果已保存至 '16S_Results' 文件夹\n")
