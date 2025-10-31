# 加载必要的R包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("phyloseq", "vegan", "ggplot2", "dplyr", "tidyr", "microbiome", 
              "pheatmap", "DESeq2", "cowplot", "RColorBrewer")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

# 1. 数据读取和预处理 ---------------------------------------------------------
# 读取OTU表格
otu_file <- "OTUTaxonAnalysis_20250307_63S.full.csv"
otu_data <- read.csv(otu_file, stringsAsFactors = FALSE, check.names = FALSE)

# 提取样本丰度矩阵
sample_cols <- grep("^S\\d+", colnames(otu_data), value = TRUE)
otu_matrix <- as.matrix(otu_data[, sample_cols])
rownames(otu_matrix) <- otu_data$otu

# 提取分类学表格
tax_cols <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
tax_table <- as.matrix(otu_data[, tax_cols])
rownames(tax_table) <- otu_data$otu

# 创建样本分组信息（根据描述修改）
sampleinfo <- list(
  group1 = c("S009", "S008"),
  group2 = c("S0067"),
  group3 = c()  # 如果有更多样本组，在此添加
)

# 创建样本元数据
sample_data <- data.frame(
  Sample = sample_cols,
  stringsAsFactors = FALSE
)

# 添加分组信息
sample_data$Group <- NA
for (group_name in names(sampleinfo)) {
  sample_data$Group[sample_data$Sample %in% sampleinfo[[group_name]]] <- group_name
}

# 移除没有分组信息的样本
sample_data <- sample_data[!is.na(sample_data$Group), ]
otu_matrix <- otu_matrix[, sample_data$Sample]

# 创建phyloseq对象
OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(tax_table)
SAM <- sample_data(sample_data)

# 过滤低丰度OTU（至少在2个样本中出现且总丰度>5）
keep_otu <- rowSums(OTU > 0) >= 2 & rowSums(OTU) > 5
physeq <- phyloseq(OTU[keep_otu, ], TAX, SAM)

# 2. Alpha多样性分析 ---------------------------------------------------------
# 计算多样性指数
alpha_div <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$Sample <- rownames(alpha_div)
alpha_div <- merge(alpha_div, as(sample_data, "data.frame"), by = "Sample")

# 保存Alpha多样性结果
write.csv(alpha_div, "alpha_diversity_results.csv", row.names = FALSE)

# 绘制Alpha多样性箱线图
p_shannon <- ggplot(alpha_div, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2) +
  theme_minimal() +
  ggtitle("Shannon Diversity Index") +
  ylab("Shannon Index")

p_observed <- ggplot(alpha_div, aes(x = Group, y = Observed, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2) +
  theme_minimal() +
  ggtitle("Observed OTUs") +
  ylab("Number of OTUs")

# 保存Alpha多样性图
ggsave("alpha_diversity_shannon.png", plot = p_shannon, width = 8, height = 6, dpi = 300)
ggsave("alpha_diversity_observed.png", plot = p_observed, width = 8, height = 6, dpi = 300)

# 3. Beta多样性分析 ---------------------------------------------------------
# 转换为相对丰度
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# 计算Bray-Curtis距离
bray_dist <- phyloseq::distance(physeq_rel, method = "bray")

# PCoA分析
pcoa_result <- ordinate(physeq_rel, method = "PCoA", distance = bray_dist)

# 提取PCoA坐标
pcoa_df <- data.frame(pcoa_result$vectors)
pcoa_df$Sample <- rownames(pcoa_df)
pcoa_df <- merge(pcoa_df, as(sample_data, "data.frame"), by = "Sample")

# 计算方差解释比例
var_exp <- round(pcoa_result$values$Relative_eig * 100, 1)

# 绘制PCoA图
p_pcoa <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCoA Plot (Bray-Curtis)",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)")
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# 保存PCoA图
ggsave("beta_diversity_pcoa.png", plot = p_pcoa, width = 8, height = 6, dpi = 300)

# PERMANOVA检验
permanova_result <- adonis2(bray_dist ~ Group, data = as(sample_data, "data.frame"))
write.csv(as.data.frame(permanova_result), "permanova_results.csv")

# 4. 物种组成分析 -----------------------------------------------------------
# 在门水平上聚合
physeq_phylum <- tax_glom(physeq_rel, "phylum")

# 转换为数据框
phylum_df <- psmelt(physeq_phylum) %>%
  group_by(Sample, phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 计算平均丰度
phylum_avg <- phylum_df %>%
  group_by(Group, phylum) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(RelAbund = Abundance / sum(Abundance))

# 绘制堆叠柱状图
p_phylum <- ggplot(phylum_avg, aes(x = Group, y = RelAbund, fill = phylum)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Phylum-level Composition",
    x = "Sample Group",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存物种组成图
ggsave("phylum_composition.png", plot = p_phylum, width = 10, height = 6, dpi = 300)

# 5. 差异物种分析 -----------------------------------------------------------
# 在属水平上分析差异
physeq_genus <- tax_glom(physeq, "genus")

# 转换为DESeq2对象
diagdds <- phyloseq_to_deseq2(physeq_genus, ~ Group)

# 设置对照组（根据实际情况修改）
diagdds$Group <- relevel(diagdds$Group, ref = "group1")

# 运行DESeq2差异分析
diff_genus <- DESeq(diagdds, test = "Wald", fitType = "parametric")

# 提取结果（比较group2 vs group1）
res <- results(diff_genus, contrast = c("Group", "group2", "group1"))
res <- res[order(res$padj), ]

# 添加分类学信息
res_df <- as.data.frame(res)
res_df$genus <- rownames(res_df)
res_df <- merge(res_df, as.data.frame(TAX)[, c("genus", "family", "order", "class", "phylum")], 
                by = "genus", all.x = TRUE)

# 保存差异分析结果
write.csv(res_df, "differential_genus_group2_vs_group1.csv", row.names = FALSE)

# 绘制差异物种火山图
volcano_data <- res_df %>%
  mutate(
    sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, 
                 ifelse(log2FoldChange > 1, "Up", "Down"), "NS")
  )

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Differential Genus (group2 vs group1)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()

# 保存火山图
ggsave("differential_genus_volcano.png", plot = p_volcano, width = 10, height = 8, dpi = 300)

# 6. 生成综合报告 -----------------------------------------------------------
# 创建多面板图
combined_plot <- plot_grid(
  p_shannon + theme(legend.position = "none"),
  p_observed + theme(legend.position = "none"),
  p_pcoa + theme(legend.position = "none"),
  p_phylum + theme(legend.position = "none"),
  labels = c("A", "B", "C", "D"),
  ncol = 2
)

# 保存综合图
ggsave("comprehensive_analysis.png", plot = combined_plot, width = 12, height = 10, dpi = 300)

# 保存会话信息
sessionInfo()
