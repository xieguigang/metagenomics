# ==============================================================================
# 16S Microbiome Comprehensive Analysis Pipeline (Revised Version)
# Date: 2025-11-03
# Input: otu_table.txt (SampleID, Group, OTUs...)
# Optional Input: clinical_data.xlsx, PICRUSt2 outputs, Tax4Fun2 outputs
# Output: A series of PDF/PNG figures and XLSX result tables.
# ==============================================================================

# --- 0. Environment Setup, Library Installation, and Data Loading ---

# 0.1 Clean environment
rm(list = ls())

# 0.2 Install and load required packages
packages <- c(
  "tidyverse", "openxlsx", "phyloseq", "microbiome", "vegan", "picante",
  "DESeq2", "randomForest", "pROC", "VennDiagram", "igraph", "ggraph",
  "Hmisc", "corrplot", "pheatmap", "ape", "ggtree", "ggpubr", "RColorBrewer",
  "reshape2", "microeco", "CoDaSeq", "factoextra", "ALDEx2"
)

# Check and install CRAN packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  library(pkg, character.only = TRUE)
}

# Bioconductor packages
bioc_packages <- c("phyloseq", "microbiome", "DESeq2", "ALDEx2")
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
              "07_Function_Prediction", "08_Phenotype_Prediction")
sapply(sub_dirs, function(dir) dir.create(file.path(main_dir, dir), showWarnings = FALSE))

# 0.4 Load and prepare data with improved error handling
tryCatch({
  # Load OTU table with proper formatting
  otu_data <- read.csv("GROUP2/species.csv", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  
  # Extract metadata
  if ("class" %in% colnames(otu_data)) {
    meta_data <- data.frame(Group = otu_data$class, row.names = rownames(otu_data))
    otu_data$class <- NULL
  } else {
    stop("Column 'class' not found in the input data")
  }
  
  for(name in colnames(otu_data)) {
    otu_data[,name] = as.integer(otu_data[,name]*10000)
  }
  
  # Create OTU matrix
  otu_matrix <- as.matrix(otu_data)
  
  # Improved taxonomy parsing function
  parse_taxonomy <- function(tax_str) {
    # Handle multiple taxonomy formats
    if (grepl(";", tax_str)) {
      tax_vec <- unlist(strsplit(tax_str, ";"))
    } else if (grepl("\\|", tax_str)) {
      tax_vec <- unlist(strsplit(tax_str, "\\|"))
    } else {
      tax_vec <- tax_str
    }
    
    tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    tax_df <- data.frame(matrix(NA, nrow = 1, ncol = length(tax_levels)))
    colnames(tax_df) <- tax_levels
    
    for (i in 1:min(length(tax_vec), length(tax_levels))) {
      if (!is.na(tax_vec[i]) & nchar(tax_vec[i]) > 0) {
        # Remove prefix like k__, p__, etc.
        cleaned_tax <- gsub("^[kpcofgs]__", "", tax_vec[i])
        if (nchar(cleaned_tax) > 0) {
          tax_df[1, i] <- cleaned_tax
        }
      }
    }
    return(tax_df)
  }
  
  # Create taxonomy table
  tax_list <- lapply(colnames(otu_matrix), function(x) parse_taxonomy(x))
  tax_table <- do.call(rbind, tax_list)
  rownames(tax_table) <- colnames(otu_matrix)
  tax_table <- as.matrix(tax_table)
  
  # Create phyloseq object
  PS <- phyloseq(otu_table(otu_matrix, taxa_are_rows = FALSE),
                 sample_data(meta_data),
                 tax_table(tax_table))
  
  # Filter out low-prevalence taxa (present in less than 20% of samples)
  PS <- filter_taxa(PS, function(x) sum(x > 0) > (0.2 * length(x)), TRUE)
  
  cat("Data loaded successfully:\n")
  cat("Samples:", nsamples(PS), "\n")
  cat("Taxa:", ntaxa(PS), "\n")
  cat("Groups:", paste(unique(meta_data$Group), collapse = ", "), "\n")
  
}, error = function(e) {
  stop("Error loading data: ", e$message)
})

# --- 1. Community Structure Analysis (Enhanced) ---

# 1.1 Enhanced Community Composition Analysis
plot_composition <- function(ps, level, top_n = 15) {
  if (!level %in% rank_names(ps)) {
    message("Taxonomic level ", level, " not available. Skipping.")
    return(NULL)
  }
  
  # Aggregate at specified taxonomic level
  ps_glom <- tax_glom(ps, taxrank = level)
  ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
  
  # Melt for plotting
  df <- psmelt(ps_rel)
  df$Abundance <- df$Abundance * 100
  
  # Identify top taxa
  top_taxa <- df %>%
    group_by(across(all_of(level))) %>%
    dplyr::summarize(mean_abund = mean(Abundance, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(mean_abund)) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= top_n) %>%
    pull(!!sym(level))
  
  # Group low abundance taxa as "Other"
  df[[level]] <- ifelse(df[[level]] %in% top_taxa, as.character(df[[level]]), "Other")
  
  # Create color palette
  n_colors <- length(unique(df[[level]]))
  color_palette <- colorRampPalette(RColorBrewer::brewer.pal(min(n_colors, 8), "Set3"))(n_colors)
  
  # Plot
  p <- ggplot(df, aes(x = Sample, y = Abundance, fill = get(level))) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~Group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_palette) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          panel.spacing = unit(0.1, "lines")) +
    labs(x = "Samples", y = "Relative Abundance (%)", fill = level,
         title = paste("Microbial Community Composition at", level, "Level"))
  
  # Save plots
  ggsave(filename = file.path(main_dir, "01_Community_Structure", 
                              paste0("Composition_", level, ".pdf")), 
         plot = p, width = 12, height = 8)
  ggsave(filename = file.path(main_dir, "01_Community_Structure", 
                              paste0("Composition_", level, ".png")), 
         plot = p, width = 12, height = 8, dpi = 300)
  
  # Save data
  write.xlsx(df, file = file.path(main_dir, "01_Community_Structure", 
                                  paste0("Composition_", level, "_data.xlsx")))
  
  return(p)
}

# Generate composition plots for different taxonomic levels
for (lvl in c("Phylum", "Class", "Order", "Family", "Genus")) {
  plot_composition(PS, lvl)
}

# 1.2 Enhanced Venn Diagram Analysis
generate_venn_diagram <- function(ps, min_abundance = 0) {
  # Extract group information
  group_list <- split(sample_names(ps), sample_data(ps)$Group)
  
  # Get OTUs present in each group
  otu_lists <- lapply(group_list, function(group_samples) {
    otu_subset <- otu_table(ps)[group_samples, ]
    present_otus <- colnames(otu_subset)[colSums(otu_subset) > min_abundance]
    return(present_otus)
  })
  
  # Create Venn diagram
  venn.plot <- venn.diagram(
    x = otu_lists,
    category.names = names(otu_lists),
    filename = NULL,
    output = TRUE,
    col = "transparent",
    fill = RColorBrewer::brewer.pal(length(otu_lists), "Set3"),
    alpha = 0.50,
    cex = 1.2,
    cat.cex = 1.1,
    cat.col = "black",
    margin = 0.1
  )
  
  # Save plots
  pdf(file.path(main_dir, "01_Community_Structure", "Venn_Diagram.pdf"), 
      width = 8, height = 8)
  grid.draw(venn.plot)
  dev.off()
  
  png(file.path(main_dir, "01_Community_Structure", "Venn_Diagram.png"), 
      width = 2400, height = 2400, res = 300)
  grid.draw(venn.plot)
  dev.off()
  
  # Calculate intersection statistics
  intersection_stats <- list()
  all_combinations <- combn(names(otu_lists), 2, simplify = FALSE)
  
  for (comb in all_combinations) {
    shared_otus <- intersect(otu_lists[[comb[1]]], otu_lists[[comb[2]]])
    intersection_stats[[paste(comb[1], "vs", comb[2])]] <- data.frame(
      Group1 = comb[1],
      Group2 = comb[2],
      Shared_OTUs = length(shared_otus),
      Unique_Group1 = length(setdiff(otu_lists[[comb[1]]], otu_lists[[comb[2]]])),
      Unique_Group2 = length(setdiff(otu_lists[[comb[2]]], otu_lists[[comb[1]]]))
    )
  }
  
  intersection_df <- do.call(rbind, intersection_stats)
  
  max_size = max(sapply(otu_lists, function(l)length(l)));
  
  for(name in names(otu_lists)) {
    v = otu_lists[[name]];
    
    if (length(v) < max_size) {
      v = append(v, rep("", times = max_size-length(v)));
      otu_lists[[name]] = v;
    }
  }
  
  print("view of the data dumps:");
  str(otu_lists)
  print(intersection_df)
  
  write.xlsx(list(OTU_Lists = otu_lists, Intersection_Stats = intersection_df),
             file = file.path(main_dir, "01_Community_Structure", "Venn_Analysis.xlsx"))
  
  return(otu_lists)
}

venn_results <- generate_venn_diagram(PS)

# 1.3 Enhanced Heatmap Analysis
generate_heatmap <- function(ps, level = "Genus", top_n = 30) {
  # Aggregate at specified taxonomic level
  ps_glom <- tax_glom(ps, taxrank = level)
  ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
  
  # Get top abundant taxa
  top_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:min(top_n, ntaxa(ps_rel))]
  ps_top <- prune_taxa(top_taxa, ps_rel)
  
  # Prepare data for heatmap
  heatmap_data <- as.matrix(otu_table(ps_top))
  rownames(heatmap_data) <- paste0("Sample_", 1:nrow(heatmap_data))
  
  # Create annotation data
  annotation_df <- data.frame(
    Group = sample_data(ps_top)$Group,
    row.names = rownames(heatmap_data)
  )
  
  # Create heatmap
  pheatmap(t(heatmap_data),
           annotation_col = annotation_df,
           show_colnames = FALSE,
           show_rownames = TRUE,
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           fontsize_row = 8,
           main = paste("Heatmap of Top", top_n, level, "Abundance"),
           filename = file.path(main_dir, "01_Community_Structure", 
                                paste0("Heatmap_", level, ".pdf")),
           width = 12, height = 10)
  
  # PNG version
  pheatmap(t(heatmap_data),
           annotation_col = annotation_df,
           show_colnames = FALSE,
           show_rownames = TRUE,
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           fontsize_row = 8,
           main = paste("Heatmap of Top", top_n, level, "Abundance"),
           filename = file.path(main_dir, "01_Community_Structure", 
                                paste0("Heatmap_", level, ".png")),
           width = 12, height = 10)
  
  # Save data
  write.xlsx(heatmap_data, 
             file = file.path(main_dir, "01_Community_Structure", 
                              paste0("Heatmap_", level, "_data.xlsx")))
}

generate_heatmap(PS, "Genus", 30)

# --- 2. Enhanced Diversity Analysis --- [3,4](@ref)

# 2.1 Comprehensive Alpha Diversity Analysis
analyze_alpha_diversity <- function(ps) {
  print(ps);
  
  # Calculate alpha diversity indices
  alpha_metrics <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
  alpha_div <- estimate_richness(ps, measures = alpha_metrics)
  
  # Add Goods Coverage
  goods_coverage <- function(x) { 1 - sum(x == 1) / sum(x) }
  coverage <- apply(otu_table(ps), 1, goods_coverage)
  alpha_div$GoodsCoverage <- coverage
  
  # Add sample information
  alpha_div$Group <- sample_data(ps)$Group
  alpha_div$SampleID <- rownames(alpha_div)
  
  # Save alpha diversity table
  write.xlsx(alpha_div, file = file.path(main_dir, "02_Diversity", "Alpha_Diversity_Table.xlsx"))
  
  # Statistical testing
  alpha_stats <- list()
  for (metric in alpha_metrics) {
    message(metric);
    # Kruskal-Wallis test
    v = alpha_div[[metric]]
    if (all(is.nan(v))) {
      alpha_stats[[metric]] <- data.frame(
        Metric = metric,
        Kruskal_Wallis_Statistic = NA,
        Kruskal_Wallis_p_value = 1
      )    
    } else {
      kruskal_test <- kruskal.test(alpha_div[[metric]] ~ alpha_div$Group)
      alpha_stats[[metric]] <- data.frame(
        Metric = metric,
        Kruskal_Wallis_Statistic = kruskal_test$statistic,
        Kruskal_Wallis_p_value = kruskal_test$p.value
      )      
    }
  }
  
  alpha_stats_df <- do.call(rbind, alpha_stats)
  write.xlsx(alpha_stats_df, file = file.path(main_dir, "02_Diversity", "Alpha_Diversity_Statistics.xlsx"))
  
  # Plot alpha diversity
  plot_alpha_diversity <- function(alpha_div, metric) {
    p <- ggplot(alpha_div, aes(x = Group, y = get(metric), fill = Group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
      stat_compare_means(method = "kruskal.test", label = "p.format") +
      theme_bw() +
      labs(x = "Group", y = metric, 
           title = paste("Alpha Diversity:", metric),
           fill = "Group") +
      theme(legend.position = "right")
    
    ggsave(file.path(main_dir, "02_Diversity", paste0("Alpha_", metric, ".pdf")), 
           plot = p, width = 8, height = 6)
    ggsave(file.path(main_dir, "02_Diversity", paste0("Alpha_", metric, ".png")), 
           plot = p, width = 8, height = 6, dpi = 300)
    
    return(p)
  }
  
  # Generate plots for all metrics
  for (metric in alpha_metrics) {
    plot_alpha_diversity(alpha_div, metric)
  }
  
  return(alpha_div)
}

alpha_results <- analyze_alpha_diversity(PS)

# 2.2 Enhanced Beta Diversity Analysis [4](@ref)
analyze_beta_diversity <- function(ps) {
  # Calculate distance matrices
  dist_methods <- c("bray", "jaccard", "unifrac", "wunifrac")
  dist_matrices <- list()
  
  for (method in dist_methods) {
    tryCatch({
      dist_matrices[[method]] <- phyloseq::distance(ps, method = method)
    }, error = function(e) {
      message("Could not calculate distance using method: ", method, ". Error: ", e$message)
    })
  }
  
  # Ordination analysis
  ordination_results <- list()
  for (method in names(dist_matrices)) {
    # PCoA
    tryCatch({
      pcoa_result <- ordinate(ps, method = "PCoA", distance = dist_matrices[[method]])
      ordination_results[[paste("PCoA", method, sep = "_")]] <- pcoa_result
      
      p <- plot_ordination(ps, pcoa_result, color = "Group", shape = "Group") +
        geom_point(size = 3, alpha = 0.8) +
        stat_ellipse(level = 0.8) +
        theme_bw() +
        labs(title = paste("PCoA Plot (", method, " distance)", sep = ""))
      
      ggsave(file.path(main_dir, "02_Diversity", paste0("PCoA_", method, ".pdf")), 
             plot = p, width = 8, height = 6)
      ggsave(file.path(main_dir, "02_Diversity", paste0("PCoA_", method, ".png")), 
             plot = p, width = 8, height = 6, dpi = 300)
    }, error = function(e) {
      message("Error in PCoA with method ", method, ": ", e$message)
    })
    
    # NMDS
    tryCatch({
      nmds_result <- ordinate(ps, method = "NMDS", distance = dist_matrices[[method]])
      ordination_results[[paste("NMDS", method, sep = "_")]] <- nmds_result
      
      p <- plot_ordination(ps, nmds_result, color = "Group", shape = "Group") +
        geom_point(size = 3, alpha = 0.8) +
        stat_ellipse(level = 0.8) +
        theme_bw() +
        labs(title = paste("NMDS Plot (", method, " distance)", sep = ""))
      
      ggsave(file.path(main_dir, "02_Diversity", paste0("NMDS_", method, ".pdf")), 
             plot = p, width = 8, height = 6)
      ggsave(file.path(main_dir, "02_Diversity", paste0("NMDS_", method, ".png")), 
             plot = p, width = 8, height = 6, dpi = 300)
    }, error = function(e) {
      message("Error in NMDS with method ", method, ": ", e$message)
    })
  }
  
  # Statistical tests
  beta_stats <- list()
  group_factor <- sample_data(ps)$Group
  
  for (method in names(dist_matrices)) {
    # ANOSIM
    anosim_res <- anosim(dist_matrices[[method]], group_factor)
    # PERMANOVA (adonis)
    adonis_res <- adonis2(dist_matrices[[method]] ~ group_factor)
    
    beta_stats[[method]] <- data.frame(
      Method = method,
      ANOSIM_R = anosim_res$statistic,
      ANOSIM_p = anosim_res$signif,
      PERMANOVA_R2 = adonis_res$R2[1],
      PERMANOVA_F = adonis_res$F[1],
      PERMANOVA_p = adonis_res$`Pr(>F)`[1]
    )
  }
  
  beta_stats_df <- do.call(rbind, beta_stats)
  write.xlsx(beta_stats_df, file = file.path(main_dir, "02_Diversity", "Beta_Diversity_Statistics.xlsx"))
  
  # UPGMA clustering
  for (method in names(dist_matrices)) {
    tryCatch({
      hc <- hclust(dist_matrices[[method]], method = "average")
      
      pdf(file.path(main_dir, "02_Diversity", paste0("UPGMA_", method, ".pdf")), 
          width = 10, height = 6)
      plot(hc, main = paste("UPGMA Clustering (", method, " distance)", sep = ""),
           xlab = "", sub = "")
      dev.off()
      
      png(file.path(main_dir, "02_Diversity", paste0("UPGMA_", method, ".png")), 
          width = 2400, height = 1800, res = 300)
      plot(hc, main = paste("UPGMA Clustering (", method, " distance)", sep = ""),
           xlab = "", sub = "")
      dev.off()
    }, error = function(e) {
      message("Error in UPGMA with method ", method, ": ", e$message)
    })
  }
  
  return(list(distances = dist_matrices, stats = beta_stats_df))
}

beta_results <- analyze_beta_diversity(PS)

# --- 3. Enhanced Differential Abundance Analysis --- [2](@ref)

# 3.1 ALDEx2 for Compositional Data Analysis
run_aldex2_analysis <- function(ps, mc.samples = 128) {
  # 提取数据
  otu_table_aldex <- t(otu_table(ps))
  conditions <- as.character(sample_data(ps)$Group)
  
  # 获取所有唯一组别
  unique_groups <- unique(conditions)
  num_groups <- length(unique_groups)
  
  cat("检测到", num_groups, "个组别:", paste(unique_groups, collapse = ", "), "\n")
  
  if (num_groups < 2) {
    stop("需要至少2个组别进行比较")
  }
  
  # 运行ALDEx2 CLR转换
  aldex_obj <- aldex.clr(otu_table_aldex, conditions, 
                         mc.samples = mc.samples, 
                         verbose = FALSE)
  
  results_list <- list()
  
  if (num_groups == 2) {
    # 两组比较：使用t检验和效应量
    aldex_ttest <- aldex.ttest(aldex_obj, conditions)
    aldex_effect <- aldex.effect(aldex_obj, conditions, 
                                 include.sample.summary = FALSE)
    aldex_results <- data.frame(aldex_ttest, aldex_effect)
    
    # 标识显著特征
    sig_features <- which(aldex_results$we.eBH < 0.05 & abs(aldex_results$effect) > 1)
    cat("发现", length(sig_features), "个显著差异特征\n")
    
    results_list[["All_Groups"]] <- aldex_results
    
  } else {
    # 多组比较：使用Kruskal-Wallis检验
    cat("进行多组Kruskal-Wallis检验\n")
    aldex_kw <- aldex.kw(aldex_obj)
    results_list[["Kruskal_Wallis"]] <- aldex_kw
    
    # 标识显著特征（基于Kruskal-Wallis检验）
    sig_features <- which(aldex_kw$kw.eBH < 0.05)
    cat("Kruskal-Wallis检验发现", length(sig_features), "个显著差异特征\n")
    
    # 进行两两比较
    cat("进行两两组间比较\n")
    pairwise_combinations <- combn(unique_groups, 2, simplify = FALSE)
    
    for (pair in pairwise_combinations) {
      group1 <- pair[1]
      group2 <- pair[2]
      comparison_name <- paste(group1, "vs", group2, sep = "_")
      
      # 选择属于这两个组的样本
      group1_samples <- which(conditions == group1)
      group2_samples <- which(conditions == group2)
      selected_samples <- c(group1_samples, group2_samples)
      selected_conditions <- conditions[selected_samples]
      
      # 提取对应的OTU数据
      selected_otu <- otu_table_aldex[selected_samples, ]
      
      # 运行两组比较
      tryCatch({
        aldex_pairwise <- aldex.clr(selected_otu, selected_conditions, 
                                    mc.samples = mc.samples, 
                                    verbose = FALSE)
        aldex_ttest_pairwise <- aldex.ttest(aldex_pairwise, selected_conditions)
        aldex_effect_pairwise <- aldex.effect(aldex_pairwise, selected_conditions, 
                                              include.sample.summary = FALSE)
        aldex_results_pairwise <- data.frame(aldex_ttest_pairwise, aldex_effect_pairwise)
        
        results_list[[comparison_name]] <- aldex_results_pairwise
        
        sig_pairwise <- which(aldex_results_pairwise$we.eBH < 0.05 & 
                                abs(aldex_results_pairwise$effect) > 1)
        cat("比较", comparison_name, "发现", length(sig_pairwise), "个显著差异特征\n")
        
      }, error = function(e) {
        cat("比较", comparison_name, "时出错:", e$message, "\n")
      })
    }
  }
  
  return(results_list)
}

# 执行ALDEx2分析
tryCatch({
  aldex2_results <- run_aldex2_analysis(PS, mc.samples = 128)
  
  # 保存结果
  if (length(aldex2_results) > 0) {
    write.xlsx(aldex2_results, 
               file = file.path(main_dir, "03_Differential_Analysis", "ALDEx2_Results.xlsx"))
    
    # 生成摘要报告
    summary_stats <- list()
    for (result_name in names(aldex2_results)) {
      result <- aldex2_results[[result_name]]
      
      if (grepl("vs", result_name)) {
        # 两两比较结果
        sig_count <- sum(result$we.eBH < 0.05 & abs(result$effect) > 1, na.rm = TRUE)
        summary_stats[[result_name]] <- data.frame(
          Comparison = result_name,
          Significant_Features = sig_count,
          Total_Features = nrow(result)
        )
      } else if (result_name == "Kruskal_Wallis") {
        # Kruskal-Wallis结果
        sig_count <- sum(result$kw.eBH < 0.05, na.rm = TRUE)
        summary_stats[[result_name]] <- data.frame(
          Comparison = "Kruskal_Wallis (Overall)",
          Significant_Features = sig_count,
          Total_Features = nrow(result)
        )
      }
    }
    
    if (length(summary_stats) > 0) {
      summary_df <- do.call(rbind, summary_stats)
      write.xlsx(summary_df, 
                 file = file.path(main_dir, "03_Differential_Analysis", "ALDEx2_Summary.xlsx"))
      
      # 绘制摘要图
      p_summary <- ggplot(summary_df, aes(x = Comparison, y = Significant_Features, fill = Comparison)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "ALDEx2差异分析摘要", 
             x = "比较组", y = "显著差异特征数量")
      
      ggsave(file.path(main_dir, "03_Differential_Analysis", "ALDEx2_Summary_Plot.pdf"), 
             plot = p_summary, width = 10, height = 6)
      ggsave(file.path(main_dir, "03_Differential_Analysis", "ALDEx2_Summary_Plot.png"), 
             plot = p_summary, width = 10, height = 6, dpi = 300)
    }
    
    cat("ALDEx2分析完成，结果已保存\n")
  }
  
}, error = function(e) {
  cat("ALDEx2分析失败:", e$message, "\n")
})

# 3.2 Enhanced Random Forest Analysis
enhanced_random_forest <- function(ps, ntree = 1000) {
  # Prepare data
  otu_mat <- t(otu_table(ps))
  group_vec <- as.factor(sample_data(ps)$Group)
  
  # Remove near-zero variance features
  nzv <- caret::nearZeroVar(otu_mat)
  if (length(nzv) > 0) {
    otu_mat <- otu_mat[, -nzv]
  }
  
  # Train random forest model
  set.seed(123)
  rf_model <- randomForest(x = otu_mat, y = group_vec, 
                           importance = TRUE, 
                           ntree = ntree,
                           proximity = TRUE)
  
  # Extract importance
  rf_importance <- importance(rf_model)
  rf_importance_df <- data.frame(
    Taxon = rownames(rf_importance),
    MeanDecreaseAccuracy = rf_importance[, "MeanDecreaseAccuracy"],
    MeanDecreaseGini = rf_importance[, "MeanDecreaseGini"]
  )
  rf_importance_df <- rf_importance_df[order(rf_importance_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
  
  # Add taxonomy information
  tax_info <- as.data.frame(tax_table(ps))
  tax_info$Taxon <- rownames(tax_info)
  rf_importance_df <- merge(rf_importance_df, tax_info, by = "Taxon", all.x = TRUE)
  
  # Save results
  write.xlsx(rf_importance_df, 
             file = file.path(main_dir, "03_Differential_Analysis", "RandomForest_Importance.xlsx"))
  
  # Plot top important features
  top_n <- min(20, nrow(rf_importance_df))
  p_rf <- ggplot(rf_importance_df[1:top_n, ], 
                 aes(x = reorder(Taxon, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_bw() +
    labs(x = "Taxon", y = "Mean Decrease in Accuracy", 
         title = paste("Random Forest Feature Importance (Top", top_n, ")"))
  
  ggsave(file.path(main_dir, "03_Differential_Analysis", "RandomForest_Importance.pdf"), 
         plot = p_rf, width = 10, height = 8)
  ggsave(file.path(main_dir, "03_Differential_Analysis", "RandomForest_Importance.png"), 
         plot = p_rf, width = 10, height = 8, dpi = 300)
  
  # Plot proximity matrix
  proximity_df <- as.data.frame(cmdscale(rf_model$proximity, k = 2))
  proximity_df$Group <- group_vec
  colnames(proximity_df)[1:2] <- c("MDS1", "MDS2")
  
  p_proximity <- ggplot(proximity_df, aes(x = MDS1, y = MDS2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.8) +
    theme_bw() +
    labs(title = "Random Forest Proximity Plot")
  
  ggsave(file.path(main_dir, "03_Differential_Analysis", "RandomForest_Proximity.pdf"), 
         plot = p_proximity, width = 8, height = 6)
  ggsave(file.path(main_dir, "03_Differential_Analysis", "RandomForest_Proximity.png"), 
         plot = p_proximity, width = 8, height = 6, dpi = 300)
  
  return(list(model = rf_model, importance = rf_importance_df))
}

rf_results <- enhanced_random_forest(PS)

# --- 4. Enhanced Network Analysis using microeco --- [1](@ref)

analyze_network <- function(ps, cor_method = "spearman", cor_cutoff = 0.6, p_cutoff = 0.05) {
  # Convert to microeco format
  tryCatch({
    # Extract components for microeco
    otu_table_meco <- as.data.frame(t(otu_table(ps)))
    sample_info_meco <- as.data.frame(sample_data(ps))
    taxonomy_table_meco <- as.data.frame(tax_table(ps))
    
    # Create microeco object
    dataset <- microtable$new(sample_table = sample_info_meco, 
                              otu_table = otu_table_meco, 
                              tax_table = taxonomy_table_meco)
    
    # Calculate correlation network
    t1 <- trans_network$new(
      dataset = dataset,
      cal_cor = "WGCNA",  # Use WGCNA for correlation calculation
      taxa_level = "OTU",
      filter_thres = 0.001,
      cor_method = cor_method
    )
    
    # Build network
    t1$cal_network(p_thres = p_cutoff, COR_cut = cor_cutoff)
    
    # Calculate modules
    t1$cal_module(method = "cluster_fast_greedy", module_name_prefix = "M")
    
    # Save network for Gephi
    t1$save_network(file.path(main_dir, "05_Network_Analysis", "network.gexf"))
    
    # Basic network statistics
    network_stats <- data.frame(
      Nodes = igraph::vcount(t1$res_network),
      Edges = igraph::ecount(t1$res_network),
      Density = igraph::edge_density(t1$res_network),
      Transitivity = igraph::transitivity(t1$res_network),
      Average_Degree = mean(igraph::degree(t1$res_network))
    )
    
    write.xlsx(network_stats, 
               file = file.path(main_dir, "05_Network_Analysis", "Network_Statistics.xlsx"))
    
    # Plot network
    pdf(file.path(main_dir, "05_Network_Analysis", "Network_Plot.pdf"), 
        width = 10, height = 8)
    plot(t1$res_network, 
         vertex.size = 5,
         vertex.label = NA,
         edge.width = 1,
         main = "Microbial Co-occurrence Network")
    dev.off()
    
    png(file.path(main_dir, "05_Network_Analysis", "Network_Plot.png"), 
        width = 2000, height = 1600, res = 300)
    plot(t1$res_network, 
         vertex.size = 5,
         vertex.label = NA,
         edge.width = 1,
         main = "Microbial Co-occurrence Network")
    dev.off()
    
    return(t1)
    
  }, error = function(e) {
    message("Network analysis failed: ", e$message)
    return(NULL)
  })
}

network_result <- analyze_network(PS)

# --- 5. Enhanced Functional Prediction with Tax4Fun2 --- [5](@ref)

analyze_functional_prediction <- function(ps, tax4fun2_ref_path = NULL) {
  # Check if Tax4Fun2 reference data is available
  if (is.null(tax4fun2_ref_path)) {
    message("Tax4Fun2 reference data path not provided. Skipping functional prediction.")
    return(NULL)
  }
  
  tryCatch({
    # Prepare data for Tax4Fun2
    otu_table_for_tax4fun <- as.data.frame(t(otu_table(ps)))
    
    # Run Tax4Fun2 prediction
    func_prediction <- Tax4Fun2::tax4fun2(
      otu_table = otu_table_for_tax4fun,
      reference_data = tax4fun2_ref_path,
      pathway_normalization = TRUE
    )
    
    # Extract and save results
    kegg_pathways <- func_prediction$pathway_abundance
    kegg_modules <- func_prediction$module_abundance
    
    write.xlsx(kegg_pathways, 
               file = file.path(main_dir, "07_Function_Prediction", "KEGG_Pathways.xlsx"))
    write.xlsx(kegg_modules, 
               file = file.path(main_dir, "07_Function_Prediction", "KEGG_Modules.xlsx"))
    
    # Create functional profiles
    func_profile <- as.data.frame(t(kegg_pathways))
    func_profile$Group <- sample_data(ps)$Group
    
    # Plot top functional pathways
    top_pathways <- names(sort(colMeans(kegg_pathways), decreasing = TRUE)[1:20])
    func_melt <- reshape2::melt(func_profile[, c(top_pathways, "Group")], id.vars = "Group")
    
    p_pathway <- ggplot(func_melt, aes(x = Group, y = value, fill = Group)) +
      geom_boxplot() +
      facet_wrap(~ variable, scales = "free_y", ncol = 5) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Top 20 KEGG Pathways by Abundance")
    
    ggsave(file.path(main_dir, "07_Function_Prediction", "KEGG_Pathways_Plot.pdf"), 
           plot = p_pathway, width = 15, height = 10)
    ggsave(file.path(main_dir, "07_Function_Prediction", "KEGG_Pathways_Plot.png"), 
           plot = p_pathway, width = 15, height = 10, dpi = 300)
    
    return(func_prediction)
    
  }, error = function(e) {
    message("Functional prediction with Tax4Fun2 failed: ", e$message)
    return(NULL)
  })
}

# Uncomment and modify the path if you have Tax4Fun2 reference data
# func_results <- analyze_functional_prediction(PS, tax4fun2_ref_path = "path/to/Tax4Fun2_ReferenceData")

# --- 6. Save Session Information ---

# Save R session information
sink(file.path(main_dir, "Session_Info.txt"))
print(sessionInfo())
sink()

# Save analysis completion message
message("Microbiome analysis completed successfully!")
message("Results are saved in: ", main_dir)
message("Analysis completed at: ", Sys.time())

# --- End of Enhanced Script ---