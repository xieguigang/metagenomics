# 加载必要的R包
required_packages <- c("vegan", "ape", "ggplot2", "dplyr", "tidyr", "reshape2",
                       "pairwiseAdonis", "microbiome", "microbiomeMarker",
                       "ANCOMBC", "ALDEx2", "igraph", "ggpubr", "GUniFrac", "VennDiagram",
                       "proxy", "ggraph", "tidygraph", "ggrepel", "corrplot", "RColorBrewer")

# 检查并安装缺失的包
for(pkg in required_packages){
    if(!require(pkg, character.only = TRUE)){
        if(pkg %in% c("ANCOMBC", "microbiomeMarker", "ALDEx2","pairwiseAdonis")){
            if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            BiocManager::install(pkg)
        } else {
            install.packages(pkg)
        }
    }
    message(pkg);
    library(pkg, character.only = TRUE)
}

process_analysis = function(data, sample_info, outputdir) {
    # 读取数据
    sample_info = read.csv(sample_info, row.names = 1, check.names = FALSE);
    otu_table <- read.csv(data, row.names = 1, check.names = FALSE);
    dir.create(outputdir, showWarnings = FALSE);

    otu_table[,"taxonomy"] = NULL;
    otu_table[,"Tree"] = NULL;
    otu_table[,"OTUs"] = NULL;
    otu_table = t(otu_table);

    otu_table = as.data.frame(otu_table[rownames(sample_info),]);

    print("view of the data for run analysis:");
    str(sample_info);
    str(otu_table);

    analyze_by_sample_type(otu_table, sample_info, outputdir);
}

# 针对特定样本类型的分析函数
analyze_by_sample_type <- function(otu_table, sample_info, output_dir) {
    message("准备数据");
    print(head(sample_info));
    otu_table = cbind(ID = rownames(sample_info), class = sample_info$sample_info, otu_table);
    sample_type = "metagenomics";
    rownames(otu_table) = rownames(sample_info);

    message("数据将要输出到这个文件夹：",output_dir);
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE);
    write.csv(otu_table, file = file.path(output_dir,"data.csv"), row.names = FALSE);

    # 数据预处理和归一化函数
    preprocess_data <- function(otu_table) {
        message("开始处理OTU数据表");

        # 提取样本信息
        sample_info <- otu_table[, 1:2]
        otu_data <- otu_table[, -c(1:2)]
        scale <- 100000000

        for(name in colnames(otu_data ) ) {            
            v = as.numeric( otu_data[, name]) * scale; 
            v = as.integer(v); 

            if (all(v == 0)) {
                v = as.integer( runif(length(v)) * 10)
            }

            otu_data[, name] = v;
        }

        print(head(otu_data ));

        # 转换OTU数据为数值矩阵
        otu_matrix <- as.matrix(otu_data)
        rownames(otu_matrix) <- sample_info$ID

        # 数据归一化步骤
        cat("进行重抽样标准化...\n")
        min_depth <- min(rowSums(otu_matrix))
        rarefied_otu <- rrarefy(otu_matrix, sample = min_depth)

        # 相对丰度标准化
        relative_abu <- sweep(otu_matrix, 1, rowSums(otu_matrix), "/")
        relative_abu[is.na(relative_abu)] <- 0

        # CSS标准化
        css_normalize <- function(otu_mat) {
            lib_sizes <- rowSums(otu_mat)
            css_factors <- lib_sizes / quantile(lib_sizes, 0.75)
            css_norm <- sweep(otu_mat, 1, css_factors, "/")
            return(css_norm)
        }

        css_otu <- css_normalize(otu_matrix)

        # 对数转换
        log_otu <- log10(otu_matrix + 1)

        return(list(
            sample_info = sample_info,
            otu_matrix = otu_matrix,
            rarefied_otu = rarefied_otu,
            relative_abu = relative_abu,
            css_otu = css_otu,
            log_otu = log_otu
        ))
    }

    # α多样性分析函数（BMI组间比较）
    alpha_diversity_analysis <- function(processed_data) {
        rarefied_otu <- processed_data$rarefied_otu
        sample_info <- processed_data$sample_info

        # 计算α多样性指数
        shannon <- vegan::diversity(rarefied_otu, index = "shannon")
        simpson <- vegan::diversity(rarefied_otu, index = "simpson")
        richness <- rowSums(rarefied_otu > 0)
        chao1 <- estimateR(rarefied_otu)[2, ]

        alpha_df <- data.frame(
            ID = rownames(rarefied_otu),
            Shannon = shannon,
            Simpson = simpson,
            Richness = richness,
            Chao1 = chao1
        )

        # 合并样本信息
        alpha_df <- merge(alpha_df, sample_info, by = "ID")

        print("alpha_diversity");
        print(head(alpha_df));

        # BMI组间比较
        results <- list()
        diversity_indices <- c("Shannon", "Simpson", "Richness", "Chao1")

        for(index in diversity_indices) {
            # 检查是否有多个BMI组别
            bmi_groups <- unique(alpha_df$class)
            if(length(bmi_groups) < 2) {
                results[[index]] <- list(
                    test_result = NULL,
                    test_type = "Insufficient groups",
                    message = "需要至少两个BMI组别进行比较"
                )
                next
            }

            # 正态性检验（对每个组别）
            normality_p <- c();

            for(group in bmi_groups) {
                group_data <- alpha_df[alpha_df$class == group, index];

                if(length(group_data) > 3) { # 需要足够样本进行正态性检验
                    print("view group data for shapiro.test:");
                    print(group_data);

                    if (all(group_data == group_data[0])) {
                        # all 'x' values are identical
                        normality_p <- c(normality_p, 0)
                    } else {
                        normality_test <- shapiro.test(group_data)
                        normality_p <- c(normality_p, normality_test$p.value)
                    }
                }
            }

            # 如果所有组别数据正态且组别数=2，使用t检验，否则使用ANOVA或Kruskal-Wallis
            if(length(bmi_groups) == 2 && all(normality_p > 0.05)) {
                # 两组比较，数据正态，使用t检验
                group1_data <- alpha_df[alpha_df$class == bmi_groups[1], index]
                group2_data <- alpha_df[alpha_df$class == bmi_groups[2], index]
                test_result <- t.test(group1_data, group2_data);
                test_type <- "Student's t-test"
            } else if(length(bmi_groups) == 2) {
                # 两组比较，数据非正态，使用Mann-Whitney检验
                group1_data <- alpha_df[alpha_df$class == bmi_groups[1], index]
                group2_data <- alpha_df[alpha_df$class == bmi_groups[2], index]
                test_result <- wilcox.test(group1_data, group2_data)
                test_type <- "Mann-Whitney U test"
            } else if(length(bmi_groups) > 2 && all(normality_p > 0.05)) {
                # 多组比较，数据正态，使用ANOVA
                anova_model <- aov(as.formula(paste(index, "~ class")), data = alpha_df)
                test_result <- summary(anova_model)
                test_type <- "ANOVA"
            } else if(length(bmi_groups) > 2) {
                # 多组比较，数据非正态，使用Kruskal-Wallis检验
                test_result <- kruskal.test(as.formula(paste(index, "~ class")), data = alpha_df)
                test_type <- "Kruskal-Wallis test"
            }

            results[[index]] <- list(
                test_result = test_result,
                test_type = test_type,
                groups = bmi_groups,
                normality_p = normality_p
            )
        }

        return(list(alpha_df = alpha_df, test_results = results))
    }

    # β多样性分析函数（BMI组间比较）
    beta_diversity_analysis <- function(processed_data) {
        rarefied_otu <- processed_data$rarefied_otu
        sample_info <- processed_data$sample_info

        # 计算距离矩阵
        bray_dist <- vegdist(rarefied_otu, method = "bray")
        jaccard_dist <- vegdist(rarefied_otu > 0, method = "jaccard")

        # PCoA分析
        pcoa_bray <- cmdscale(bray_dist, k = 3, eig = TRUE)
        pcoa_jaccard <- cmdscale(jaccard_dist, k = 3, eig = TRUE)

        # PERMANOVA检验（BMI组间比较）
        if(length(unique(sample_info$class)) > 1) {
            permanova_bray <- adonis2(bray_dist ~ class,
                                      data = sample_info, permutations = 999)

            permanova_jaccard <- adonis2(jaccard_dist ~ class,
                                         data = sample_info, permutations = 999)
        } else {
            permanova_bray <- NULL
            permanova_jaccard <- NULL
        }

        return(list(
            bray_dist = bray_dist,
            jaccard_dist = jaccard_dist,
            pcoa_bray = pcoa_bray,
            pcoa_jaccard = pcoa_jaccard,
            permanova_bray = permanova_bray,
            permanova_jaccard = permanova_jaccard
        ))
    }

    # 群落相似性分析函数（BMI组内和组间相似性）
    community_similarity_analysis <- function(processed_data) {
        relative_abu <- processed_data$relative_abu
        sample_info <- processed_data$sample_info

        # 计算相似性矩阵
        bray_sim <- as.matrix(1 - vegdist(relative_abu, method = "bray"))

        # 计算组内和组间相似性
        bmi_groups <- unique(sample_info$class)
        similarity_results <- data.frame()

        for(i in 1:length(bmi_groups)) {
            for(j in i:length(bmi_groups)) {
                group1 <- bmi_groups[i]
                group2 <- bmi_groups[j]

                group1_samples <- sample_info$ID[sample_info$class == group1]
                group2_samples <- sample_info$ID[sample_info$class == group2]

                # 提取相似性值
                similarities <- bray_sim[group1_samples, group2_samples]

                if(group1 == group2) {
                    # 组内相似性：取上三角矩阵（避免对角线）
                    similarities_vec <- similarities[upper.tri(similarities)]
                    comparison_type <- "within_group"
                } else {
                    # 组间相似性
                    similarities_vec <- as.vector(similarities)
                    comparison_type <- "between_group"
                }

                if(length(similarities_vec) > 0) {
                    similarity_results <- rbind(similarity_results, data.frame(
                        group1 = group1,
                        group2 = group2,
                        comparison = comparison_type,
                        similarity = similarities_vec,
                        mean_similarity = mean(similarities_vec, na.rm = TRUE)
                    ))
                }
            }
        }

        return(similarity_results)
    }

    # 核心菌群分析函数（按BMI组别分别分析）
    core_microbiome_analysis <- function(processed_data, prevalence = 0.5, abundance = 0.001) {
        relative_abu <- processed_data$relative_abu
        sample_info <- processed_data$sample_info

        core_results <- list()
        bmi_groups <- unique(sample_info$class)

        for(group in bmi_groups) {
            group_samples <- sample_info$ID[sample_info$class == group]
            group_rel <- relative_abu[group_samples, ]

            # 筛选核心菌群
            core_otus <- colnames(group_rel)[
                colMeans(group_rel > abundance) >= prevalence
            ]

            core_results[[group]] <- list(
                core_otus = core_otus,
                core_abundance = colMeans(group_rel[, core_otus, drop = FALSE]),
                sample_size = length(group_samples)
            )
        }

        # 计算共享核心菌群
        all_core_otus <- lapply(core_results, function(x) x$core_otus)
        shared_core <- Reduce(intersect, all_core_otus)

        return(c(list(
            shared_core = shared_core,
            all_core_otus = all_core_otus
        ), core_results))
    }

    # 物种共现网络分析函数（按BMI组别分别分析，过滤孤立节点）
    co_occurrence_network <- function(processed_data, bmi_group, cor_threshold = 0.6) {
        css_otu <- processed_data$css_otu
        sample_info <- processed_data$sample_info

        # 筛选特定BMI组别的样本
        target_samples <- sample_info$ID[sample_info$class == bmi_group]
        target_otu <- css_otu[target_samples, ]

        if(length(target_samples) < 5) {
            return(list(network = NULL, properties = NULL, cor_matrix = NULL))
        }

        # 过滤低丰度OTU
        prevalence <- colMeans(target_otu > 0)
        target_otu <- target_otu[, prevalence > 0.2 & colSums(target_otu) > 0]

        if(ncol(target_otu) < 3) {
            return(list(network = NULL, properties = NULL, cor_matrix = NULL))
        }

        # 计算相关性矩阵
        cor_matrix <- cor(target_otu, method = "spearman")

        # 构建网络
        adj_matrix <- ifelse(abs(cor_matrix) > cor_threshold, abs(cor_matrix), 0)
        diag(adj_matrix) <- 0

        network <- graph.adjacency(
            adj_matrix,
            mode = "undirected",
            weighted = TRUE,
            diag = FALSE
        )

        # 过滤孤立节点（度为0的节点）
        degrees <- degree(network)
        non_isolated_vertices <- which(degrees > 0)

        if(length(non_isolated_vertices) < 2) {
            # 如果过滤后节点数少于2个，返回空结果
            return(list(network = NULL, properties = NULL, cor_matrix = NULL))
        }

        # 提取非孤立节点的子图
        filtered_network <- induced_subgraph(network, vids = non_isolated_vertices)

        # 更新相关性矩阵，只保留非孤立节点对应的行和列
        filtered_cor_matrix <- cor_matrix[non_isolated_vertices, non_isolated_vertices]

        # 计算过滤后的网络属性
        network_properties <- list(
            nodes = vcount(filtered_network),
            edges = ecount(filtered_network),
            density = edge_density(filtered_network),
            transitivity = transitivity(filtered_network),
            average_degree = mean(degree(filtered_network)),
            original_nodes = vcount(network),  # 原始节点数
            isolated_nodes_removed = vcount(network) - vcount(filtered_network)  # 移除的孤立节点数
        )

        return(list(
            network = filtered_network,
            properties = network_properties,
            cor_matrix = filtered_cor_matrix
        ))
    }

    # ALDEx2差异分析函数（BMI组间比较）
    aldex2_analysis <- function(processed_data) {
        otu_matrix <- processed_data$otu_matrix
        sample_info <- processed_data$sample_info
        # 获取唯一的组别
        unique_groups <- unique(sample_info$class)

        # 检查是否有多个BMI组别
        if(length(unique_groups) < 2) {
            return(NULL)
        }

        # 准备ALDEx2输入数据
        conditions <- as.character(sample_info$class)

        # 根据组别数量选择不同的分析方法
        if (length(unique_groups) == 2) {
            # --- 两组比较 ---
            message("检测到两组数据，使用 ALDEx2 t-test...")

            aldex_result <- aldex(
                reads = t(otu_matrix),
                conditions = conditions,  # aldex() 包装函数支持 conditions 参数
                test = "t",
                effect = TRUE,
                include.sample.summary = FALSE,
                verbose = FALSE
            )

        } else {
            # --- 多组比较 (3组及以上) ---
            message("检测到", length(unique_groups), "组数据，使用 ALDEx2 Kruskal-Wallis test...")

            # 使用 test="kw" 进行 Kruskal-Wallis 检验
            # 注意：effect 参数通常只用于 t-test，这里去掉或设为默认
            aldex_result <- aldex(
                reads = t(otu_matrix),
                conditions = conditions, # 参数传给 aldex()，它会自动分发给内部函数
                test = "kw",             # 核心修改：指定使用 Kruskal-Wallis 检验
                denom = "all",           # 标准化方法
                include.sample.summary = FALSE,
                verbose = FALSE
            )
        }

        return(aldex_result)
    }

    # 主分析流程
    main_analysis <- function(otu_table) {
        # 数据预处理和归一化
        processed_data <- preprocess_data(otu_table)

        cat("=== 16s微生物组BMI组间分析开始 ===\n")
        cat("样本类型:", sample_type, "\n")
        cat("样本数量:", nrow(processed_data$sample_info), "\n")
        cat("BMI组别:", paste(unique(processed_data$sample_info$class), collapse = ", "), "\n")
        cat("OTU数量:", ncol(processed_data$otu_matrix), "\n")
        cat("重抽样深度:", min(rowSums(processed_data$otu_matrix)), "\n\n")

        # 1. α多样性分析
        cat("1. 进行α多样性分析...\n")
        alpha_results <- alpha_diversity_analysis(processed_data)

        # 2. β多样性分析
        cat("2. 进行β多样性分析...\n")
        beta_results <- beta_diversity_analysis(processed_data)

        # 3. ALDEx2差异分析
        cat("3. 进行ALDEx2差异分析...\n")
        aldex2_results <- aldex2_analysis(processed_data)

        # 4. 群落相似性分析
        cat("4. 进行群落相似性分析...\n")
        similarity_results <- community_similarity_analysis(processed_data)

        # 5. 核心菌群分析
        cat("5. 进行核心菌群分析...\n")
        core_results <- core_microbiome_analysis(processed_data)

        # 6. 网络分析（按BMI组别）
        cat("6. 进行物种共现网络分析...\n")
        network_results <- list()
        bmi_groups <- unique(processed_data$sample_info$class)
        for(group in bmi_groups) {
            network_results[[group]] <- co_occurrence_network(processed_data, group)
        }

        return(list(
            processed_data = processed_data,
            alpha_diversity = alpha_results,
            beta_diversity = beta_results,
            community_similarity = similarity_results,
            core_microbiome = core_results,
            network_analysis = network_results,
            aldex2_results = aldex2_results
        ))
    }

    save_dir = function(output_dir, name) {
        output_dir = file.path(output_dir, name);
        dir.create(output_dir, showWarnings = FALSE);
        return(output_dir);
    }

    # 结果导出和可视化函数
    export_and_visualize_results <- function(results, output_dir) {
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }

        theme_set(theme_bw(base_size = 12))

        # 导出ALDEx2结果
        if (!is.null(results$aldex2_results)) {
            write.csv(results$aldex2_results, file = file.path(output_dir, "ALDEx2_results.csv"))
        }

        # 1. α多样性结果导出和可视化
        export_alpha_diversity(results, save_dir(output_dir,"alpha_diversity"))

        # 2. β多样性结果导出和可视化
        export_beta_diversity(results, save_dir(output_dir,"beta_diversity"))

        # 3. 群落相似性结果导出和可视化
        export_community_similarity(results, save_dir(output_dir,"community_similarity"))

        # 4. 核心菌群结果导出和可视化
        export_core_microbiome(results, save_dir(output_dir,"core_microbiome"))

        # 5. 网络分析结果导出和可视化
        export_network_analysis(results, save_dir(output_dir,"network_analysis"))

        # 6. 导出预处理数据
        export_processed_data(results, save_dir(output_dir,"pre-processing"))

        cat("所有结果已导出到目录:", output_dir, "\n")
    }

    # α多样性结果导出和可视化
    export_alpha_diversity <- function(results, output_dir) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        alpha_df <- results$alpha_diversity$alpha_df
        test_results <- results$alpha_diversity$test_results

        write.csv(alpha_df, file.path(output_dir, "alpha_diversity_indices.csv"), row.names = FALSE)

        # 导出统计检验结果
        alpha_stats <- data.frame()
        for (index in names(test_results)) {
            if (!is.null(test_results[[index]]$test_result)) {
                test_result <- test_results[[index]]$test_result
                stats_df <- data.frame(
                    Index = index,
                    Test_Type = test_results[[index]]$test_type,
                    P_value = ifelse(!is.null(test_result$p.value), test_result$p.value, NA),
                    Statistic = ifelse(!is.null(test_result$statistic), test_result$statistic, NA),
                    Groups = paste(test_results[[index]]$bmi_groups, collapse = ", ")
                )
                alpha_stats <- rbind(alpha_stats, stats_df)
            }
        }

        if (nrow(alpha_stats) > 0) {
            write.csv(alpha_stats, file.path(output_dir, "alpha_diversity_statistics.csv"), row.names = FALSE)
        }

        # α多样性可视化
        diversity_indices <- c("Shannon", "Simpson", "Richness", "Chao1")

        for (index in diversity_indices) {
            p <- ggplot(alpha_df, aes(x = class, y = .data[[index]], fill = class)) +
                geom_boxplot(alpha = 0.7) +
                geom_jitter(width = 0.2, size = 1) +
                labs(title = paste("Alpha diversity index:", index, "-", sample_type),
                     x = "Group", y = index) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                    base_size = 13
                )

            # 添加显著性标记（如果适用）
            if (index %in% names(test_results) && !is.null(test_results[[index]]$test_result)) {
                p_value <- test_results[[index]]$test_result$p.value
                if ((!is.null(p_value)) && (!is.na(p_value))) {
                    significance <- ifelse(p_value < 0.001, "***",
                                           ifelse(p_value < 0.01, "**",
                                                  ifelse(p_value < 0.05, "*", "NS")))

                    # 简单的显著性标记
                    if (require(ggpubr) && significance != "NS") {
                        # 比较所有组别
                        comparisons <- list()
                        bmi_groups <- unique(alpha_df$class)
                        if(length(bmi_groups) > 1) {
                            for(i in 1:(length(bmi_groups)-1)) {
                                for(j in (i+1):length(bmi_groups)) {
                                    comparisons <- c(comparisons, list(c(bmi_groups[i], bmi_groups[j])))
                                }
                            }
                            p <- p + stat_compare_means(comparisons = comparisons, method = "wilcox.test")
                        }
                    }
                }
            }

            ggsave(file.path(output_dir, paste0("alpha_diversity_", tolower(index), ".png")),
                   p, width = 6, height = 5, dpi = 300)
            ggsave(file.path(output_dir, paste0("alpha_diversity_", tolower(index), ".pdf")),
                   p, width = 6, height = 5)
        }
    }

    # β多样性结果导出和可视化
    export_beta_diversity <- function(results, output_dir) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        beta_results <- results$beta_diversity
        sample_info <- results$processed_data$sample_info

        # 导出距离矩阵
        if(!is.null(beta_results$bray_dist)) {
            write.csv(as.matrix(beta_results$bray_dist),
                      file.path(output_dir, "bray_curtis_distance_matrix.csv"))
        }

        if(!is.null(beta_results$jaccard_dist)) {
            write.csv(as.matrix(beta_results$jaccard_dist),
                      file.path(output_dir, "jaccard_distance_matrix.csv"))
        }

        # PCoA可视化
        if(!is.null(beta_results$pcoa_bray)) {
            pcoa_scores <- as.data.frame(beta_results$pcoa_bray$points)
            colnames(pcoa_scores) <- paste0("PC", 1:ncol(pcoa_scores))
            pcoa_scores$ID <- rownames(pcoa_scores)
            pcoa_scores <- merge(pcoa_scores, sample_info, by = "ID")
            write.csv(pcoa_scores, file.path(output_dir, "pcoa_coordinates_bray.csv"), row.names = FALSE)

            # PCoA图
            p <- ggplot(pcoa_scores, aes(x = PC1, y = PC2, color = class)) +
                geom_point(size = 3, alpha = 0.8) +
                stat_ellipse(aes(fill = class), alpha = 0.1, geom = "polygon") +
                labs(x = paste("PCoA 1 (", round(beta_results$pcoa_bray$eig[1]/sum(beta_results$pcoa_bray$eig)*100, 1), "%)", sep = ""),
                     y = paste("PCoA 2 (", round(beta_results$pcoa_bray$eig[2]/sum(beta_results$pcoa_bray$eig)*100, 1), "%)", sep = ""),
                     title = paste("PCoA based on Bray-Curtis distance -", sample_type)) +
                theme_minimal(base_size = 14)

            ggsave(file.path(output_dir, "pcoa_plot_bray.png"), p, width = 8, height = 6, dpi = 300)
            ggsave(file.path(output_dir, "pcoa_plot_bray.pdf"), p, width = 8, height = 6)
        }

        # 导出PERMANOVA结果
        if(!is.null(beta_results$permanova_bray)) {
            write.csv(as.data.frame(beta_results$permanova_bray),
                      file.path(output_dir, "permanova_bray_results.csv"))
        }
    }

    # 群落相似性结果导出和可视化
    export_community_similarity <- function(results, output_dir) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        similarity_df <- results$community_similarity

        write.csv(similarity_df, file.path(output_dir, "community_similarity_results.csv"), row.names = FALSE)

        # 相似性可视化
        if(nrow(similarity_df) > 0) {
            p <- ggplot(similarity_df, aes(x = comparison, y = similarity, fill = comparison)) +
                geom_boxplot(alpha = 0.7) +
                labs(x = "Comparison Type", y = "Bray-Curtis Similarity",
                    title = paste("Community Similarity Within and Between BMI Groups -", sample_type)) +
                theme(legend.position = "none")

            ggsave(file.path(output_dir, "community_similarity.png"), p, width = 6, height = 5, dpi = 300)
            ggsave(file.path(output_dir, "community_similarity.pdf"), p, width = 6, height = 5)
        }

        # 热图可视化
        if(nrow(similarity_df) > 0) {
            # 计算组内和组间平均相似性
            avg_similarity <- similarity_df %>%
                group_by(group1, group2, comparison) %>%
                summarise(avg_sim = mean(similarity, na.rm = TRUE), .groups = 'drop')

            # 创建热图数据
            heatmap_data <- reshape2::dcast(avg_similarity, group1 ~ group2, value.var = "avg_sim")
            rownames(heatmap_data) <- heatmap_data$group1
            heatmap_matrix <- as.matrix(heatmap_data[, -1])

            # 绘制热图
            p_heatmap <- ggplot(melt(heatmap_matrix), aes(Var2, Var1, fill = value)) +
                geom_tile() +
                scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                     midpoint = 0.5, name = "Similarity") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                labs(title = paste("Community Similarity Heatmap Between BMI Groups -", sample_type),
                     x = "Group", y = "Group")


            ggsave(file.path(output_dir, "community_similarity_heatmap.png"),
                   p_heatmap, width = 6, height = 5, dpi = 300)
            ggsave(file.path(output_dir, "community_similarity_heatmap.pdf"),
                   p_heatmap, width = 6, height = 5)
        }
    }

    # 核心菌群结果导出和可视化
    export_core_microbiome <- function(results, output_dir) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        core_results <- results$core_microbiome

        # 导出核心菌群摘要
        core_summary <- data.frame()
        for(group in names(core_results)) {
            if(!group %in% c("shared_core", "all_core_otus")) {
                core_summary <- rbind(core_summary, data.frame(
                    Group = group,
                    Core_OTUs = length(core_results[[group]]$core_otus),
                    Sample_Size = core_results[[group]]$sample_size
                ))
            }
        }

        write.csv(core_summary, file.path(output_dir, "core_microbiome_summary.csv"), row.names = FALSE)

        # 导出各组的核心菌群
        for(group in names(core_results)) {
            if(!group %in% c("shared_core", "all_core_otus")) {
                group_val = core_results[[group]];

                message(sprintf("group value of %s", group));
                print(group_val);

                try({
                    core_otus_df <- data.frame(
                        OTU_ID = group_val$core_otus,
                        Mean_Abundance = group_val$core_abundance
                    );
                    write.csv(core_otus_df,
                              file.path(output_dir, paste0("core_microbiome_", group, ".csv")),
                              row.names = FALSE)
                });
            }
        }

        # 导出共享核心菌群
        shared_core_df <- data.frame(OTU_ID = core_results$shared_core)
        write.csv(shared_core_df, file.path(output_dir, "shared_core_microbiome.csv"), row.names = FALSE)

        # Venn图（如果有多个组别）
        if(length(core_results$all_core_otus) >= 2) {
            try({
                # 1. 定义一个足够长的颜色向量（这里在原有3个颜色基础上扩充了5个高分文章常用色）
                my_colors <- c("#BD3C29", "#0172B6", "#FF000F", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91")

                # 2. 动态获取分组数量
                n_groups <- length(core_results$all_core_otus)

                # 3. 绘图
                venn_plot <- venn.diagram(
                    x = core_results$all_core_otus,
                    filename = NULL,
                    category.names = names(core_results$all_core_otus),
                    fill = my_colors[1:n_groups],  # 自动根据分组数量提取颜色
                    cat.cex = 1.2,
                    cex = 1.2
                )

                png(file.path(output_dir, "core_microbiome_venn.png"), width = 1000, height = 1000, res = 300)
                grid.draw(venn_plot)
                dev.off()
                png(file.path(output_dir, "core_microbiome_venn.pdf"), width = 8, height = 8, res = 300)
                grid.draw(venn_plot)
                dev.off()
            });
        }

        # 核心菌群丰度柱状图
        # if(length(core_results) > 2) { # 除了shared_core和all_core_otus
        # 准备数据
        core_abundance_data <- data.frame()
        for(group in names(core_results)) {
            # if(!group %in% c("shared_core", "all_core_otus")) {
            group_val = core_results[[group]];

            message(sprintf("group value of '%s':", group));
            print(group_val);

            try({
                group_data <- data.frame(
                    Group = group,
                    OTU_ID = names(group_val$core_abundance),
                    Abundance = group_val$core_abundance
                )
                core_abundance_data <- rbind(core_abundance_data, group_data)
            });
            # }
        }

        # 选择丰度最高的前10个OTU进行可视化
        top_otus <- core_abundance_data %>%
            group_by(OTU_ID) %>%
            summarise(Max_Abundance = max(Abundance), .groups = 'drop') %>%
            arrange(desc(Max_Abundance)) %>%
            head(10) %>%
            pull(OTU_ID)

        plot_data <- core_abundance_data[core_abundance_data$OTU_ID %in% top_otus, ]

        p <- ggplot(plot_data, aes(x = OTU_ID, y = Abundance, fill = Group)) +
            geom_bar(stat = "identity", position = "dodge") +
            labs(title = paste("Core Microbiota Abundance Comparison -", sample_type),
                 x = "OTU", y = "Relative Abundance") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))


        ggsave(file.path(output_dir, "core_microbiome_abundance.png"),
               p, width = 8, height = 6, dpi = 300)
        ggsave(file.path(output_dir, "core_microbiome_abundance.pdf"),
               p, width = 10, height = 6)
        # }
    }

    # 网络分析结果导出和可视化
    export_network_analysis <- function(results, output_dir) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        network_results <- results$network_analysis

        # 导出网络属性
        network_properties <- data.frame()
        for(group in names(network_results)) {
            if(!is.null(network_results[[group]]$properties)) {
                props <- network_results[[group]]$properties
                network_properties <- rbind(network_properties, data.frame(
                    Group = group,
                    Nodes = props$nodes,
                    Edges = props$edges,
                    Density = props$density,
                    Transitivity = props$transitivity,
                    Average_Degree = props$average_degree
                ))
            }
        }

        if(nrow(network_properties) > 0) {
            write.csv(network_properties, file.path(output_dir, "network_properties.csv"), row.names = FALSE)
        }

        # 网络可视化
        for(group in names(network_results)) {
            if(!is.null(network_results[[group]]$network)) {
                network <- network_results[[group]]$network
                if(vcount(network) > 0) {
                    # 绘制网络图
                    net_plot <- ggraph(network, layout = "fr") +
                        geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE,
                                       arrow = arrow(length = unit(2, 'mm')), end_cap = circle(0.05, 'inches')) +
                        geom_node_point(aes(size = degree(network), color = degree(network))) +
                        # geom_node_text(aes(label = name), vjust = 1, hjust = 1, size = 2) +
                        scale_color_gradient(low = "lightblue", high = "darkblue") +
                        scale_size_continuous(range = c(3, 12)) +
                        theme_void() +
                        labs(title = paste("Species co-occurrence network -", group, "-", sample_type)) +
                        # 修改2：增大标题字体大小（例如设置为16）
                        theme(plot.title = element_text(size = 24)) +
                        # 修改3：增大图例（颜色条）的标题字体大小（例如设置为12）
                        guides(color = guide_legend(title.theme = element_text(size = 16)))

                    ggsave(file.path(output_dir, paste0("network_plot_", group, ".png")),
                           net_plot, width = 12, height = 9, dpi = 300)
                    ggsave(file.path(output_dir, paste0("network_plot_", group, ".pdf")),
                           net_plot, width = 12, height = 9)
                }

                # 相关性热图（优化版 - 包含聚类和布局调整）
                if(!is.null(network_results[[group]]$cor_matrix)) {
                    cor_matrix <- network_results[[group]]$cor_matrix
                    if(nrow(cor_matrix) > 2) {
                        try({
                            # 计算距离矩阵用于聚类
                            dist_matrix <- as.dist(1 - abs(cor_matrix))

                            # 进行层次聚类
                            hc_rows <- hclust(dist_matrix, method = "complete")
                            hc_cols <- hclust(dist_matrix, method = "complete")  # 对于对称矩阵，行列聚类相同

                            # 重新排序相关性矩阵
                            ordered_matrix <- cor_matrix[hc_rows$order, hc_cols$order]
                            max_w <- 45;

                            # 保存高质量PDF版本
                            pdf(file.path(output_dir, paste0("correlation_heatmap_", group, ".pdf")),
                                width = max(max_w , nrow(ordered_matrix) * 0.5),  # 根据矩阵大小动态调整宽度
                                height = max(max_w , nrow(ordered_matrix) * 0.5),
                                pointsize = 36)

                            # 绘制带聚类的上三角热图
                            corrplot(ordered_matrix,
                                     method = "color",
                                     type = "upper",
                                     order = "hclust",           # 使用聚类顺序
                                     hclust.method = "complete", # 聚类方法
                                     tl.cex = 0.6,               # 标签字体大小
                                     tl.col = "black",
                                     tl.srt = 45,                # 标签旋转角度
                                     title = paste("Species correlation heatmap -", group, "-", sample_type),
                                     title.cex = 9, # 新增：放大标题字体，默认值为1，可调整此数值
                                     mar = c(2, 2, 3, 2),        # 调整边距防止标题溢出
                                     addgrid.col = "white",      # 添加白色网格线分隔
                                     col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(100), # 自定义颜色
                                     cl.position = "r",          # 颜色条位置在右侧
                                     number.cex = 0,           # 数字大小
                                     # addCoef.col = NA,      # 显示相关系数
                                     diag = TRUE                 # 显示对角线
                            )
                            dev.off()
                        })
                    }
                }

            }
        }
    }

    # 导出预处理数据
    export_processed_data <- function(results, output_dir) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        processed_data <- results$processed_data

        write.csv(processed_data$sample_info, file.path(output_dir, "sample_info.csv"), row.names = FALSE)
        write.csv(processed_data$otu_matrix, file.path(output_dir, "raw_otu_matrix.csv"))
        write.csv(processed_data$rarefied_otu, file.path(output_dir, "rarefied_otu_matrix.csv"))
        write.csv(processed_data$relative_abu, file.path(output_dir, "relative_abundance_matrix.csv"))
        write.csv(processed_data$css_otu, file.path(output_dir, "css_normalized_matrix.csv"))
        write.csv(processed_data$log_otu, file.path(output_dir, "log_transformed_matrix.csv"))
    }

    # 执行分析
    cat("开始执行", sample_type, "样本分析...\n")
    results <- main_analysis(otu_table)

    # 导出结果
    export_and_visualize_results(results, output_dir)

    # 生成摘要报告
    generate_summary_report(results, output_dir)

    cat("=== ", sample_type, "样本分析完成! ===\n")
}

# 生成摘要报告函数
generate_summary_report <- function(results, output_dir) {
    report_file <- file.path(output_dir, "analysis_summary.txt")
    alpha_test <- results$alpha_diversity$test_results$Shannon$test_result

    sink(report_file)
    cat("=== 16S微生物组BMI组间分析报告 ===\n\n")
    cat("分析日期:", date(), "\n\n")

    # 样本信息
    cat("1. 样本信息:\n")
    cat("   总样本数:", nrow(results$processed_data$sample_info), "\n")
    cat("   BMI组别:", paste(unique(results$processed_data$sample_info$class), collapse = ", "), "\n\n")

    # α多样性摘要
    cat("2. α多样性分析:\n")
    if (!is.null(alpha_test)) {
        if (is.data.frame(alpha_test[[1]])) {
            cat("   Shannon指数检验p值:", round(alpha_test[[1]]$`Pr(>F)`[1], 4), "\n")
        } else {
            cat("   Shannon指数检验p值:", round(alpha_test$p.value[1], 4), "\n")
        }
    }
    cat("\n")

    # β多样性摘要
    cat("3. β多样性分析:\n")
    if (!is.null(results$beta_diversity$permanova_bray)) {
        cat("   Bray-Curtis PERMANOVA p值:", round(results$beta_diversity$permanova_bray$`Pr(>F)`[1], 4), "\n")
    }
    cat("\n")

    # 群落相似性摘要
    cat("4. 群落相似性分析:\n")
    similarity_df <- results$community_similarity
    if(nrow(similarity_df) > 0) {
        within_sim <- similarity_df[similarity_df$comparison == "within_group", "similarity"]
        between_sim <- similarity_df[similarity_df$comparison == "between_group", "similarity"]
        if(length(within_sim) > 0) {
            cat("   组内平均相似性:", round(mean(within_sim, na.rm = TRUE), 4), "\n")
        }
        if(length(between_sim) > 0) {
            cat("   组间平均相似性:", round(mean(between_sim, na.rm = TRUE), 4), "\n")
        }
    }
    cat("\n")

    # 核心菌群摘要
    cat("5. 核心菌群分析:\n")
    core_results <- results$core_microbiome
    if(!is.null(core_results$shared_core)) {
        cat("   共享核心菌群OTU数量:", length(core_results$shared_core), "\n")
        for(group in names(core_results)) {
            if(!group %in% c("shared_core", "all_core_otus")) {
                cat("   ", group, "组核心菌群OTU数量:", length(core_results[[group]]$core_otus), "\n")
            }
        }
    }
    cat("\n")

    # 网络分析摘要
    cat("6. 网络分析:\n")
    network_results <- results$network_analysis
    for(group in names(network_results)) {
        if(!is.null(network_results[[group]]$properties)) {
            props <- network_results[[group]]$properties
            cat("   ", group, "组网络节点数:", props$nodes, ", 边数:", props$edges, "\n")
        }
    }
    cat("\n")

    sink()
    cat("分析报告已保存到:", report_file, "\n")
}


