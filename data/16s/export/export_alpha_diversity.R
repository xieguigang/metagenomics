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