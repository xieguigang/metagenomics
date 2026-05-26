
#' 生成摘要报告函数
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


