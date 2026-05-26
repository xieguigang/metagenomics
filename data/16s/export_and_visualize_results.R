#' 结果导出和可视化函数
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