#' 主分析流程
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