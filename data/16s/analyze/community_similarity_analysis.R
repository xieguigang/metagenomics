#' 群落相似性分析函数（BMI组内和组间相似性）
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