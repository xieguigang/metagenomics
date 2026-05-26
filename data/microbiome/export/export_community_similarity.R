#' 群落相似性结果导出和可视化
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