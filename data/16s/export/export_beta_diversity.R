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