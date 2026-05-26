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