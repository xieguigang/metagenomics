

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
