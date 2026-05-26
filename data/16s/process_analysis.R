

process_analysis = function(data, sample_info, output_dir) {
    # 读取数据
    sample_info = read.csv(sample_info, row.names = 1, check.names = FALSE);
    otu_table <- read.csv(data, row.names = 1, check.names = FALSE);
    dir.create(output_dir, showWarnings = FALSE);

    otu_table[,"taxonomy"] = NULL;
    otu_table[,"Tree"] = NULL;
    otu_table[,"OTUs"] = NULL;
    otu_table = t(otu_table);

    otu_table = as.data.frame(otu_table[rownames(sample_info),]);

    print("view of the data for run analysis:");
    str(sample_info);
    str(otu_table);

    message("准备数据");
    print(head(sample_info));
    otu_table = cbind(ID = rownames(sample_info), class = sample_info$sample_info, otu_table);
    sample_type = "metagenomics";
    rownames(otu_table) = rownames(sample_info);

    message("数据将要输出到这个文件夹：",output_dir);
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE);
    write.csv(otu_table, file = file.path(output_dir,"data.csv"), row.names = FALSE);

    # 执行分析
    cat("开始执行", sample_type, "样本分析...\n")
    results <- main_analysis(otu_table)

    # 导出结果
    export_and_visualize_results(results, output_dir)

    # 生成摘要报告
    generate_summary_report(results, output_dir)

    cat("=== ", sample_type, "样本分析完成! ===\n")
}
