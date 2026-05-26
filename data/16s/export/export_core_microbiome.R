# 核心菌群结果导出和可视化
export_core_microbiome <- function(results, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    core_results <- results$core_microbiome

    # 导出核心菌群摘要
    core_summary <- data.frame()
    for(group in names(core_results)) {
        if(!group %in% c("shared_core", "all_core_otus")) {
            core_summary <- rbind(core_summary, data.frame(
                Group = group,
                Core_OTUs = length(core_results[[group]]$core_otus),
                Sample_Size = core_results[[group]]$sample_size
            ))
        }
    }

    write.csv(core_summary, file.path(output_dir, "core_microbiome_summary.csv"), row.names = FALSE)

    # 导出各组的核心菌群
    for(group in names(core_results)) {
        if(!group %in% c("shared_core", "all_core_otus")) {
            group_val = core_results[[group]];

            message(sprintf("group value of %s", group));
            print(group_val);

            try({
                core_otus_df <- data.frame(
                    OTU_ID = group_val$core_otus,
                    Mean_Abundance = group_val$core_abundance
                );
                write.csv(core_otus_df,
                            file.path(output_dir, paste0("core_microbiome_", group, ".csv")),
                            row.names = FALSE)
            });
        }
    }

    # 导出共享核心菌群
    shared_core_df <- data.frame(OTU_ID = core_results$shared_core)
    write.csv(shared_core_df, file.path(output_dir, "shared_core_microbiome.csv"), row.names = FALSE)

    # Venn图（如果有多个组别）
    if(length(core_results$all_core_otus) >= 2) {
        try({
            # 1. 定义一个足够长的颜色向量（这里在原有3个颜色基础上扩充了5个高分文章常用色）
            my_colors <- c("#BD3C29", "#0172B6", "#FF000F", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91")

            # 2. 动态获取分组数量
            n_groups <- length(core_results$all_core_otus)

            # 3. 绘图
            venn_plot <- venn.diagram(
                x = core_results$all_core_otus,
                filename = NULL,
                category.names = names(core_results$all_core_otus),
                fill = my_colors[1:n_groups],  # 自动根据分组数量提取颜色
                cat.cex = 1.2,
                cex = 1.2
            )

            png(file.path(output_dir, "core_microbiome_venn.png"), width = 1000, height = 1000, res = 300)
            grid.draw(venn_plot)
            dev.off()
            png(file.path(output_dir, "core_microbiome_venn.pdf"), width = 8, height = 8, res = 300)
            grid.draw(venn_plot)
            dev.off()
        });
    }

    # 核心菌群丰度柱状图
    # if(length(core_results) > 2) { # 除了shared_core和all_core_otus
    # 准备数据
    core_abundance_data <- data.frame()
    for(group in names(core_results)) {
        # if(!group %in% c("shared_core", "all_core_otus")) {
        group_val = core_results[[group]];

        message(sprintf("group value of '%s':", group));
        print(group_val);

        try({
            group_data <- data.frame(
                Group = group,
                OTU_ID = names(group_val$core_abundance),
                Abundance = group_val$core_abundance
            )
            core_abundance_data <- rbind(core_abundance_data, group_data)
        });
        # }
    }

    # 选择丰度最高的前10个OTU进行可视化
    top_otus <- core_abundance_data %>%
        group_by(OTU_ID) %>%
        summarise(Max_Abundance = max(Abundance), .groups = 'drop') %>%
        arrange(desc(Max_Abundance)) %>%
        head(10) %>%
        pull(OTU_ID)

    plot_data <- core_abundance_data[core_abundance_data$OTU_ID %in% top_otus, ]

    p <- ggplot(plot_data, aes(x = OTU_ID, y = Abundance, fill = Group)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = paste("Core Microbiota Abundance Comparison -", sample_type),
                x = "OTU", y = "Relative Abundance") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


    ggsave(file.path(output_dir, "core_microbiome_abundance.png"),
            p, width = 8, height = 6, dpi = 300)
    ggsave(file.path(output_dir, "core_microbiome_abundance.pdf"),
            p, width = 10, height = 6)
    # }
}