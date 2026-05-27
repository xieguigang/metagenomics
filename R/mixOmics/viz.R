const viz = function(data_raw) {
    library(tidyverse);
    library(igraph);
    library(ggraph);
    library(ComplexHeatmap);
    library(circlize);
    library(ggrepel);

    if (!is.data.frame(data_raw)) {
        # ==========================================
        # 1. 数据读取与预处理
        # ==========================================
        # 读取数据（处理可能的UTF-8 BOM头问题）
        data_raw <- read.csv(data_raw, check.names = FALSE, fileEncoding = "UTF-8-BOM");
    }

    # 过滤掉 association 为 None 的行
    df <- data_raw %>% 
    filter(association != "None") %>% 
    mutate(
        association = factor(association, levels = c("Monotonic", "NonMonotonic")),
        # 将spearman-rho的符号转化为正负相关分类，用于网络边颜色
        corr_dir = ifelse(`spearman-rho` >= 0, "Positive", "Negative"),
        # 计算pvalue的-log10值，用于散点大小映射，避免无限值加一个小常数
        neg_log10_pval = -log10(pvalue + 1e-20)
    );

    # ==========================================
    # 2. 微生物与代谢物的关联网络图
    # ==========================================
    # 构建igraph对象
    g <- graph_from_data_frame(df, directed = FALSE);

    # 提取节点属性并分类
    nodes <- tibble(name = V(g)$name) %>% mutate(type = ifelse(name %in% df$source, "Microbe", "Metabolite"));

    V(g)$type <- nodes$type;

    # 绘制网络图
    p_net <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(
        width = score,
        color = corr_dir,
        linetype = association
    ),
    alpha = 0.7
    ) +
    scale_edge_color_manual(values = c("Positive" = "#D64541", "Negative" = "#4B77BE")) +
    scale_edge_linetype_manual(values = c("Monotonic" = "solid", "NonMonotonic" = "dashed")) +
    scale_edge_width(range = c(0.5, 2.5), guide = "legend") +
    geom_node_point(aes(shape = type, color = type), size = 6, fill = "white", stroke = 1.5) +
    scale_shape_manual(values = c("Microbe" = 21, "Metabolite" = 24)) +
    scale_color_manual(values = c("Microbe" = "#2C3E50", "Metabolite" = "#F39C12", "Positive" = "#D64541", "Negative" = "#4B77BE"), guide = "none") +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5, max.overlaps = 20) +
    theme_void() +
    labs(
        title = "Microbe-Metabolite Association Network",
        edge_width = "Score",
        edge_color = "Correlation",
        edge_linetype = "Association Type",
        shape = "Node Type"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"));

    # 导出网络图
    ggsave("01_Network.pdf", p_net, width = 10, height = 8, dpi = 300);
    ggsave("01_Network.png", p_net, width = 10, height = 8, dpi = 300);


    # ==========================================
    # 3. 微生物与代谢物的关联热图
    # ==========================================
    # 转换为宽格式矩阵，缺失值填充为0
    mat_wide <- df %>%
    select(source, target, score) %>%
    pivot_wider(names_from = target, values_from = score, values_fill = 0) %>%
    column_to_rownames(var = "source") %>%
    as.matrix();

    # 设定颜色映射（由于score通常是正值，采用从浅到深的连续颜色）
    col_fun <- colorRamp2(c(0, quantile(mat_wide, 0.75), max(mat_wide)), c("#FCFDBF", "#FE9F6D", "#B73779"));

    # 绘制ComplexHeatmap
    pdf("02_Heatmap.pdf", width = 10, height = 8);
    png("02_Heatmap.png", width = 10, height = 8, units = "in", res = 300);

    ht <- Heatmap(mat_wide,
                name = "Score",
                col = col_fun,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_row_dend = TRUE,
                show_column_dend = TRUE,
                row_names_side = "left",
                column_names_side = "top",
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(title = "Score"),
                column_title = "Microbe-Metabolite Association Heatmap",
                column_title_gp = gpar(fontsize = 12, fontface = "bold")
    );
    draw(ht);
    dev.off();


    # ==========================================
    # 4. 微生物与代谢物的关联曼哈顿图
    # ==========================================
    # 寻找每个微生物score Top1的代谢物用于打标签
    top1_labels <- df %>%
    group_by(source) %>%
    slice_max(order_by = score, n = 1) %>%
    ungroup();

    # 绘制曼哈顿jitter散点图
    p_manh <- ggplot(df, aes(x = source, y = score)) +
    geom_jitter(aes(size = neg_log10_pval, shape = association, fill = association), 
                width = 0.25, alpha = 0.8, color = "black", stroke = 0.5) +
    scale_shape_manual(values = c("Monotonic" = 21, "NonMonotonic" = 22)) +
    scale_fill_manual(values = c("Monotonic" = "#E74C3C", "NonMonotonic" = "#3498DB")) +
    geom_text_repel(
        data = top1_labels,
        aes(label = target),
        size = 3, 
        box.padding = 0.5,
        max.overlaps = 15,
        segment.color = "grey50"
    ) +
    theme_bw() +
    labs(
        title = "Manhattan Plot of Microbe-Metabolite Associations",
        x = "Microbe (Source)",
        y = "Association Score",
        size = expression(-log[10](pvalue)),
        shape = "Association Type",
        fill = "Association Type"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        panel.grid.major.x = element_blank()
    );

    # 导出曼哈顿图
    ggsave("03_Manhattan.pdf", p_manh, width = 10, height = 6, dpi = 300);
    ggsave("03_Manhattan.png", p_manh, width = 10, height = 6, dpi = 300);


    # ==========================================
    # 5. 微生物与代谢物的关联得分的火山图
    # ==========================================
    # 找出Score降序排列的Top5关联
    top5_labels <- df %>%
    slice_max(order_by = score, n = 5);

    # 绘制火山图
    p_vol <- ggplot(df, aes(x = `spearman-rho`, y = score)) +
    geom_point(aes(size = neg_log10_pval, shape = association, fill = association), 
                alpha = 0.8, color = "black", stroke = 0.5) +
    scale_shape_manual(values = c("Monotonic" = 21, "NonMonotonic" = 22)) +
    scale_fill_manual(values = c("Monotonic" = "#E74C3C", "NonMonotonic" = "#3498DB")) +
    geom_text_repel(
        data = top5_labels,
        aes(label = paste0(source, " -> ", target)),
        size = 3,
        box.padding = 0.5,
        max.overlaps = 10,
        segment.color = "grey50"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_bw() +
    labs(
        title = "Volcano Plot of Association Score vs Spearman Rho",
        x = expression(Spearman~rho),
        y = "Association Score",
        size = expression(-log[10](pvalue)),
        shape = "Association Type",
        fill = "Association Type"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
    );

    # 导出火山图
    ggsave("04_Volcano.pdf", p_vol, width = 9, height = 7, dpi = 300);
    ggsave("04_Volcano.png", p_vol, width = 9, height = 7, dpi = 300);
}