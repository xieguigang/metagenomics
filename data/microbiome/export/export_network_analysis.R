#' 网络分析结果导出和可视化
export_network_analysis <- function(results, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    network_results <- results$network_analysis

    # 导出网络属性
    network_properties <- data.frame()
    for(group in names(network_results)) {
        if(!is.null(network_results[[group]]$properties)) {
            props <- network_results[[group]]$properties
            network_properties <- rbind(network_properties, data.frame(
                Group = group,
                Nodes = props$nodes,
                Edges = props$edges,
                Density = props$density,
                Transitivity = props$transitivity,
                Average_Degree = props$average_degree
            ))
        }
    }

    if(nrow(network_properties) > 0) {
        write.csv(network_properties, file.path(output_dir, "network_properties.csv"), row.names = FALSE)
    }

    # 网络可视化
    for(group in names(network_results)) {
        if(!is.null(network_results[[group]]$network)) {
            network <- network_results[[group]]$network
            if(vcount(network) > 0) {
                # 绘制网络图
                net_plot <- ggraph(network, layout = "fr") +
                    geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE,
                                    arrow = arrow(length = unit(2, 'mm')), end_cap = circle(0.05, 'inches')) +
                    geom_node_point(aes(size = degree(network), color = degree(network))) +
                    # geom_node_text(aes(label = name), vjust = 1, hjust = 1, size = 2) +
                    scale_color_gradient(low = "lightblue", high = "darkblue") +
                    scale_size_continuous(range = c(3, 12)) +
                    theme_void() +
                    labs(title = paste("Species co-occurrence network -", group, "-", sample_type)) +
                    # 修改2：增大标题字体大小（例如设置为16）
                    theme(plot.title = element_text(size = 24)) +
                    # 修改3：增大图例（颜色条）的标题字体大小（例如设置为12）
                    guides(color = guide_legend(title.theme = element_text(size = 16)))

                ggsave(file.path(output_dir, paste0("network_plot_", group, ".png")),
                        net_plot, width = 12, height = 9, dpi = 300)
                ggsave(file.path(output_dir, paste0("network_plot_", group, ".pdf")),
                        net_plot, width = 12, height = 9)
            }

            # 相关性热图（优化版 - 包含聚类和布局调整）
            if(!is.null(network_results[[group]]$cor_matrix)) {
                cor_matrix <- network_results[[group]]$cor_matrix
                if(nrow(cor_matrix) > 2) {
                    try({
                        # 计算距离矩阵用于聚类
                        dist_matrix <- as.dist(1 - abs(cor_matrix))

                        # 进行层次聚类
                        hc_rows <- hclust(dist_matrix, method = "complete")
                        hc_cols <- hclust(dist_matrix, method = "complete")  # 对于对称矩阵，行列聚类相同

                        # 重新排序相关性矩阵
                        ordered_matrix <- cor_matrix[hc_rows$order, hc_cols$order]
                        max_w <- 45;

                        # 保存高质量PDF版本
                        pdf(file.path(output_dir, paste0("correlation_heatmap_", group, ".pdf")),
                            width = max(max_w , nrow(ordered_matrix) * 0.5),  # 根据矩阵大小动态调整宽度
                            height = max(max_w , nrow(ordered_matrix) * 0.5),
                            pointsize = 36)

                        # 绘制带聚类的上三角热图
                        corrplot(ordered_matrix,
                                    method = "color",
                                    type = "upper",
                                    order = "hclust",           # 使用聚类顺序
                                    hclust.method = "complete", # 聚类方法
                                    tl.cex = 0.6,               # 标签字体大小
                                    tl.col = "black",
                                    tl.srt = 45,                # 标签旋转角度
                                    title = paste("Species correlation heatmap -", group, "-", sample_type),
                                    title.cex = 9, # 新增：放大标题字体，默认值为1，可调整此数值
                                    mar = c(2, 2, 3, 2),        # 调整边距防止标题溢出
                                    addgrid.col = "white",      # 添加白色网格线分隔
                                    col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(100), # 自定义颜色
                                    cl.position = "r",          # 颜色条位置在右侧
                                    number.cex = 0,           # 数字大小
                                    # addCoef.col = NA,      # 显示相关系数
                                    diag = TRUE                 # 显示对角线
                        )
                        dev.off()
                    })
                }
            }

        }
    }
}