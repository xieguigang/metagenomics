#' 物种共现网络分析函数（按BMI组别分别分析，过滤孤立节点）
co_occurrence_network <- function(processed_data, bmi_group, cor_threshold = 0.6) {
    css_otu <- processed_data$css_otu
    sample_info <- processed_data$sample_info

    # 筛选特定BMI组别的样本
    target_samples <- sample_info$ID[sample_info$class == bmi_group]
    target_otu <- css_otu[target_samples, ]

    if(length(target_samples) < 5) {
        return(list(network = NULL, properties = NULL, cor_matrix = NULL))
    }

    # 过滤低丰度OTU
    prevalence <- colMeans(target_otu > 0)
    target_otu <- target_otu[, prevalence > 0.2 & colSums(target_otu) > 0]

    if(ncol(target_otu) < 3) {
        return(list(network = NULL, properties = NULL, cor_matrix = NULL))
    }

    # 计算相关性矩阵
    cor_matrix <- cor(target_otu, method = "spearman")

    # 构建网络
    adj_matrix <- ifelse(abs(cor_matrix) > cor_threshold, abs(cor_matrix), 0)
    diag(adj_matrix) <- 0

    network <- graph.adjacency(
        adj_matrix,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
    )

    # 过滤孤立节点（度为0的节点）
    degrees <- degree(network)
    non_isolated_vertices <- which(degrees > 0)

    if(length(non_isolated_vertices) < 2) {
        # 如果过滤后节点数少于2个，返回空结果
        return(list(network = NULL, properties = NULL, cor_matrix = NULL))
    }

    # 提取非孤立节点的子图
    filtered_network <- induced_subgraph(network, vids = non_isolated_vertices)

    # 更新相关性矩阵，只保留非孤立节点对应的行和列
    filtered_cor_matrix <- cor_matrix[non_isolated_vertices, non_isolated_vertices]

    # 计算过滤后的网络属性
    network_properties <- list(
        nodes = vcount(filtered_network),
        edges = ecount(filtered_network),
        density = edge_density(filtered_network),
        transitivity = transitivity(filtered_network),
        average_degree = mean(degree(filtered_network)),
        original_nodes = vcount(network),  # 原始节点数
        isolated_nodes_removed = vcount(network) - vcount(filtered_network)  # 移除的孤立节点数
    )

    return(list(
        network = filtered_network,
        properties = network_properties,
        cor_matrix = filtered_cor_matrix
    ))
}