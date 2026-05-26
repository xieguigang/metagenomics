
#' ALDEx2差异分析函数（BMI组间比较）
aldex2_analysis <- function(processed_data) {
    otu_matrix <- processed_data$otu_matrix
    sample_info <- processed_data$sample_info
    # 获取唯一的组别
    unique_groups <- unique(sample_info$class)

    # 检查是否有多个BMI组别
    if(length(unique_groups) < 2) {
        return(NULL)
    }

    # 准备ALDEx2输入数据
    conditions <- as.character(sample_info$class)

    # 根据组别数量选择不同的分析方法
    if (length(unique_groups) == 2) {
        # --- 两组比较 ---
        message("检测到两组数据，使用 ALDEx2 t-test...")

        aldex_result <- aldex(
            reads = t(otu_matrix),
            conditions = conditions,  # aldex() 包装函数支持 conditions 参数
            test = "t",
            effect = TRUE,
            include.sample.summary = FALSE,
            verbose = FALSE
        )

    } else {
        # --- 多组比较 (3组及以上) ---
        message("检测到", length(unique_groups), "组数据，使用 ALDEx2 Kruskal-Wallis test...")

        # 使用 test="kw" 进行 Kruskal-Wallis 检验
        # 注意：effect 参数通常只用于 t-test，这里去掉或设为默认
        aldex_result <- aldex(
            reads = t(otu_matrix),
            conditions = conditions, # 参数传给 aldex()，它会自动分发给内部函数
            test = "kw",             # 核心修改：指定使用 Kruskal-Wallis 检验
            denom = "all",           # 标准化方法
            include.sample.summary = FALSE,
            verbose = FALSE
        )
    }

    return(aldex_result)
}