#' 核心菌群分析函数（按BMI组别分别分析）
core_microbiome_analysis <- function(processed_data, prevalence = 0.5, abundance = 0.001) {
    relative_abu <- processed_data$relative_abu
    sample_info <- processed_data$sample_info

    core_results <- list()
    bmi_groups <- unique(sample_info$class)

    for(group in bmi_groups) {
        group_samples <- sample_info$ID[sample_info$class == group]
        group_rel <- relative_abu[group_samples, ]

        # 筛选核心菌群
        core_otus <- colnames(group_rel)[
            colMeans(group_rel > abundance) >= prevalence
        ]

        core_results[[group]] <- list(
            core_otus = core_otus,
            core_abundance = colMeans(group_rel[, core_otus, drop = FALSE]),
            sample_size = length(group_samples)
        )
    }

    # 计算共享核心菌群
    all_core_otus <- lapply(core_results, function(x) x$core_otus)
    shared_core <- Reduce(intersect, all_core_otus)

    return(c(list(
        shared_core = shared_core,
        all_core_otus = all_core_otus
    ), core_results))
}