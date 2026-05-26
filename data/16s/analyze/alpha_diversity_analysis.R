#' α多样性分析函数（BMI组间比较）
alpha_diversity_analysis <- function(processed_data) {
    rarefied_otu <- processed_data$rarefied_otu
    sample_info <- processed_data$sample_info

    # 计算α多样性指数
    shannon <- vegan::diversity(rarefied_otu, index = "shannon")
    simpson <- vegan::diversity(rarefied_otu, index = "simpson")
    richness <- rowSums(rarefied_otu > 0)
    chao1 <- estimateR(rarefied_otu)[2, ]

    alpha_df <- data.frame(
        ID = rownames(rarefied_otu),
        Shannon = shannon,
        Simpson = simpson,
        Richness = richness,
        Chao1 = chao1
    )

    # 合并样本信息
    alpha_df <- merge(alpha_df, sample_info, by = "ID")

    print("alpha_diversity");
    print(head(alpha_df));

    # BMI组间比较
    results <- list()
    diversity_indices <- c("Shannon", "Simpson", "Richness", "Chao1")

    for(index in diversity_indices) {
        # 检查是否有多个BMI组别
        bmi_groups <- unique(alpha_df$class)
        if(length(bmi_groups) < 2) {
            results[[index]] <- list(
                test_result = NULL,
                test_type = "Insufficient groups",
                message = "需要至少两个BMI组别进行比较"
            )
            next
        }

        # 正态性检验（对每个组别）
        normality_p <- c();

        for(group in bmi_groups) {
            group_data <- alpha_df[alpha_df$class == group, index];

            if(length(group_data) > 3) { # 需要足够样本进行正态性检验
                print("view group data for shapiro.test:");
                print(group_data);

                if (all(group_data == group_data[0])) {
                    # all 'x' values are identical
                    normality_p <- c(normality_p, 0)
                } else {
                    normality_test <- shapiro.test(group_data)
                    normality_p <- c(normality_p, normality_test$p.value)
                }
            }
        }

        # 如果所有组别数据正态且组别数=2，使用t检验，否则使用ANOVA或Kruskal-Wallis
        if(length(bmi_groups) == 2 && all(normality_p > 0.05)) {
            # 两组比较，数据正态，使用t检验
            group1_data <- alpha_df[alpha_df$class == bmi_groups[1], index]
            group2_data <- alpha_df[alpha_df$class == bmi_groups[2], index]
            test_result <- t.test(group1_data, group2_data);
            test_type <- "Student's t-test"
        } else if(length(bmi_groups) == 2) {
            # 两组比较，数据非正态，使用Mann-Whitney检验
            group1_data <- alpha_df[alpha_df$class == bmi_groups[1], index]
            group2_data <- alpha_df[alpha_df$class == bmi_groups[2], index]
            test_result <- wilcox.test(group1_data, group2_data)
            test_type <- "Mann-Whitney U test"
        } else if(length(bmi_groups) > 2 && all(normality_p > 0.05)) {
            # 多组比较，数据正态，使用ANOVA
            anova_model <- aov(as.formula(paste(index, "~ class")), data = alpha_df)
            test_result <- summary(anova_model)
            test_type <- "ANOVA"
        } else if(length(bmi_groups) > 2) {
            # 多组比较，数据非正态，使用Kruskal-Wallis检验
            test_result <- kruskal.test(as.formula(paste(index, "~ class")), data = alpha_df)
            test_type <- "Kruskal-Wallis test"
        }

        results[[index]] <- list(
            test_result = test_result,
            test_type = test_type,
            groups = bmi_groups,
            normality_p = normality_p
        )
    }

    return(list(alpha_df = alpha_df, test_results = results))
}