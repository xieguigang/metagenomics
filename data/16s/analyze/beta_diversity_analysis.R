#' β多样性分析函数（BMI组间比较）
beta_diversity_analysis <- function(processed_data) {
    rarefied_otu <- processed_data$rarefied_otu
    sample_info <- processed_data$sample_info

    # 计算距离矩阵
    bray_dist <- vegdist(rarefied_otu, method = "bray")
    jaccard_dist <- vegdist(rarefied_otu > 0, method = "jaccard")

    # PCoA分析
    pcoa_bray <- cmdscale(bray_dist, k = 3, eig = TRUE)
    pcoa_jaccard <- cmdscale(jaccard_dist, k = 3, eig = TRUE)

    # PERMANOVA检验（BMI组间比较）
    if(length(unique(sample_info$class)) > 1) {
        permanova_bray <- adonis2(bray_dist ~ class,
                                    data = sample_info, permutations = 999)

        permanova_jaccard <- adonis2(jaccard_dist ~ class,
                                        data = sample_info, permutations = 999)
    } else {
        permanova_bray <- NULL
        permanova_jaccard <- NULL
    }

    return(list(
        bray_dist = bray_dist,
        jaccard_dist = jaccard_dist,
        pcoa_bray = pcoa_bray,
        pcoa_jaccard = pcoa_jaccard,
        permanova_bray = permanova_bray,
        permanova_jaccard = permanova_jaccard
    ))
}