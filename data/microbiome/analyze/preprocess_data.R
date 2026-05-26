#' 数据预处理和归一化函数
preprocess_data <- function(otu_table) {
    message("开始处理OTU数据表");

    # 提取样本信息
    sample_info <- otu_table[, 1:2]
    otu_data <- otu_table[, -c(1:2)]
    scale <- 100000000

    for(name in colnames(otu_data ) ) {            
        v = as.numeric( otu_data[, name]) * scale; 
        v = as.integer(v); 

        if (all(v == 0)) {
            v = as.integer( runif(length(v)) * 10)
        }

        otu_data[, name] = v;
    }

    print(head(otu_data ));

    # 转换OTU数据为数值矩阵
    otu_matrix <- as.matrix(otu_data)
    rownames(otu_matrix) <- sample_info$ID

    # 数据归一化步骤
    cat("进行重抽样标准化...\n")
    min_depth <- min(rowSums(otu_matrix))
    rarefied_otu <- rrarefy(otu_matrix, sample = min_depth)

    # 相对丰度标准化
    relative_abu <- sweep(otu_matrix, 1, rowSums(otu_matrix), "/")
    relative_abu[is.na(relative_abu)] <- 0

    # CSS标准化
    css_normalize <- function(otu_mat) {
        lib_sizes <- rowSums(otu_mat)
        css_factors <- lib_sizes / quantile(lib_sizes, 0.75)
        css_norm <- sweep(otu_mat, 1, css_factors, "/")
        return(css_norm)
    }

    css_otu <- css_normalize(otu_matrix)

    # 对数转换
    log_otu <- log10(otu_matrix + 1)

    return(list(
        sample_info = sample_info,
        otu_matrix = otu_matrix,
        rarefied_otu = rarefied_otu,
        relative_abu = relative_abu,
        css_otu = css_otu,
        log_otu = log_otu
    ))
}