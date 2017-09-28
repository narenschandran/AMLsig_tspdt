library(data.table)

expression_data_folder <- "../preprocessed/expression_data"
sample_description_folder <- "../preprocessed/sample_description"
rnk_matrix_folder <- "../preprocessed/rnk_matrix_files"
rnk_matrix_sample_labels_folder <- "../preprocessed/rnk_matrix_sample_labels"

expression_data_files <- list.files(expression_data_folder)
sample_description_files <- list.files(sample_description_folder)


normalized_expression_data <- lapply(expression_data_files, function(expression_data_file) {
    print(paste("Reading in", expression_data_file))
    data.frame(fread(file.path(expression_data_folder, expression_data_file), stringsAsFactors=F), row.names=1)
})

names(normalized_expression_data) <- sapply(expression_data_files, function(expression_data_file) {
    strsplit(expression_data_file, "[.]")[[1]][1]
})

sample_description <- lapply(sample_description_files, function(sample_description_file) {
    print(paste("Reading in", sample_description_file))
    data.frame(fread(file.path(sample_description_folder, sample_description_file), stringsAsFactors=F), row.names=1)
})

names(sample_description) <- sapply(sample_description_files, function(sample_description_file) {
    strsplit(sample_description_file, "_sample_description")[[1]][1]
})

# names(normalized_expression_data) == names(sample_description) # Sanity check
lapply(names(normalized_expression_data), function(dataset) {
    norm_data <- normalized_expression_data[[dataset]]
    cytogenetic_risk_info <- sample_description[[dataset]]
    if (grepl("GSE", dataset)) {
        unknown_samples <- cytogenetic_risk_info$GEO_accession[(cytogenetic_risk_info$cytogenetic_risk=="unknown") | (cytogenetic_risk_info$cytogenetic_risk=="skip")]
    } else {
        unknown_samples <- cytogenetic_risk_info$sample_id[cytogenetic_risk_info$cytogenetic_risk=="unknown"]
    }
    unknown_sample_index <- colnames(norm_data) %in% unknown_samples
    subset_norm_data <- norm_data[,(!(unknown_sample_index))]
    subset_cytogenetic_risk <- sapply(colnames(subset_norm_data), function(cname) {
        if (grepl("GSE", dataset)) {
            cytogenetic_risk_info$cytogenetic_risk[cytogenetic_risk_info$GEO_accession==cname]
        } else {
            cytogenetic_risk_info$cytogenetic_risk[cytogenetic_risk_info$sample_id==cname]
        }
    })
    subset_sample_label <- sapply(subset_cytogenetic_risk, function(crisk) {
        if (crisk=="good") {
            "good"
        } else if (crisk=="intermediate") {
            "intr_poor"
        } else if (crisk=="poor") {
            "intr_poor"
        } else {
            "check_for_errors"
        }
    })
    sample_labels <- data.frame(sample_id=names(subset_sample_label), label=subset_sample_label)
    rnk_matrix <- apply(subset_norm_data, 2, function(exp_data) {
        rank(as.numeric(as.character(exp_data)))
    })
    rownames(rnk_matrix) <- rownames(subset_norm_data)
    print(paste("Writing rnk_matrix and sample_labels for", dataset))
    write.csv(sample_labels, file.path(rnk_matrix_sample_labels_folder, paste(dataset, "_sample_labels.csv", sep="")))
    write.csv(rnk_matrix, file.path(rnk_matrix_folder, paste(dataset, "_rnk_matrix.csv", sep="")))
})

common_genes <- Reduce(intersect, lapply(normalized_expression_data, function(x) { rownames(x) }))
write.csv(data.frame(common_genes=common_genes), "../preprocessed/common_genes.csv")
