library(data.table)
library(rpart)

nclusters <- 5
source("common_functions.R")
rnk_matrix_folder <- "../preprocessed/rnk_matrix_files"
label_folder <- "../preprocessed/rnk_matrix_sample_labels/"
training_data_file <- "GSE6891_rnk_matrix.csv"
training_label_file <- "GSE6891_sample_labels.csv"
gene_pairs_path <- "../output_files/gene_pairs"
common_genes_file <- "../preprocessed/common_genes.csv"
gene_pairs_files <- list.files(gene_pairs_path)

test_data_files <- c(
"GSE10358_rnk_matrix.csv",
"GSE12417_GPL570_rnk_matrix.csv",
"GSE13159_AMLonly_rnk_matrix.csv",
"GSE15434_rnk_matrix.csv",
"GSE16015_rnk_matrix.csv",
"GSE22845_rnk_matrix.csv",
"GSE30285_rnk_matrix.csv",
"GSE39730_rnk_matrix.csv",
"GSE61804_rnk_matrix.csv",
"LAML_TCGA_rnk_matrix.csv",
"SRP050272_rnk_matrix.csv"
)


test_label_files <- c(
"GSE10358_sample_labels.csv",
"GSE12417_GPL570_sample_labels.csv",
"GSE13159_AMLonly_sample_labels.csv",
"GSE15434_sample_labels.csv",
"GSE16015_sample_labels.csv",
"GSE22845_sample_labels.csv",
"GSE30285_sample_labels.csv",
"GSE39730_sample_labels.csv",
"GSE61804_sample_labels.csv",
"LAML_TCGA_sample_labels.csv",
"SRP050272_sample_labels.csv"
)

common_genes <- as.character(read.csv(common_genes_file, row.names=1)$common_genes)

training_data <- data.frame(fread(file.path(rnk_matrix_folder, training_data_file), stringsAsFactors=F), row.names=1)
training_labels <- data.frame(fread(file.path(label_folder, training_label_file), stringsAsFactors=F), row.names=1)

test_data <- lapply(test_data_files, function(test_data_file) {
    data.frame(fread(file.path(rnk_matrix_folder, test_data_file), stringsAsFactors=F), row.names=1)
})

names(test_data) <- sapply(test_data_files, function(test_data_file) {
    strsplit(test_data_file, "_rnk")[[1]][1]
})

test_labels <- lapply(test_label_files, function(test_label_file) {
    data.frame(fread(file.path(label_folder, test_label_file), stringsAsFactors=F), row.names=1)
})

names(test_labels) <- sapply(test_label_files, function(test_label_file) {
    strsplit(test_label_file, "_sample")[[1]][1]
})

gene_pairs <- load_gene_pairs(gene_pairs_path, gene_pairs_files, common_genes)

all_gene_pairs_all_cp_vals_results <- lapply(names(gene_pairs), function(gp) {
    print(gp)
    gene_pair_df <- gene_pairs[[gp]]
    rpart_result <- derive_uncut_decision_tree(training_data, training_labels$label, gene_pair_df, nclusters)
    rframe <- rpart_result$frame$var
    used_gene_pairs <- do.call("rbind", lapply(rframe[rframe!="<leaf>"], function(gp) {
        g1 = strsplit(as.character(gp), "_")[[1]][1]
        g2 = strsplit(as.character(gp), "_")[[1]][2]
        data.frame(gene1=g1, gene2=g2)
    }))
    used_gene_pairs <- used_gene_pairs[!duplicated(used_gene_pairs),]
    rpart_result <- derive_uncut_decision_tree(training_data, training_labels$label, used_gene_pairs, ncluster)
    cp_vals <- c(round((as.numeric(as.character(printcp(rpart_result)[-1,"CP"])) + 0.0001), 5))
    print(paste("cp_vals to use for pruning:", paste(cp_vals, collapse=", ")))
    all_cp_vals_predictions <- lapply(cp_vals, function(cp_val) {
        print(gp)
        print(paste("Pruning using cp_val", cp_val))
        rpart_pruned <- prune(rpart_result, cp=cp_val) 
        gene_pairs_frame <- unique(rpart_pruned$frame$var)
        ngene_pairs <- length(gene_pairs_frame[gene_pairs_frame!="<leaf>"])
        cbind(decision_tree_validation(rpart_pruned, test_data, test_labels, used_gene_pairs, nclusters), ngene_pairs)
    })
    cbind(cp_vals, do.call("rbind", all_cp_vals_predictions))
})

names(all_gene_pairs_all_cp_vals_results) <- gene_pairs_files

final_results <- do.call("rbind", lapply(names(all_gene_pairs_all_cp_vals_results), function(hvals) {
    final_df <- cbind(rep(hvals, nrow(all_gene_pairs_all_cp_vals_results[[hvals]])), all_gene_pairs_all_cp_vals_results[[hvals]])
    colnames(final_df) <- c("gene_pairs_file", "cp_vals", "intrpoor_sensitivity", "intrpoor_specificity", "good_sensitivity", "good_specificity", "num_gene_pairs")
    final_df
}))

write.csv(final_results, "../output_files/final_results.csv")
final_results$intrpoor_sensitivity <- round(as.numeric(as.character(final_results$intrpoor_sensitivity))*100, 2)
final_results$good_sensitivity <- round(as.numeric(as.character(final_results$good_sensitivity))*100, 2)
final_results$intrpoor_specificity <- round(as.numeric(as.character(final_results$intrpoor_specificity))*100, 2)
final_results$good_specificity <- round(as.numeric(as.character(final_results$good_specificity))*100, 2)


write.csv(final_results, "../output_files/final_sensitivity_specificity_results.csv")
