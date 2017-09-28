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

gene_pairs <- load_gene_pairs(gene_pairs_path, "ELN_GSE6891derived_h1_150_h2_150.csv", common_genes=common_genes)[[1]]

rpart_result <- derive_uncut_decision_tree(training_data, training_labels$label, gene_pairs, nclusters)
rframe <- rpart_result$frame$var
used_gene_pairs <- do.call("rbind", lapply(rframe[rframe!="<leaf>"], function(gp) {
        g1 = strsplit(as.character(gp), "_")[[1]][1]
        g2 = strsplit(as.character(gp), "_")[[1]][2]
        data.frame(gene1=g1, gene2=g2)
    }))
used_gene_pairs <- used_gene_pairs[!duplicated(used_gene_pairs),]
rpart_result <- derive_uncut_decision_tree(training_data, training_labels$label, used_gene_pairs, ncluster)
cp_val <- 0.05
print(paste("Pruning using cp_val", cp_val))
rpart_pruned <- prune(rpart_result, cp=cp_val) 
gene_pairs_frame <- unique(rpart_pruned$frame$var)
ngene_pairs <- length(gene_pairs_frame[gene_pairs_frame!="<leaf>"])
final_df <- decision_tree_validation_generate_dataset_results(rpart_pruned, test_data, test_labels, used_gene_pairs, nclusters)
colnames(final_df) <- c("Total Intermediate and Poor prognosis samples", "Predicted Intermediate and Poor prognosis samples", "Correctly predicted Intermediate and Poor prognosis samples",
"Total Good prognosis samples", "Predicted Good prognosis samples", "Correctly predicted Good prognosis samples")


training_prediction <- as.character(predict(rpart_pruned, type="class"))
training_predicted_intrpoor <- sum(training_prediction=="intr_poor")
training_total_intrpoor <- sum(training_labels$label=="intr_poor")
training_correct_intrpoor <- sum(training_prediction=="intr_poor" & training_labels$label=="intr_poor")
training_predicted_good <- sum(training_prediction=="good")
training_correct_good <- sum(training_prediction=="good" & training_labels$label=="good")
training_total_good <- sum(training_labels$label=="good")

final_df_w_training <- rbind(final_df, c(training_total_intrpoor, training_predicted_intrpoor, training_correct_intrpoor, training_total_good, training_predicted_good, training_correct_good))
rownames(final_df_w_training) <- c(rownames(final_df), "GSE6891_training_dataset")


write.csv(final_df, "../output_files/ELN_h1_150_h2_150_performance_validation.csv")
write.csv(final_df_w_training, "../output_files/ELN_h1_150_h2_150_performance_validation_incl_training.csv")
