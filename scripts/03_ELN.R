source("common_functions.R")
library(data.table)

nclusters <- 5
h1_min <- 25
h1_max <- 300
h1_step <- 25
h2_min <- 25
h2_max <- 300
h2_step <- 25

expression_data_folder <- "../preprocessed/expression_data"
sample_description_folder <- "../preprocessed/sample_description"
training_expression_data_file <- "GSE6891.csv"
training_sample_description_file <- "GSE6891_sample_description.csv"
gene_pair_output_folder <- "../output_files/gene_pairs"

training_norm_data <- data.frame(fread(file.path(expression_data_folder, training_expression_data_file)), row.names=1)
training_sample_description <- data.frame(fread(file.path(sample_description_folder, training_sample_description_file)), row.names=1)

TCGA_expression_data <- data.frame(fread(file.path(expression_data_folder, "LAML_TCGA.csv"), stringsAsFactors=F), row.names=1)
# The expression is in character and keeps refusing to be coereced into numeric
TCGA_genes_mean_expression <- apply(TCGA_expression_data, 1, function(x) { mean(as.numeric(as.character(x))) })
TCGA_genes <- rownames(TCGA_expression_data)[TCGA_genes_mean_expression>10] # Getting well expressed genes only

unknown_samples <- training_sample_description$GEO_accession[training_sample_description$cytogenetic_risk=="unknown"]
good_samples <- training_sample_description$GEO_accession[training_sample_description$cytogenetic_risk=="good"]
intr_samples <- training_sample_description$GEO_accession[training_sample_description$cytogenetic_risk=="intermediate"]
poor_samples <- training_sample_description$GEO_accession[training_sample_description$cytogenetic_risk=="poor"]

subset_training_norm_data <- training_norm_data[ ,(!(colnames(training_norm_data) %in% unknown_samples))]

training_labels <- rep("", ncol(subset_training_norm_data))
training_labels[colnames(subset_training_norm_data) %in% good_samples] <- "good"
training_labels[colnames(subset_training_norm_data) %in% intr_samples] <- "intr_poor"
training_labels[colnames(subset_training_norm_data) %in% poor_samples] <- "intr_poor"

simple_training_labels <- rep(0, length(training_labels))
simple_training_labels[training_labels=="good"] <- 1

full_rnk_matrix <- apply(subset_training_norm_data, 2, function(x) { rank(as.numeric(as.character(x))) })
rownames(full_rnk_matrix) <- rownames(subset_training_norm_data)
rnk_matrix <- full_rnk_matrix[rownames(full_rnk_matrix) %in% TCGA_genes, ]
rnk_mean <- apply(rnk_matrix, 1, FUN=function(x) {mean(as.numeric(as.character(x)))})
rnk_sd <- apply(rnk_matrix, 1, FUN=function(x) {sd(as.numeric(as.character(x)))})
rnk_df <- as.data.frame(cbind(rnk_mean, rnk_sd))
quants <- quantile(rnk_df$rnk_mean)

max_iteration_df <- generate_iteration_df(rnk_df, quants, h1_max, h2_max)
gene_pair_rnk_indicator_matrix <- generate_gene_pair_rnk_indicator_matrix(rnk_matrix, max_iteration_df, nclusters)
write.csv(gene_pair_rnk_indicator_matrix, paste("../preprocessed/gene_pair_rnk_indicator_matrix_h1_", h1_max, "h2_", h2_max, ".csv", sep=""))
colsd <- apply(gene_pair_rnk_indicator_matrix, 2, function(x) { sd(as.numeric(as.character(x))) })

h1_vals <- seq(h1_min, h1_max, h1_step)
h2_vals <- seq(h2_min, h2_max, h2_step)
h1_h2_vals_df <- expand.grid(h1_vals, h2_vals)
colnames(h1_h2_vals_df) <- c("h1_vals", "h2_vals")

if (FALSE %in% (h1_h2_vals_df$h1_vals <= h1_max)) {
    print("h1 val greater than that was used to generate the gpri_matrix detected. Caution")
}
if (FALSE %in% (h1_h2_vals_df$h2_vals <= h2_max)) {
   print("h2 val greater than that was used to generate the gpri_matrix detected. Caution")
}
all_gene_pairs <- apply(h1_h2_vals_df, 1, function(h_vals) {
    h1_val <- h_vals[1]
    h2_val <- h_vals[2]
    print(paste("Getting gene pairs for h1=", h1_val, "h2=", h2_val))
    iteration_df <- generate_iteration_df(rnk_df, quants, h1_val, h2_val)
    subset_gpri_matrix <- generate_subset_gene_pair_rnk_indicator_matrix(gene_pair_rnk_indicator_matrix[,colsd!=0], iteration_df)
    gene_pairs <- ELN_get_gene_pairs(subset_gpri_matrix, simple_training_labels)
})

names(all_gene_pairs) <- apply(h1_h2_vals_df, 1, function(h_vals) {
    h1_val <- h_vals[1]
    h2_val <- h_vals[2]
    paste("h1", h1_val, "h2", h2_val, sep="_")
})


lapply(names(all_gene_pairs), function(h_val_combination) {
    write.csv(all_gene_pairs[[h_val_combination]], file.path(gene_pair_output_folder, paste("ELN_GSE6891derived_", h_val_combination, ".csv", sep="")))
})
