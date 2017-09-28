library(tximport)

expression_data_output_path <- "../preprocessed/expression_data"
sample_description_output_path <- "../preprocessed/sample_description"

# I'm removing the first row, because in this case, it is just a row with gene_id in the first column and "normalized count" in all the rest of the columns
TCGA_expression_data <- read.csv("../input_files/TCGA/LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", sep="\t", row.names=1, header=T, stringsAsFactors=F)[-1,]

# I'm removing all genes with no annotation. Also, I'ma also removing SLC35E2 because there are two entries of it, which raises confusion.
removal_index <- !(grepl("[?]", rownames(TCGA_expression_data)) | grepl("SLC35E2", rownames(TCGA_expression_data)))

TCGA_known_genes_expression <- TCGA_expression_data[removal_index,]
rownames(TCGA_known_genes_expression) <- sapply(rownames(TCGA_known_genes_expression), function(gene) {
    strsplit(gene, "[|]")[[1]][1]
})

# The following line makes the colnames into the same format as the SAMPLE_ID column in the clinical data
colnames(TCGA_known_genes_expression) <- sapply(colnames(TCGA_known_genes_expression), function(cname) {
    paste(strsplit(cname, "[.]03[A/B]")[[1]][1], ".03", sep="")
})
rownames(TCGA_known_genes_expression) <- make.names(rownames(TCGA_known_genes_expression))

# The table starts from row 6. The other rows just have descriptions of the columns
TCGA_clinical_data <- read.csv("../input_files/TCGA/data_clinical.txt", skip=5, sep="\t")

TCGA_cytogenetic_risk_info <- as.data.frame(t(sapply(colnames(TCGA_known_genes_expression), function(sample_id) {
index <- grepl(sample_id, TCGA_clinical_data$SAMPLE_ID)
cytogenetic_risk <- if (TCGA_clinical_data$RISK_CYTO[index] == "N.D.") { "unknown" } else { tolower(TCGA_clinical_data$RISK_CYTO[index]) }
c(sample_id, cytogenetic_risk)
})))

colnames(TCGA_cytogenetic_risk_info) <- c("sample_id", "cytogenetic_risk")

write.csv(TCGA_known_genes_expression, file.path(expression_data_output_path, "LAML_TCGA.csv"))
write.csv(TCGA_cytogenetic_risk_info, file.path(sample_description_output_path, "LAML_TCGA_sample_description.csv"))

# SRP050272 processing

SRP050272_input_path <- "../input_files/SRP050272"
tx2gene <- read.csv(file.path(SRP050272_input_path, "tx2gene.csv"))
files <- file.path(SRP050272_input_path, "salmon_output", list.files(path=file.path(SRP050272_input_path, "salmon_output"), pattern="SRX"), "quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
SRP050272_expression <- txi$abundance

colnames(SRP050272_expression) <- list.files(path=file.path(SRP050272_input_path, "salmon_output"), pattern="SRX")
rownames(SRP050272_expression) <- make.names(rownames(SRP050272_expression))

SRP050272_cytogenetic_risk_info <- data.frame(sample_id=colnames(SRP050272_expression), cytogenetic_risk=rep("intermediate", ncol(SRP050272_expression)))

write.csv(SRP050272_expression, file.path(expression_data_output_path, "SRP050272.csv"))
write.csv(SRP050272_cytogenetic_risk_info, file.path(sample_description_output_path, "SRP050272_sample_description.csv"))
