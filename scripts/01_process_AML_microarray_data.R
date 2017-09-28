raw_files_dir <- "../input_files/RAW_files"
array_info_file_path <- "../gse_array_info.csv"
expression_file_output_path <- "../preprocessed/expression_data"

library(oligo)
library(hgu133plus2hsentrezg.db)
library(huex10sthsentrezg.db)

array_info <- read.csv(array_info_file_path, colClasses="character")

probe_annotations <- list()
probe_annotations[["U133plus2"]] <- toTable(hgu133plus2hsentrezgSYMBOL)
probe_annotations[["HuEx10ST"]]  <- toTable(huex10sthsentrezgSYMBOL)

gse_RAW_folders <- list.files(raw_files_dir)

input_folders_paths <- sapply(gse_RAW_folders, function(gse_RAW_folder) {
    file.path(raw_files_dir, gse_RAW_folder)
})
names(input_folders_paths) <- gse_RAW_folders

normalized_expression_data <- lapply(gse_RAW_folders, function(gse_RAW_folder) {
    probe_annotation <- probe_annotations[[array_info$Array[array_info$GSE_folders==gse_RAW_folder]]]
    pkgname <- as.character(array_info$pkgname[array_info$GSE_folders==gse_RAW_folder])
    cel_file_paths <- file.path(input_folders_paths[[gse_RAW_folder]], list.celfiles(input_folders_paths[[gse_RAW_folder]]))
    norm_data <- exprs(rma(read.celfiles(cel_file_paths, pkgname=pkgname)))
    annotation <- sapply(rownames(norm_data), function(probe_id) {
                    if (probe_id %in% probe_annotation$probe_id) {
                        probe_annotation$symbol[probe_annotation$probe_id==probe_id]
                    } else {
                        probe_id
                    }
    })
    rownames(norm_data) <- annotation
    duplicated_genes <- rownames(norm_data)[duplicated(rownames(norm_data))]
    norm_data <- norm_data[!(rownames(norm_data) %in% duplicated_genes),]
    if (gse_RAW_folder=="GSE6891_RAW") {
        print(paste("Reformatting colnames for", gse_RAW_folder)) 
        cnames1 <- sapply(colnames(norm_data)[1:460], function(cname) { strsplit(cname, "[.]")[[1]][1] })
        cnames2 <- sapply(colnames(norm_data)[461:536], function(cname) { strsplit(cname, "[_]")[[1]][1] })
        cnames <- c(cnames1, cnames2)
    } else if (gse_RAW_folder=="GSE10358_RAW") {
        print(paste("Reformatting colnames for", gse_RAW_folder))
        cnames1 <- sapply(colnames(norm_data)[1:25], function(cname1) {
            strsplit(cname1, "_")[[1]][1]
        })
        cnames2 <- sapply(colnames(norm_data)[26:209], function(cname2) {
            strsplit(cname2, "[.]")[[1]][1]
        })
        cnames3 <- sapply(colnames(norm_data)[210:300], function(cname3) {
            strsplit(cname3, "_")[[1]][1]
        })
        cnames <- c(cnames1, cnames2, cnames3)
    } else if (gse_RAW_folder %in% c("GSE12417_GPL570_RAW", "GSE13159_AMLonly_RAW", "GSE15434_RAW", "GSE16015_RAW", "GSE22845_RAW", "GSE30285_RAW")) {
        print(paste("Reformatting colnames for", gse_RAW_folder))
        cnames <- sapply(colnames(norm_data), function(cname) { strsplit(cname, "[.]")[[1]][1] })
    } else if (gse_RAW_folder %in% c("GSE39730_RAW", "GSE61804_RAW")) {
        print(paste("Reformatting colnames for", gse_RAW_folder))
        cnames <- sapply(colnames(norm_data), function(cname) { strsplit(cname, "_")[[1]][1] })
    } else {
        print("No information on how to reformat colnames given. Using default names")
        cnames <- colnames(norm_data)
    }
    colnames(norm_data) <- cnames
    outfile_name <- if (gse_RAW_folder=="GSE12417_GPL570_RAW") { 
        file.path(expression_file_output_path, "GSE12417_GPL570.csv")
    } else {
        file.path(expression_file_output_path, paste(strsplit(gse_RAW_folder, "_RAW")[[1]][1], ".csv", sep=""))
    }
    write.csv(norm_data, outfile_name) 
    norm_data
    rownames(norm_data) <- make.names(rownames(norm_data))
})
names(normalized_expression_data) <- sapply(gse_RAW_folders, function(base_name) {
    strsplit(base_name, "_")[[1]][1]
})
soft_file_directory <- "../input_files/soft_files"
sample_description_output_path <- "../preprocessed/sample_description"

library(GEOquery)

source("common_functions.R")
source("../input_files/prognosis/GSE10358_prognosis.R")
source("../input_files/prognosis/GSE30285_prognosis.R")
source("../input_files/prognosis/GSE61804_prognosis.R")
source("../input_files/prognosis/GSE13159_prognosis.R")


# There is an error which goes like:
# not all columns named in 'colClasses' exist
# It's apparently not a problem. Should go when updated
# Source: https://support.bioconductor.org/p/86596/

soft_files <- list.files(soft_file_directory)

sample_description_files <- lapply(soft_files, function(soft_file) {

soft_file_path <- file.path(soft_file_directory, soft_file)
gse <- load_gse(soft_file_path)

if (soft_file=="GSE6891_family.soft") {
    gse_extract <- sample_description_extractor(gse, "characteristics_ch1")

    risk <- lapply(gse_extract, function(extract) {
        extract_condition <- grepl("risk", extract)
        if (TRUE %in% extract_condition) {
            risk_string <- extract[grepl("risk", extract)]
            strsplit(risk_string, "risk: ")[[1]][2]
        } else {
            "unknown"
        }
    })

    cytogenetic_risk <- rep("", length(risk))
    cytogenetic_risk[risk=="cytogenetic good"]  <- "good"
    cytogenetic_risk[risk=="cytogenetic intermediate"] <- "intermediate"
    cytogenetic_risk[risk=="cytogenetic poor"] <- "poor"
    cytogenetic_risk[risk=="unknown"] <- "unknown"

    cytogenetic_risk_info <- data.frame(GEO_accession=names(risk), cytogenetic_risk=unlist(cytogenetic_risk))

} else if (soft_file=="GSE13159_family.soft") {
    gse_extract <- t(sample_description_extractor(gse, "characteristics_ch1"))
    karyotype_data <- data.frame(GEO_accession=rownames(gse_extract), karyotype=gse_extract[,2])

    cytogenetic_risk <- rep("", nrow(karyotype_data))
    skip_index <- !grepl("AML", karyotype_data$karyotype)
    cytogenetic_risk[karyotype_data$karyotype %in% GSE13159_good_karyotypes] <- "good"
    cytogenetic_risk[karyotype_data$karyotype %in% GSE13159_intr_karyotypes] <- "intermediate"
    cytogenetic_risk[karyotype_data$karyotype %in% GSE13159_poor_karyotypes] <- "poor"
    cytogenetic_risk[skip_index] <- "skip"
    cytogenetic_risk_info <- cbind(karyotype_data, cytogenetic_risk)

} else if (soft_file=="GSE10358_family.soft") {
    print("Processing GSE10358.soft")
    gse_extract <- sample_description_extractor(gse, "characteristics_ch1")
    description_table1 <- lapply(gse_extract[1:280], function(extract) {
        risk_string <- extract[grepl("cytogenetic", extract)]
        strsplit(risk_string, ": ")[[1]][2]
    })
    karyotype_data1 <- data.frame(GEO_accession=names(unlist(description_table1)), karyotype=unlist(description_table1))
    description_table2 <- lapply(gse_extract[281:304], function(extract) {
        risk_string <- extract[grepl("karyotype", extract)]
        strsplit(risk_string, ": ")[[1]][2]
    })
    karyotype_data2 <- data.frame(GEO_accession=names(unlist(description_table2)), karyotype=unlist(description_table2))
    karyotype_data <- rbind(karyotype_data1, karyotype_data2)
    cytogenetic_risk <- rep("", nrow(karyotype_data))
    cytogenetic_risk[karyotype_data$karyotype %in% GSE10358_good_karyotypes] <- "good"
    cytogenetic_risk[karyotype_data$karyotype %in% GSE10358_intr_karyotypes] <- "intermediate"
    cytogenetic_risk[karyotype_data$karyotype %in% GSE10358_poor_karyotypes] <- "poor"
    cytogenetic_risk[karyotype_data$karyotype %in% GSE10358_not_determined] <- "unknown"
    cytogenetic_risk[karyotype_data$karyotype %in% GSE10358_unknown_karyotypes] <- "unknown"
    cytogenetic_risk[cytogenetic_risk==""] <- "unknown"
    cytogenetic_risk_info <- cbind(karyotype_data, cytogenetic_risk) 

} else if (soft_file=="GSE12417_family.soft") {
    platform_id_extract <- sample_description_extractor(gse, "platform_id")
    GPL570_index <- grepl("GPL570", platform_id_extract)
    full_gse_extract <- sample_description_extractor(gse, "characteristics_ch1")
    karyotype <- sapply(full_gse_extract[GPL570_index], function(karyotype_string) {
        strsplit(karyotype_string, "; ")[[1]][1]
    })
    karyotype_data <- data.frame(GEO_accession=names(platform_id_extract)[GPL570_index], karyotype=karyotype)
    cytogenetic_risk <- rep("", nrow(karyotype_data))
    cytogenetic_risk[grepl("normal karyotype", karyotype_data$karyotype)] <- "intermediate"
    cytogenetic_risk[grepl("MDS", karyotype_data$karyotype)] <- "skip"
    cytogenetic_risk_info <- cbind(karyotype_data, cytogenetic_risk)

} else if (soft_file=="GSE15434_family.soft") {
    gse_extract <- sample_description_extractor(gse, "characteristics_ch1")
    karyotype_data <- data.frame(GEO_accession=colnames(gse_extract), karyotype=sapply(as.character(gse_extract[5,]), function(k) {
        strsplit(k, "diagnosis: ")[[1]][2]
    }))
    cytogenetic_risk <- rep("", nrow(karyotype_data))
    cytogenetic_risk[karyotype_data$karyotype=="AML with normal karyotype"] <- "intermediate"
    cytogenetic_risk_info <- cbind(karyotype_data, cytogenetic_risk)

} else if (soft_file=="GSE16015_family.soft") {
    gse_extract <- sample_description_extractor(gse, "characteristics_ch1")
    karyotype_data <- data.frame(GEO_accession=as.character(colnames(gse_extract)), karyotype=sapply(as.character(gse_extract[2,]), function(k) {
        strsplit(k, "karyotype: ")[[1]][2]
    }))
    cytogenetic_risk <- rep("", nrow(karyotype_data))
    cytogenetic_risk[karyotype_data$karyotype=="normal karyotype"] <- "intermediate"
    cytogenetic_risk[karyotype_data$karyotype=="other chromosomal aberrations"] <- "unknown"
    cytogenetic_risk_info <- cbind(karyotype_data, cytogenetic_risk)

} else if (soft_file=="GSE22845_family.soft") {
    gse_extract <- sample_description_extractor(gse, "title")
    karyotype_data <- data.frame(GEO_accession=names(gse_extract), karyotype=gse_extract)
    cytogenetic_risk <- rep("", nrow(karyotype_data))
    cytogenetic_risk[grepl("normal karyotype", karyotype_data$karyotype)] <- "intermediate"
    cytogenetic_risk[cytogenetic_risk==""] <- "unknown"
    cytogenetic_risk_info <- cbind(karyotype_data, cytogenetic_risk)

} else if (soft_file=="GSE30285_family.soft") {

    gse_extract <- sample_description_extractor(gse, "characteristics_ch1")
    karyotype_data <- data.frame(GEO_accession=as.character(colnames(gse_extract)), karyotype=sapply(as.character(gse_extract[6,]), function(k) {
        strsplit(k, "abnormality: ")[[1]][2]
    }))
    cytogenetic_risk <- rep("", nrow(karyotype_data))
    cytogenetic_risk[karyotype_data$karyotype %in% GSE30285_good_karyotypes] <- "good"
    cytogenetic_risk[karyotype_data$karyotype %in% GSE30285_intr_karyotypes] <- "intermediate"
    cytogenetic_risk_info <- cbind(karyotype_data, cytogenetic_risk)

} else if (soft_file=="GSE39730_family.soft") {
    gse_extract <- sample_description_extractor(gse, "characteristics_ch1")
    karyotype_data <- data.frame(GEO_accession=colnames(gse_extract), karyotype=sapply(gse_extract[2,], function(karyotype_string) { strsplit(karyotype_string, "karyotype: ")[[1]][2] }))
    cytogenetic_risk <- rep("", nrow(karyotype_data))
    cytogenetic_risk[karyotype_data$karyotype=="complex karyotype"] <- "poor"
    cytogenetic_risk_info <- cbind(karyotype_data, cytogenetic_risk)

} else if (soft_file=="GSE61804_family.soft") {
    gse_extract <- sample_description_extractor(gse, "characteristics_ch1")
    karyotype_data <- data.frame(GEO_accession=as.character(colnames(gse_extract)), karyotype=sapply(as.character(gse_extract[4,]), function(k) {
        strsplit(k, "condition: ")[[1]][2]
    }))
    skip_index <- !grepl("AML", karyotype_data$karyotype)
    cytogenetic_risk <- rep("", nrow(karyotype_data))
    cytogenetic_risk[karyotype_data$karyotype %in% GSE61804_good_karyotypes] <- "good"
    cytogenetic_risk[karyotype_data$karyotype %in% GSE61804_poor_karyotypes] <- "poor"
    cytogenetic_risk[karyotype_data$karyotype %in% GSE61804_intr_karyotypes] <- "intermediate"
    cytogenetic_risk[skip_index] <- "skip"
    cytogenetic_risk[cytogenetic_risk==""] <- "unknown"
    cytogenetic_risk_info <- cbind(karyotype_data, cytogenetic_risk)

} else {
    print(paste("There are no insturctions on how to process", soft_file))
}
output_file_path <- if (soft_file=="GSE12417_family.soft") {
    file.path(sample_description_output_path, "GSE12417_GPL570_sample_description.csv")
} else if (soft_file=="GSE13159_family.soft") {
  file.path(sample_description_output_path, "GSE13159_AMLonly_sample_description.csv")
} else {
    file.path(sample_description_output_path, paste(strsplit(soft_file, "_")[[1]][1], "_sample_description.csv", sep=""))
}
write.csv(cytogenetic_risk_info, output_file_path)
})

names(sample_description_files) <- soft_files
