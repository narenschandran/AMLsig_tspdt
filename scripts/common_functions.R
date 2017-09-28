##### EXTRACT FUNCTIONS ########


# There is an error which goes like:
# not all columns named in 'colClasses' exist
# It's apparently not a problem. Should go when updated
# Source: https://support.bioconductor.org/p/86596/

library(GEOquery)


load_gse <- function(gse_file) {
    print(paste("Reading file", gse_file))
    gse_time <- system.time(gse <- getGEO(filename=gse_file))
    print(paste("Finished reading file", gse_file))
    print(gse_time)
    gse
}


sample_description_extractor <- function(gse, extraction_feature) {
    gsmlist <- GSMList(gse)
    sapply(gsmlist, function(gsm) {
        Meta(gsm)[[extraction_feature]]
    })
}

##### ITERATION DF ########

generate_iteration_df <- function(rnk_df, quants, h1, h2) {

    quantile_list1 <- lapply(1:4, function(n) {
        if (n==1) {
            index <- (rnk_df$rnk_mean >= quants[n]) & (rnk_df$rnk_mean < quants[n+1])
        } else {
            index <- (rnk_df$rnk_mean > quants[n]) & (rnk_df$rnk_mean <= quants[n+1])
        }
        subset_df <- rnk_df[index,]
        subset_df <- subset_df[order(subset_df$rnk_sd),]
        nr <- nrow(subset_df)
        return_df <- as.data.frame(cbind(rownames(subset_df)[(nr-h1+1):nr], rownames(subset_df)[(nr-h1+1):nr]))
        colnames(return_df) <- c("g1", "g2")
        return_df
    })


    iteration_df1 <- do.call('rbind', lapply(quantile_list1, function(q) {
        iteration_df <- expand.grid(q$g1, q$g2)
        colnames(iteration_df) <- c("g1", "g2")
        iteration_df
        })) 



    quantile_list2 <- lapply(1:4, function(n) {
        if (n==1) {
                index <- (rnk_df$rnk_mean >= quants[n]) & (rnk_df$rnk_mean < quants[n+1])
        } else {
                index <- (rnk_df$rnk_mean > quants[n]) & (rnk_df$rnk_mean <= quants[n+1])
        }
        subset_df <- rnk_df[index,]
        subset_df <- subset_df[order(subset_df$rnk_sd),]
        nr <- nrow(subset_df)
        return_df <- as.data.frame(cbind(rownames(subset_df)[1:h2], rownames(subset_df)[(nr-h2+1):nr]))
        colnames(return_df) <- c("g1", "g2")
        return_df
    })

    iteration_df2 <- do.call('rbind', lapply(quantile_list2, function(q) {
        iteration_df <- expand.grid(q$g1, q$g2)
        colnames(iteration_df) <- c("g1", "g2")
        iteration_df
    }))

    iteration_df <- rbind(iteration_df1, iteration_df2)

    iteration_df[,1] <- as.character(iteration_df[,1])
    iteration_df[,2] <- as.character(iteration_df[,2])
    iteration_df
}


### GENE PAIR MATRIX ###

library(parallel)


generate_gene_pair_rnk_indicator_matrix <- function(rnk_matrix, iteration_df, nclusters) {
    # Clusters functionality is to extend into parallelization. It is turned off for now because I can't reliably handle long runs of programs with parallelization
    print("Making the gene pair indicator matrix")
#    print(paste("Initiating", nclusters, "clusters"))
#    cl <- makeCluster(nclusters, type="FORK")
#    gene_pair_rnk_matrix <- parApply(cl, iteration_df, 1, function(genes) {
     gene_pair_rnk_matrix <- apply(iteration_df, 1, function(genes) {
        print(genes)
        g1_index <- rownames(rnk_matrix) %in% genes[1]
        print(which(g1_index))
        g2_index <- rownames(rnk_matrix) %in% genes[2]
        g1rnk <- rnk_matrix[g1_index,]
        g2rnk <- rnk_matrix[g2_index,]
        rnk_vector <- g1rnk > g2rnk # Gene 1 has higer rank, and hence, lower expression than gene 2
    })
#    stopCluster(cl)
    print("Created gene pair rnk matrix.")
    colnames(gene_pair_rnk_matrix) <- mapply(FUN=function(g1, g2) {
        paste(g1, g2, sep="_")
    }, iteration_df[,1], iteration_df[,2])
    gene_pair_rnk_matrix * 1
}

generate_subset_gene_pair_rnk_indicator_matrix <- function(gpri_matrix, iteration_df) {
    print("Subsetting the gene pair indicator matrix")
    gene_pairs <- paste(as.character(iteration_df[,1]), as.character(iteration_df[,2]), sep="_")
    subset_gpri_matrix <- gpri_matrix[ ,colnames(gpri_matrix) %in% gene_pairs]
}
### ELASTIC NET REGRESSION ###

library(glmnet)

ELN_get_gene_pairs <- function(gpri_matrix, simple_labels) {
    # gpri_matrix is gene_pair_rnk_indicator_matrix
    # simple_labels have 0 or 1 for class
    ELN_model <- cv.glmnet(gpri_matrix, simple_labels, intercept=F)
    coef_matrix <- coef(ELN_model, s="lambda.min")
    gene_pairs <- coef_matrix[coef_matrix[,1]!=0,1]
    gene_pair_df <- data.frame(gene_pair=names(gene_pairs), coef=gene_pairs)
}

### GENE PAIRS LOAD ###

load_gene_pairs <- function(gene_pairs_path, gene_pairs_files, common_genes) {
    gene_pairs <- lapply(gene_pairs_files, function(gene_pairs_file) {
        print(gene_pairs_file)
        gene_pairs <- fread(file.path(gene_pairs_path, gene_pairs_file))$gene_pair
        gene_pairs_df <- t(sapply(gene_pairs, function(gene_pair) { unlist(strsplit(gene_pair, "_")[[1]]) }))
        colnames(gene_pairs_df) <- c("gene1", "gene2")
        common_genes_index <- apply(gene_pairs_df, 1, function(genes) { (genes[1] %in% common_genes) & (genes[2] %in% common_genes) })
#        gene_pairs_df <- gene_pairs_df[common_genes_index,]
        gene_pairs_df
    })
    names(gene_pairs) <- gene_pairs_files
    gene_pairs
}

### rpart ###

library(rpart)

derive_uncut_decision_tree <- function(rnk_matrix, sample_labels, gene_pairs_df, nclusters) {
    gpri_matrix <- cbind(generate_gene_pair_rnk_indicator_matrix(rnk_matrix, gene_pairs_df, nclusters), sample_labels)
    rpart(sample_labels ~ ., method="class", data=as.data.frame(gpri_matrix))
}

decision_tree_validation <- function(rpart_result, test_data, test_labels, gene_pairs_df, nclusters) {
    prediction_results <- do.call("rbind", lapply(names(test_data), function(tdata) {
        print(paste("Testing on", tdata))
        test_indicator_matrix <- as.data.frame(cbind(generate_gene_pair_rnk_indicator_matrix(test_data[[tdata]], gene_pairs_df, nclusters), test_labels[[tdata]]$label))
        colnames(test_indicator_matrix)[ncol(test_indicator_matrix)] <- "labels"
        tprediction <- predict(rpart_result, type="class", newdata=as.data.frame(test_indicator_matrix)) 
        predicted_intrpoor <- sum(tprediction=="intr_poor")
        predicted_good <- sum(tprediction=="good")
        correct_intrpoor <- sum(which(tprediction=="intr_poor") %in% which(test_indicator_matrix$labels=="intr_poor"))
        correct_good <-  sum(which(tprediction=="good") %in% which(test_indicator_matrix$labels=="good"))
        total_intrpoor <- sum(test_indicator_matrix$labels=="intr_poor")
        total_good <- sum(test_indicator_matrix$labels=="good")
        c(total_intrpoor, predicted_intrpoor, correct_intrpoor, total_good, predicted_good, correct_good)
    }))
    rownames(prediction_results) <- names(test_data)
    colnames(prediction_results) <- c("total_intrpoor", "predicted_intrpoor", "correct_intrpoor", "total_good", "predicted_good", "correct_good")
    prediction_results <- as.data.frame(prediction_results)
    data.frame(intrpoor_sensitivity=sum(prediction_results$correct_intrpoor)/sum(prediction_results$total_intrpoor), intrpoor_specificity=sum(prediction_results$correct_intrpoor)/sum(prediction_results$predicted_intrpoor), good_sensitivity=sum(prediction_results$correct_good)/sum(prediction_results$total_good), good_specificity=sum(prediction_results$correct_good)/sum(prediction_results$predicted_good))
}


### Rough functions


decision_tree_validation_generate_dataset_results <- function(rpart_result, test_data, test_labels, gene_pairs_df, nclusters) {
    prediction_results <- do.call("rbind", lapply(names(test_data), function(tdata) {
        print(paste("Testing on", tdata))
        test_indicator_matrix <- as.data.frame(cbind(generate_gene_pair_rnk_indicator_matrix(test_data[[tdata]], gene_pairs_df, nclusters), test_labels[[tdata]]$label))
        colnames(test_indicator_matrix)[ncol(test_indicator_matrix)] <- "labels"
        tprediction <- predict(rpart_result, type="class", newdata=as.data.frame(test_indicator_matrix)) 
        predicted_intrpoor <- sum(tprediction=="intr_poor")
        predicted_good <- sum(tprediction=="good")
        correct_intrpoor <- sum(which(tprediction=="intr_poor") %in% which(test_indicator_matrix$labels=="intr_poor"))
        correct_good <-  sum(which(tprediction=="good") %in% which(test_indicator_matrix$labels=="good"))
        total_intrpoor <- sum(test_indicator_matrix$labels=="intr_poor")
        total_good <- sum(test_indicator_matrix$labels=="good")
        c(total_intrpoor, predicted_intrpoor, correct_intrpoor, total_good, predicted_good, correct_good)
    }))
    rownames(prediction_results) <- names(test_data)
    colnames(prediction_results) <- c("total_intrpoor", "predicted_intrpoor", "correct_intrpoor", "total_good", "predicted_good", "correct_good")
    prediction_results <- as.data.frame(prediction_results)
#    data.frame(intrpoor_sensitivity=sum(prediction_results$correct_intrpoor)/sum(prediction_results$total_intrpoor), intrpoor_specificity=sum(prediction_results$correct_intrpoor)/sum(prediction_results$predicted_intrpoor), good_sensitivity=sum(prediction_results$correct_good)/sum(prediction_results$total_good), good_specificity=sum(prediction_results$correct_good)/sum(prediction_results$predicted_good))
}
