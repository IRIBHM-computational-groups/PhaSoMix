suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(xgboost)
})
 
 
source("src/utils/utility_functions.r")
 
# =============================================================================
# FONCTIONS 
# =============================================================================

select_snps_by_variance <- function(genotype_matrix, top_n) {
  message("   -> Selecting top ", top_n, " SNPs by variance...")
  
  snp_variance <- apply(genotype_matrix, 2, var, na.rm = TRUE)
  
  top_snps <- order(snp_variance, decreasing = TRUE)[1:min(top_n, length(snp_variance))]
  
  selected_matrix <- genotype_matrix[, top_snps, drop = FALSE]
  
  message("   -> Selected ", ncol(selected_matrix), " SNPs")
  
  return(selected_matrix)
}

extract_genotype_matrix <- function(vcf_data, sample_list, top_n = 10000) {
 
    meta_cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    sample_cols <- setdiff(colnames(vcf_data), meta_cols)
 
    message("   -> Total samples in VCF: ", length(sample_cols))
    message("   -> Total variants in VCF: ", nrow(vcf_data))
 
    samples_to_keep <- intersect(sample_list, sample_cols)
 
 
    if (length(samples_to_keep) == 0) {
        stop("No matching samples found in VCF!")
    }
 
    message("   -> Samples to retain: ", length(samples_to_keep))
 
    gt_data <- vcf_data[, samples_to_keep]
    rownames(gt_data) = vcf_data$POS
 
    gt_matrix <- as.matrix(gt_data)
 
    gt_numeric <- t(matrix(as.numeric(gt_matrix), nrow = nrow(gt_matrix), ncol = ncol(gt_matrix)))
    rownames(gt_numeric) <- colnames(gt_matrix)
    colnames(gt_numeric) <- rownames(gt_matrix)
 
    na_count <- colSums(is.na(gt_numeric))
    if (any(na_count > 0)) {
        message("   -> Warning: ", sum(na_count > 0), " variants with missing data")
 
        # Supprimer les variants avec >10% de missing
        keep_vars <- na_count < (nrow(gt_numeric) * 0.1)
        gt_numeric <- gt_numeric[, keep_vars, drop = FALSE]
 
        # Imputer les NAs restants avec la valeur modale
        for (j in 1:ncol(gt_numeric)) {
            if (any(is.na(gt_numeric[, j]))) {
                mode_val <- as.numeric(names(sort(table(gt_numeric[, j]), decreasing = TRUE))[1])
                gt_numeric[is.na(gt_numeric[, j]), j] <- mode_val
            }
        }
    }
 
    # Supprimer les variants monomorphes (tous 0 ou tous 1)
    variant_sd <- apply(gt_numeric, 2, sd)
    polymorphic <- variant_sd > 0
    if (sum(!polymorphic) > 0) {
        message("   -> Removing ", sum(!polymorphic), " monomorphic variants")
        gt_numeric <- gt_numeric[, polymorphic, drop = FALSE]
    }
    
    if (ncol(gt_numeric) > top_n) gt_numeric = select_snps_by_variance(gt_numeric, top_n)
 
    message("   -> Final matrix: ", nrow(gt_numeric), " samples x ", ncol(gt_numeric), " variants")
  
  return(gt_numeric)
}

run_xgboost <- function(genotype_matrix, sample_info, output_path) {   
  
    message(">> Training XGBoost model ...")
    
    snps_used = colnames(genotype_matrix)
  
    if (!all(sample_info$category %in% c(0,1))) {
        stop("category column must contain only 0 (train) and 1 (test)")
    }
    
    sample_ids <- rownames(genotype_matrix)
    data_df <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
    data_df <- merge(data_df, sample_info, by = "sample", all.x = TRUE)

    train_samples <- data_df[data_df$category == 0, "sample"]
    test_samples  <- data_df[data_df$category == 1, "sample"]

    X_train <- genotype_matrix[train_samples, , drop = FALSE]
    y_train <- data_df[match(train_samples, data_df$sample), "ancestry"]
    
    classes       <- sort(unique(y_train))
    label_mapping <- setNames(seq_along(classes) - 1, classes)
    y_train_encoded <- label_mapping[y_train]
    num_classes   <- length(classes)

    message("   -> Training samples: ", length(train_samples))
    message("   -> Classes: ", paste(names(label_mapping), collapse = ", "))

    dtrain <- xgb.DMatrix(data = X_train, label = y_train_encoded)

    params <- list(
        objective        = "binary:logistic",
        eval_metric      = "logloss",
        max_depth        = 6,
        eta              = 0.3,
        subsample        = 0.8,
        colsample_bytree = 0.8,
        min_child_weight = 1
    )
    if (num_classes > 2) {
        params$objective   <- "multi:softprob"
        params$eval_metric <- "mlogloss"
        params$num_class   <- num_classes
    }

    # -----------------------------
    # Training (toujours)
    # -----------------------------
    message("   -> Training XGBoost...")

    if (length(test_samples) == 0) {
        cv <- xgb.cv(
            params = params, data = dtrain,
            nrounds = 100, nfold = 5,
            early_stopping_rounds = 10, verbose = 0
        )
        best_nrounds <- cv$best_iteration
    } else {
        X_test        <- genotype_matrix[test_samples, , drop = FALSE]
        y_test        <- data_df[match(test_samples, data_df$sample), "ancestry"]
        y_test_encoded <- label_mapping[y_test]
        dtest  <- xgb.DMatrix(data = X_test, label = y_test_encoded)
        best_nrounds <- 100   # early stopping via watchlist
    }

    watchlist <- if (length(test_samples) > 0) list(train = dtrain, test = dtest) else list(train = dtrain)

    model <- xgb.train(
        params    = params,
        data      = dtrain,
        nrounds   = best_nrounds,
        watchlist = watchlist,
        early_stopping_rounds = if (length(test_samples) > 0) 10 else NULL,
        verbose   = 0
    )

    if (length(test_samples) == 0) {
        xgb_object <- list(
            model         = model,
            snps_used     = snps_used,
            num_classes   = num_classes,
            label_mapping = label_mapping
        )
        output_model <- file.path(output_path, "xgb_model.rds")
        saveRDS(xgb_object, output_model)
        message("   -> Model saved: ", output_model)
    }

    # -----------------------------
    # Prédiction (seulement si test samples)
    # -----------------------------
    if (length(test_samples) == 0) return(output_model)

    message("   -> Prediction XGBoost of test samples: ", length(test_samples))

    pred_probs <- predict(model, dtest)

    if (num_classes == 2) {
        label_0 <- names(label_mapping)[which(label_mapping == 0)]
        label_1 <- names(label_mapping)[which(label_mapping == 1)]
        results <- data.frame(
            sample          = test_samples,
            true_label      = y_test,
            prob_class_1    = pred_probs,
            predicted_label = ifelse(pred_probs >= 0.5, label_1, label_0),
            stringsAsFactors = FALSE
        )
    } else {
        pred_probs_matrix <- matrix(pred_probs, ncol = num_classes, byrow = TRUE)
        colnames(pred_probs_matrix) <- names(label_mapping)
        predicted_label <- names(label_mapping)[max.col(pred_probs_matrix)]
        results <- cbind(
            data.frame(sample = test_samples, true_label = y_test,
                       predicted_label = predicted_label, stringsAsFactors = FALSE),
            as.data.frame(pred_probs_matrix)
        )
    }

    accuracy <- mean(results$true_label == results$predicted_label)
    message("   -> Test Accuracy: ", round(accuracy * 100, 2), "%")

    output_file <- file.path(output_path, "xgboost_prediction.txt")
    write.table(results, output_file, col.names = TRUE, row.names = FALSE,
                sep = "\t", quote = FALSE)
    message("   -> Results saved: ", output_file)

    return(list(results = results, accuracy = accuracy, output_file = output_file))
}



get_xgboost_prediction <- function(sample, input_vcf, chr, xgboost_model_path) {
    
    input_vcf = sub("chrCHR",chr,input_vcf)
    
    if (!file.exists(input_vcf)){
        stop("File doesn't exist: ", input_vcf)
    }
        
    cols = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample)
 
    vcf_cols <- fread(input_vcf, header = TRUE, skip = "#CHROM",
                      sep = "\t", nrow = 0, data.table = FALSE)
 
    if (! sample %in% colnames(vcf_cols)) return(NA)
    
    vcf = fread(input_vcf, header=T, skip = "#CHROM", sep = "\t", select = cols, data.table=F)
    
    model_path = file.path(xgboost_model_path,chr,"xgb_model.rds")
    xgb_object  <- readRDS(model_path)
    model       <- xgb_object$model
    snps_used   <- xgb_object$snps_used
    num_classes <- xgb_object$num_classes
    label_mapping <- xgb_object$label_mapping
    
    if (length(setdiff(snps_used,vcf$POS)) > 0){
        stop("No match between VCF SNP and model SNP")
    }
    
    gt_data <- vcf[match(snps_used,vcf$POS), sample]
    gt_data <- as.numeric(gt_data)
    dtest <- xgb.DMatrix(data = matrix(gt_data, nrow = 1))
    pred_probs <- predict(model, dtest)
    if (num_classes == 2) {
        predicted_label <- ifelse(pred_probs >= 0.5, xgb_object$label_1, xgb_object$label_0)
    } else {
        pred_probs_matrix <- matrix(pred_probs, ncol = num_classes, byrow = TRUE)
        colnames(pred_probs_matrix) <- names(label_mapping)
        predicted_label <- names(label_mapping)[max.col(pred_probs_matrix)]
    }
    return(predicted_label)
}


get_xgboost_prediction_prob <- function(sample, input_vcf, chr, xgboost_model_path) {
    
    input_vcf = sub("chrCHR",chr,input_vcf)
    
    if (!file.exists(input_vcf)){
        stop("File doesn't exist: ", input_vcf)
    }
        
    cols = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample)
 
    vcf_cols <- fread(input_vcf, header = TRUE, skip = "#CHROM",
                      sep = "\t", nrow = 0, data.table = FALSE)
 
    if (! sample %in% colnames(vcf_cols)) return(NA)
    
    vcf = fread(input_vcf, header=T, skip = "#CHROM", sep = "\t", select = cols, data.table=F)
    
    model_path = file.path(xgboost_model_path,chr,"xgb_model.rds")
    xgb_object  <- readRDS(model_path)
    model       <- xgb_object$model
    snps_used   <- xgb_object$snps_used
    num_classes <- xgb_object$num_classes
    label_mapping <- xgb_object$label_mapping
    
    if (length(setdiff(snps_used,vcf$POS)) > 0){
        stop("No match between VCF SNP and model SNP")
    }
    
    gt_data <- vcf[match(snps_used,vcf$POS), sample]
    gt_data <- as.numeric(gt_data)
    dtest <- xgb.DMatrix(data = matrix(gt_data, nrow = 1))
    pred_probs <- predict(model, dtest)
    if (num_classes == 2) {
        predicted_label <- ifelse(pred_probs >= 0.5, xgb_object$label_1, xgb_object$label_0)
    } else {
        pred_probs_matrix <- matrix(pred_probs, ncol = num_classes, byrow = TRUE)
        colnames(pred_probs_matrix) <- names(label_mapping)
        predicted_label <- names(label_mapping)[max.col(pred_probs_matrix)]
    }
    return(pred_probs_matrix)
}


# =============================================================================
# FONCTION PRINCIPALE : Pipeline complet
# =============================================================================

run_xgboost_pipeline <- function(vcf, sample_info, output_path, chr_label) {
  
    message("Processing ", chr_label)
 
    genotype_matrix <- extract_genotype_matrix(vcf, sample_info$sample)
 
    pred_xgboost = run_xgboost(genotype_matrix, sample_info, output_path)
 
    return(pred_xgboost)
}