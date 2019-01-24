validate <- function(train_file, validate_file, outdir, prefix, feature_count, sep) {
    dir.create(outdir, recursive = TRUE)

    t_data <-  loadData(inpath = train_file, no_probe_anno = TRUE, sepstr = sep) 
    validate_data <- loadValidateData(inpath = validate_file, sep = sep)

    # get the appropriate number of features
    features <- getFeatures(validate_data, feature_count)

    print("features")
    print("....")
    print(features)
    print("....")

    # do_train
    trainC <- getTrainC_with_all(t_data)
    trainData <- t_data$cdata
    trainExp <- getExpData(trainData, features, trainC)
    print("....")
    print("trainC")
    print(trainC)
    print("....")
    modT <- do_train(trainExp, trainC) 
    print("modT")
    print("....")
    print(modT)
    print("....")

    testExp <- get_exp_for_validate(validate_data, features)
    prediction <- do_predict(testExp, modT)

    print("prediction")
    print("....")
    print(prediction)
    print("....")

    pred_file <- paste0(outdir, "/", prefix, "_pred_table.txt")
    sink(pred_file)
    print(prediction)
    sink()
    print("end of prediction")
    
}

validate_BS <- function(train_file, validate_file, outdir, prefix, feature_count, sep) {

    seed1 <- 1000
    set.seed(seed1)

    dir.create(outdir, recursive = TRUE)

    t_data <-  loadData(inpath = train_file, no_probe_anno = TRUE, sepstr = sep) 
    validate_data <- loadValidateData(inpath = validate_file, sep = sep)

    # get the appropriate number of features
    features <- getFeatures(validate_data, feature_count)

    print("features")
    print("....")
    print(features)
    print("....")

    # do_train
    print(".....")
    print("t_data")
    print("....")
    print(t_data)
    print("....")

    trainC <- getTrainC_with_all(t_data)
    trainData <- t_data$cdata
    trainExp <- getExpData(trainData, features, trainC)

    print("trainExp")
    print("....")
    print(trainExp)
    print("....")

    print("....")
    print("trainC")
    print(trainC)
    print("....")

    
    trainC_names <- names(trainC)
    trainC_sus <- subset(trainC, trainC == "S")
    trainC_res <- subset(trainC, trainC == "R")
    trainC_sus_names <- names(trainC_sus)
    trainC_res_names <- names(trainC_res)

    print("trainC_names")
    print("....")
    print(trainC_names)
    print("....")

    trainC_sus_lst <- createResample(trainC_sus_names, times = 10)
    trainC_res_lst <- createResample(trainC_res_names, times = 10)

    single_data_lim <- 10
    resample_count <- single_data_lim ^2
    validate_sample_counts <- dim(validate_data)[2]
    m <- matrix(0, ncol = validate_sample_counts, nrow = resample_count)
    pred_df <- data.frame(m)
    colnames(pred_df) <- colnames(validate_data) 

    arr1 <- outer(1:single_data_lim, 1:single_data_lim, FUN = "paste", sep = "_")
    resample_ids <- t(arr1)
    dim(resample_ids) <- NULL
    rownames(pred_df) <- resample_ids


    testExp <- get_exp_for_validate(validate_data, features)

    for (sus_ind in 1:single_data_lim) {
        for (res_ind in 1:single_data_lim) {
            sus_loc <- trainC_sus_lst[[sus_ind]]
            res_loc <- trainC_res_lst[[res_ind]]
            sus_samples <- trainC_sus_names[sus_loc]
            res_samples <- trainC_res_names[res_loc]
            ltrainC_names <- c(sus_samples, res_samples)
            ltrainC <- trainC[ltrainC_names]
            ltrainExp <- trainExp[ltrainC_names,]

            lmodT <- do_train(ltrainExp, ltrainC)
            lprediction <- do_predict(testExp, lmodT)
            lindex <- paste0(sus_ind, "_", res_ind)
            lpred_rnames <- rownames(lprediction[[2]])
            pred_df[lindex, lpred_rnames] <-lprediction[[2]]$R 

            print(paste0("Validation done for resample id: ", lindex))
        }
    }

    R_call <- apply(pred_df, 2, function(x) {sum(x >= 0.5)})
    S_call <- apply(pred_df, 2, function(x) {sum(x < 0.5)})
    mean_val <- apply(pred_df, 2, mean)
    sd_val <- apply(pred_df, 2, sd)
    se_val <- sd_val /sqrt(validate_sample_counts)

    bootstraping_df <- data.frame(mean_val, sd_val, se_val, R_call, S_call)

    print("Bootstrapping samples")
    print("....")
    print(bootstraping_df)
    print("....")

    BS_raw_file <- paste0(outdir, "/", prefix, "_BS_raw_table.txt")
    write.table(pred_df, BS_raw_file, sep = "\t", quote = FALSE)

    BS_metrics_file <- paste0(outdir, "/", prefix, "_BS_metrics_table.txt")
    write.table(bootstraping_df, BS_metrics_file, sep = "\t", quote = FALSE)
    print("end of Bootstrapping")


}

get_exp_for_validate <- function(vdata, features) {
    lExp1 <- vdata[features, ]
    lExp <- t(lExp1)
    return (lExp)
}

getFeatures <- function(cdata, feature_count) {

    features1 <- rownames(cdata)
    flen1 <- length(features1)
    flen <- NA
    if (feature_count == -1) {
        # default
        flen <- flen1
    } else if (feature_count == 0) {
        flen <- flen1
        print("WARN: featuren_count is 0.")
    } else if (feature_count > flen1) {
        flen <- flen1
        print(paste0("WARN: feature_count is: ", feature_count))
    } else {
        flen <- feature_count
    }

    features <- features1[1:flen]
    return (features)

}

