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

    # The save the model to the outfile
    model_file <- paste0(outdir, "/", prefix, "_model.rds")
    saveRDS(modT, model_file)

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
    testExp <- get_exp_for_validate(validate_data, features)

    print("trainExp")
    print("....")
    print(trainExp)
    print("....")

    print("....")
    print("trainC")
    print(trainC)
    print("....")

    print("Starting bootstrap")
    exe_bootstrap(trainExp, trainC, testExp, outdir, prefix)

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

