featureSelect <- function(infile, outdir, prefix, feature_count, use_all, no_probe_anno, sep) {
    # Load the data
    alldata <- loadData(inpath = infile, no_probe_anno = no_probe_anno, sepstr = sep)
    
    if (use_all != TRUE) {
        exe_without_all_data(alldata, feature_count, outdir, prefix)
    } else {
        # Use the entire data for training and feature selection and
        # No validation.
        exe_with_all_data(alldata, feature_count, outdir, prefix)

    }
}


exe_with_all_data <- function(alldata, feature_count, outdir, prefix) {

	trainC <- getTrainC_with_all(alldata)
	cdata <- alldata$cdata
	features = getFeaturesReliefF(cdata = cdata, trainC = trainC, featureCount = feature_count)
    print("Features with all data:")
    print("-----------------------")
    print(features)
    print("....")
    print("trainC")
    print("....")
    print(trainC)
    print("....")

    features_file <- paste0(outdir, "/", prefix, "_features.txt")
    features_df <- as.data.frame(features)
    write.table(features_df, file = features_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
}

exe_without_all_data <- function(alldata, feature_count, outdir, prefix) {
	# Partition the data, use one half for feature selection and training
	# and other half for validation, generate a figure
	alldata2 <- doPartitionAlternate(alldata)
	cdata <- alldata2$cdata
	trainC <- alldata2$trainC
	features = getFeaturesReliefF(cdata = cdata, trainC = trainC, featureCount = feature_count)
    print("Features without all data:")
    print("-----------------------")
    print(features)
    print("....")
	# get the expression for the training
	trainExp <- getExpData(cdata, features, trainC)
	modT <- do_train(trainExp, trainC)
	# Do the test
	testC <- alldata2$testC

    print("trainC")
    print("....")
    print(trainC)
    print("....")
    print("testC")
    print("....")
    print(testC)
    print("....")

	testExp <- getExpData(cdata, features, testC)
	prediction <- do_predict(testExp, modT) 
    print("Prediction")
    print("----------")
    print(prediction)

    resClasses <- prediction$resClasses
    lcmat <- confusionMatrix(resClasses, testC)
    print("Confusion Matrix")
    print("....")
    print(lcmat)
    print("....")

    
    features_file <- paste0(outdir, "/", prefix, "_features.txt")
    features_df <- as.data.frame(features)
    write.table(features_df, file = features_file, col.names = FALSE, row.names = FALSE, quote = FALSE)


    prediction$"confusionMatrix" <- lcmat
    pred_file <- paste0(outdir, "/", prefix, "_pred.txt")
    sink(pred_file)
    print(prediction)
    sink()

}


doPartitionAlternate <- function(alldata) {
    lclass <- alldata$lclass
    odd <- seq(1, length(lclass), 2)
    even <- seq(2, length(lclass), 2)
    trainC <- lclass[odd]
    testC <- lclass[even]
    alldata2 <- c(alldata, list(trainC = trainC, testC = testC))
    return (alldata2)
}


