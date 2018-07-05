featureSelect <- function(infile, outdir, feature_count, use_all, no_probe_anno, sep) {
    # Load the data
    alldata <- loadData(inpath = infile, no_probe_anno = no_probe_anno, sepstr = sep)
    
    if (use_all != TRUE) {
        exe_without_all_data(alldata, feature_count)
    } else {
        # Use the entire data for training and feature selection and
        # No validation.
        exe_with_all_data(alldata, feature_count)

    }
}


exe_with_all_data <- function(alldata, feature_count) {

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
}

exe_without_all_data <- function(alldata, feature_count) {
	# Partition the data, use one half for feature selection and training
	# and other half for validation, generate a figure
	alldata2 <- doPartitionAlternate(alldata)
	cdata <- alldata2$cdata
	trainC <- alldata2$trainC
	features = getFeaturesReliefF(cdata = cdata, trainC = trainC, featureCount = feature_count)
    print("Features with all data:")
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


