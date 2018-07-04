featureSelect <- function(infile, outdir, feature_count, use_all, sep) {
    # Load the data
    alldata <- loadData(inpath = infile, sepstr = sep)
    
    if (use_all != TRUE) {
        exe_without_all_data(alldata, feature_count)
    } else {
        # Use the entire data for training and feature selection and
        # No validation.
        exe_with_all_data(alldata, feature_count)

    }
}

getTrainC_with_all <- function(alldata) {
    lclass <- alldata$lclass
    trainC <- NA
    if ("I" %in% lclass) {
        trainC1 <- lclass
        trainC1[trainC1 == "I"] <- "R"
        trainC <- droplevels(trainC1)
    } else {
        trainC <- lclass
    }
    return (trainC)
    
}

exe_with_all_data <- function(alldata, feature_count) {

	trainC <- getTrainC_with_all(alldata)
	cdata <- alldata$cdata
	features = getFeaturesReliefF(cdata = cdata, trainC = trainC, featureCount = feature_count)
    print("Features with all data:")
    print("-----------------------")
    print(features)
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
	# get the expression for the training
	trainExp <- getExpData(cdata, features, trainC)
	modT <- do_train(trainExp, trainC)
	# Do the test
	testC <- alldata2$testC
	testExp <- getExpData(cdata, features, testC)
	prediction <- do_predict(testExp, modT) 
    print("Prediction")
    print("----------")
    print(prediction)

}

getExpData <- function(cdata, features, labels) {

	labels_names <- names(labels)
    expdata1 <- cdata[features, labels_names]
	expdata <- t(expdata1)
    return(expdata)
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


