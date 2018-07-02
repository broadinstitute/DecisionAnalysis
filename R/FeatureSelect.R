doPartitionAlternate <- function(alldata) {
    lclass <- alldata$lclass
    odd <- seq(1, length(lclass), 2)
    even <- seq(2, length(lclass), 2)
    trainC <- lclass[odd]
    testC <- lclass[even]
    alldata2 <- c(alldata, list(trainC = trainC, testC = testC))
    return (alldata2)
}

getFeaturesReliefF <- function(alldata, ltimes = 5, featureCount = 5, lEst = "ReliefFexpRank") {

    trainC <- alldata$trainC
    cdata <- alldata$cdata
    
    seed1 <- 100
    set.seed(seed1)
    
    finalCVCount <- getCVCount(trainC)
    index <- createMultiFolds(trainC, k = finalCVCount, times = ltimes)

    all_list <- c()
    for (j in 1:length(index)) {
        train_index <- index[[j]]
        ydata <- trainC[train_index]
        xdata1 <- cdata[, names(ydata)]
        xdata2 <- t(xdata1)
        tdata <- data.frame(xdata2, ydata)
        estReliefF <- attrEval(ydata ~ ., tdata, estimator= lEst)
        all_list <- c(all_list, estReliefF)
        
    }

    all_list_df <- data.frame(names(all_list), as.numeric(all_list))
    colnames(all_list_df) <- c("Variables", "Values")
    finalList1 <- ddply(all_list_df,  .(Variables), function(x) mean(x$Values))
    names(finalList1)[2] <- "Values"
    finalList <- finalList1[order(finalList1$Values, decreasing = TRUE), ]
    final_genes <- as.character(finalList$Variables)
    
    features <- final_genes[1:featureCount]
    alldata2 <- c(alldata, list(features = features))   

}

validation_F <- function (alldata, featureCount = 5) {

    # Set the seed for multifold creation and random forest
    seed2 <- 200
    set.seed(seed2)

    trainC <- alldata$trainC
    finalCVCount <- getCVCount(trainC)
    indexT <- createMultiFolds(trainC, k = finalCVCount, times = 5)

    ctrlT <- trainControl(method = "repeatedcv", number = finalCVCount,
        repeats = 5, returnResamp = "all", savePredictions = "all",
        classProbs = TRUE, index = indexT)


    features <- alldata$features[1:featureCount]
    cdata <- alldata$cdata
    xdata1 <- cdata[, names(trainC)]
    xdata2 <- t(xdata1)
    print("features")
    print(features)
    print("colnames start")
    print(colnames(xdata2))
    print("colnames end")
    trainData <- data.frame(xdata2[, features], trainC)
    modT <- train( trainC ~ ., data = trainData, method = "rf", 
            trControl = ctrlT)

    testC <-  alldata$testC
    ydata1 <- cdata[, names(testC)]
    ydata2 <- t(ydata1)
    testData <- data.frame(ydata2[, features])
    resClasses <- predict(modT, newdata = testData)
    resProbs <- predict(modT, newdata = testData, type = "prob")
    lmat <- confusionMatrix(resClasses, testC)
    accuracy <- as.numeric(lmat$overall["Accuracy"])

    alldata2 <- c(alldata, list(modT = modT, lmat = lmat, accuracy = accuracy, 
                resClasses = resClasses, resProbs = resProbs))
    
    return (alldata2)

}

featureSelect <- function(alldata, featureCount) {

    alldata2 <- doPartitionAlternate(alldata)
    alldata3 <- getFeaturesReliefF(alldata2, featureCount)

    plotCombinedMetric_F(alldata3, lplotname, ltitle, outpath, lmethod, featureCount)
    

}
