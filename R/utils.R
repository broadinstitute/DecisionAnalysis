doPartitionAlternate <- function(alldata) {
    # If there are intermediate data, we shall not train on
    # them, rather we shall use them for test and shall 
    # consider them resistant during test.

    lclass <- alldata$lclass

    lclass_next <- NA
    inter_names <- NA    
    if ("I" %in% lclass) {
        inter1 <- subset(lclass, "I" == lclass)
        inter <- droplevels(inter)
        lclass_names <- names(lclass)
        inter_names <- names(inter)
        lclass_no_inter_names <- setdiff(lclass_names, inter_names)
        lclass_no_inter1 <- lclass[lclass_no_inter_names]
        lclass_next <- droplevels(lclass_no_inter1)
    } else {
        lclass_next <- lclass
    }
    
    testC_final <- NA
    if (!is.na(inter_names)) {
        # Put the inter_names into the test set
        testC_names1 <- names(testC)
        testC_vals1 <- as.character(testC)
        testC_names <- c(testC_names1, inter_names)
        inter_vals <- rep("R", length(inter_names))
        testC_vals <- c(testC_vals1, inter_vals)
        names(testC_vals) <- testC_names
        testC_final <- factor(testC_vals)        
    } else {
        testC_final <- testC
    }

    odd <- seq(1, length(lclass_next), 2)
    even <- seq(2, length(lclass_next), 2)
    trainC <- lclass_next[odd]
    testC <- lclass_next[even]

    alldata2 <- c(alldata, list(trainC = trainC, testC = testC_final))
    return (alldata2)
}

getCVCount <- function(classLabels) {

    standardCVCount <- 5
    smallClassCount <- min(table(classLabels))
    finalCVCount <- min(standardCVCount, smallClassCount)
    return (finalCVCount)
}


getFeaturesReliefF <- function(cdata, trainC, ltimes = 5, featureCount = 5, lEst = "ReliefFexpRank") {

    #trainC <- alldata$trainC
    #cdata <- alldata$cdata
    
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
    return (features)
    #alldata2 <- c(alldata, list(features = features))   

}

do_train <- function(trainExp, trainC) {

    seed2 <- 200
    set.seed(seed2)

    finalCVCount <- getCVCount(trainC)
    indexT <- createMultiFolds(trainC, k = finalCVCount, times = 5)

    ctrlT <- trainControl(method = "repeatedcv", number = finalCVCount,
        repeats = 5, returnResamp = "all", savePredictions = "all",
        classProbs = TRUE, index = indexT)

    trainData <- data.frame(trainExp, trainC)
    modT <- train( trainC ~ ., data = trainData, method = "rf", 
            trControl = ctrlT)

    return (modT)
}




do_predict <- function(testData, modT) {

    resClasses <- predict(modT, newdata = testData)
    resProbs <- predict(modT, newdata = testData, type = "prob")

    res <- list()
    res$"resClasses" <- resClasses
    res$"resProbs" <- resProbs

    return (res)
}


