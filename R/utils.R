doPartitionAlternate <- function(alldata) {
    # If there are intermediate data, we shall not train on
    # them, rather we shall use them for test and shall 
    # consider them resistant during test.

    lclass <- alldata$lclass

    lclass_next <- NA
    inter_names <- c() 
    if ("I" %in% lclass) {
        inter1 <- subset(lclass, "I" == lclass)
        inter <- droplevels(inter1)
        lclass_names <- names(lclass)
        inter_names <- names(inter)
        lclass_no_inter_names <- setdiff(lclass_names, inter_names)
        lclass_no_inter1 <- lclass[lclass_no_inter_names]
        lclass_next <- droplevels(lclass_no_inter1)
    } else {
        lclass_next <- lclass
    }

    odd <- seq(1, length(lclass_next), 2)
    even <- seq(2, length(lclass_next), 2)
    trainC <- lclass_next[odd]
    testC <- lclass_next[even]
   
    testC_final <- NA
    if (length(inter_names) > 0) {
        # Put the inter_names into the test set
        testC_names1 <- names(testC)
        testC_vals1 <- as.character(testC)
        testC_names <- c(testC_names1, inter_names)
        inter_vals <- rep("R", length(inter_names))
        testC_vals <- c(testC_vals1, inter_vals)
        names(testC_vals) <- testC_names

        # One has to sort the testC_final with respect to MIC
        MIC <- alldata$MIC
        lMIC <- MIC[testC_names]
        lMIC_o <- order(lMIC)
        testC_final1 <- testC_vals[lMIC_o] 

        testC_final <- factor(testC_final1)  
    } else {
        testC_final <- testC
    }

    alldata2 <- c(alldata, list(trainC = trainC, testC = testC_final))

    return (alldata2)
}

getCVCount <- function(classLabels) {

    standardCVCount <- 5
    smallClassCount <- min(table(classLabels))
    finalCVCount <- min(standardCVCount, smallClassCount)
    return (finalCVCount)
}


getFeaturesReliefF <- function(cdata, trainC, featureCount = 10, ltimes = 5, lEst = "ReliefFexpRank") {

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

getExpData <- function(cdata, features, labels) {

    labels_names <- names(labels)
    expdata1 <- cdata[features, labels_names]
    expdata <- t(expdata1)
    return(expdata)
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


exe_bootstrap <- function(trainExp, trainC, testExp, outdir, prefix) {

    seed1 <- 1000
    set.seed(seed1)


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
    validate_sample_counts <- dim(testExp)[1]
    print("validate_sample_counts")
    print(validate_sample_counts)
    print("....")
    m <- matrix(0, ncol = validate_sample_counts, nrow = resample_count)
    pred_df <- data.frame(m)
    colnames(pred_df) <- rownames(testExp)

    arr1 <- outer(1:single_data_lim, 1:single_data_lim, FUN = "paste", sep = "_")
    resample_ids <- t(arr1)
    dim(resample_ids) <- NULL
    rownames(pred_df) <- resample_ids



    for (sus_ind in 1:single_data_lim) {
        for (res_ind in 1:single_data_lim) {
            lindex <- paste0(sus_ind, "_", res_ind)
            sus_loc <- trainC_sus_lst[[sus_ind]]
            res_loc <- trainC_res_lst[[res_ind]]
            
            sus_samples <- trainC_sus_names[sus_loc]
            res_samples <- trainC_res_names[res_loc]
            ltrainC_names <- c(sus_samples, res_samples)
            ltrainC <- trainC[ltrainC_names]
            ltrainExp <- trainExp[ltrainC_names,]
            lmodT <- do_train(ltrainExp, ltrainC)
            lprediction <- do_predict(testExp, lmodT)
            lpred_rnames <- rownames(lprediction[[2]])
            pred_df[lindex, lpred_rnames] <-lprediction[[2]]$R
            print(paste0("Validation starts for resample id: ", lindex))

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




