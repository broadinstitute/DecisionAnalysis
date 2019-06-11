featureSelect <- function(infile, outdir, prefix, feature_count, use_all, no_probe_anno, sep) {
     dir.create(outdir, recursive = TRUE)
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

featureSelect_BS <- function(infile, outdir, prefix, feature_count, no_probe_anno, sep) {
     dir.create(outdir, recursive = TRUE)
    # Load the data
    alldata <- loadData(inpath = infile, no_probe_anno = no_probe_anno, sepstr = sep)
    
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

	trainExp <- getExpData(cdata, features, trainC)
	testC <- alldata2$testC
	testExp <- getExpData(cdata, features, testC)

    exe_bootstrap(trainExp, trainC, testExp, outdir, prefix)

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

Leave_One_Out_CV <- function(infile, outdir, prefix, no_probe_anno = TRUE, sep = ',') {
    dir.create(outdir, recursive = TRUE)
    alldata <- loadData(inpath = infile, no_probe_anno = no_probe_anno, sepstr = sep)    

    # Tran and test in a loop
    lsamples <- colnames(alldata$cdata)
    sample_count <- length(lsamples)
    test_sample_vec <- c()
    resClasses_vec <- c()
    resProbs_mat <- c()
    for (lsample in lsamples) {
        #print(lsample)
        # Get training samples
        test_S <- lsample
        train_S <- setdiff(lsamples, test_S)
        
        print("train_S")
        print(train_S)
        print("test_S")
        print(test_S)
        print("-------------")

        # Now get accuracy based on training and test for each iteration
        lclass <- alldata$lclass
        ltrain_C <- factor(lclass[train_S])
        ltest_C <- factor(lclass[test_S])
        cdata <- alldata$cdata
        features <- rownames(cdata)
        trainExp <- getExpData(cdata, features, ltrain_C)
        testExp <- getExpData(cdata, features, ltest_C)

        modT <- do_train(trainExp, ltrain_C)
        # Now do the test
        lprediction <- do_predict(testExp, modT)
        resClasses <- factor(lprediction$resClasses)
        resProbs <- lprediction$resProbs
        test_sample_vec <- c(test_sample_vec, lsample)
        resClasses_vec <- c(resClasses_vec, as.character(resClasses))
        resProbs_mat <- rbind(resProbs_mat, resProbs)
    }

    resClasses_vec <- factor(resClasses_vec)
    resProbs_df <- data.frame(resProbs_mat)
    rownames(resProbs_df) <- test_sample_vec
    names(resClasses_vec) <- test_sample_vec
    prediction <- list()
    prediction$"resClasses" <- factor(resClasses_vec)
    prediction$"resProbs" <- resProbs_df

    resClasses <- prediction$resClasses
    testC <- alldata$lclass[test_sample_vec]
    lcmat <- confusionMatrix(resClasses, testC)
    
    alldata$testC <- testC
    prediction$"confusionMatrix" <- lcmat

    pred_file <- paste0(outdir, "/", prefix, "_pred.txt")
    sink(pred_file)
    print(prediction)
    sink()
    plot_file <- paste0(outdir, "/", prefix, "_plot.pdf")
    print(plot_file)
    draw_plot(alldata, prediction, features, prefix, plot_file)

    
}

Leave_Two_Out_CV <- function(infile, outdir, prefix, no_probe_anno = TRUE, sep = ',') {
    dir.create(outdir, recursive = TRUE)
    alldata <- loadData(inpath = infile, no_probe_anno = no_probe_anno, sepstr = sep)

    # Get number of susceptible samples and number of resistant samples
    lclass <- alldata$lclass
    # get the S_samples and R_samples
    all_samples <- names(lclass)
    S_samples <- names(lclass[lclass == "S"])
    R_samples <- setdiff(all_samples, S_samples)

    all_count <- length(lclass)
    S_count <- length(S_samples)
    R_count <- length(R_samples)
    
    # Go over all the susceptible and all the resistant samples to 
    # create (S -1) * (R -1) training sets; this would be a leave two out 
    # cross validation

    test_sample_vec <- c()
    resClasses_vec <- c()
    resProbs_mat <- c()

    for (S_lsam in S_samples) {
        for (R_lsam in R_samples) {

            # Build up the training and test set
            test_S <- c(S_lsam, R_lsam)
            train_S <- setdiff(all_samples, test_S)
            ltest_C <- lclass[test_S]
            ltrain_C <- lclass[train_S]

            print("train_S")
            print(train_S)
            print("test_S")
            print(test_S)
            print("-------------")

            cdata <- alldata$cdata
            features <- rownames(cdata)
            trainExp <- getExpData(cdata, features, ltrain_C)
            testExp <- getExpData(cdata, features, ltest_C)
            modT <- do_train(trainExp, ltrain_C)
            # Now do the test
            lprediction <- do_predict(testExp, modT)
            resClasses <- factor(lprediction$resClasses)
            resProbs <- lprediction$resProbs
    

            # How to store the data
            test_sample_vec <- c(test_sample_vec, test_S)
            resClasses_vec <- c(resClasses_vec, as.character(resClasses))
            resProbs_mat <- rbind(resProbs_mat, resProbs, make.row.names = FALSE)

        }
    }
    resProbs_mat_f <- cbind(test_samples = test_sample_vec, resProbs_mat)    
    resProbs_mat_f$test_samples<- factor(resProbs_mat_f$test_samples, levels = names(alldata$MIC))

    resClasses <- resClasses_vec
    names(resClasses) <- test_sample_vec
    resClasses <- factor(resClasses)

    testC <- lclass[test_sample_vec]
    lcmat <- confusionMatrix(resClasses, testC)

    prediction <- list()
    prediction$"resClasses" <- resClasses
    prediction$"resProbs" <- resProbs_mat
    
    alldata$testC <- testC
    prediction$"confusionMatrix" <- lcmat
    prediction$"resProbs_mat_f" <- resProbs_mat_f

    pred_file <- paste0(outdir, "/", prefix, "_pred.txt")
    sink(pred_file)
    print(prediction)
    sink()

    #p <- ggplot(resProbs_mat_f, aes(test_sample_vec, R))
    #p + geom_boxplot()

    point_plot_file <- paste0(outdir, "/", prefix, "_point_plot.pdf")
    print(point_plot_file)
    box_plot_file <- paste0(outdir, "/", prefix, "_box_plot.pdf")
    print(box_plot_file)
    draw_plot_LTOCV(alldata, prediction, features, prefix, point_plot_file,
        box_plot_file)

}

draw_plot_LTOCV <- function(alldata, prediction, features, prefix, point_plot_file, box_plot_file) {
    con_mat <- prediction$confusionMatrix 
    accuracy <- as.numeric(con_mat$overall["Accuracy"])
    resProbs <- prediction$resProbs
    resProbs_mat_f <- prediction$resProbs_mat_f
    
    predSample <- resProbs_mat_f$test_samples
    labels <- alldata$testC[predSample]
    predSample_u <- names(alldata$MIC)
    predMIC_u <- as.numeric(alldata$MIC[predSample_u])
    predMIC <- as.numeric(alldata$MIC[predSample])
    lgroups1 <- paste0(predSample, "_", predMIC)
    lgroups1_u <- paste0(predSample_u, "_", predMIC_u)
    lgroups <- factor(lgroups1, levels = lgroups1_u)

    mid_point <- max(alldata$MIC[alldata$lclass == 'S'])
    facet_var1 <- ifelse (predMIC <=mid_point , c('Sus'), c('Res'))
    facet_var <- factor(facet_var1, levels = c("Sus", "Res"))

    probRes <- data.frame(lnames = predSample, probs = resProbs[,1], 
			types = as.character(labels), lgroups = lgroups, 
			facet_var = facet_var)

    Palette1 <- c('forestgreen', 'red')

    fMap <- alldata$fMap
    ltitle <- get_title(features, accuracy, fMap, prefix)
    resProbs_mat_f <- prediction$resProbs_mat_f

    plt <- ggplot(probRes, aes(x = lgroups, y = probs, colour = facet_var)) +
        geom_point(size = 2) + 
        facet_grid(. ~ facet_var, scales = "free", space = "free") + 
        scale_colour_manual(values=Palette1) + 
        theme(axis.text.x = element_text(size=10,angle= 45))  + 
        xlab("Strain_MIC") + ylab("Probablity of resistance") +
		labs(colour = 'Strain\ngroups') + ggtitle(ltitle) 

    ggsave(point_plot_file)

    plt <- ggplot(probRes, aes(x = lgroups, y = probs, colour = facet_var)) +
        geom_point(size = 2) + geom_boxplot() + 
        facet_grid(. ~ facet_var, scales = "free", space = "free") + 
        scale_colour_manual(values=Palette1) + 
        theme(axis.text.x = element_text(size=10,angle= 45))  + 
        xlab("Strain_MIC") + ylab("Probablity of resistance") +
		labs(colour='Strain\ngroups') + ggtitle(ltitle) 

    ggsave(box_plot_file)
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
    print("mtry")
    print("....")
    print(modT$finalModel$mtry)
    print("....")
    print("importance")
    print("....")
    print(modT$finalModel$importance)
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

    plot_file <- paste0(outdir, "/", prefix, "_plot.pdf") 
    print(plot_file)
    draw_plot(alldata2, prediction, features, prefix, plot_file)

}

draw_plot <- function(alldata, prediction, features, prefix, plot_file) {
    con_mat <- prediction$confusionMatrix 
    accuracy <- as.numeric(con_mat$overall["Accuracy"])
    resProbs <- prediction$resProbs
    
    labels <- alldata$testC[rownames(resProbs)]
    predSample <- names(alldata$testC)
    predMIC <- as.numeric(alldata$MIC[predSample])
    lgroups1 <- paste0(predSample, "_", predMIC)
    lgroups <- factor(lgroups1, levels = lgroups1)

    mid_point <- max(alldata$MIC[alldata$lclass == 'S'])
    facet_var1 <- ifelse (predMIC <=mid_point , c('Sus'), c('Res'))
    facet_var <- factor(facet_var1, levels = c("Sus", "Res"))

    probRes <- data.frame(lnames = rownames(resProbs), probs = resProbs[,1], 
			types = as.character(labels), lgroups = lgroups, 
			facet_var = facet_var)

    Palette1 <- c('forestgreen', 'red')

    fMap <- alldata$fMap
    ltitle <- get_title(features, accuracy, fMap, prefix)

    plt <- ggplot(probRes, aes(x = lgroups, y = probs, colour = facet_var)) +
        geom_point(size = 3) +
        facet_grid(. ~ facet_var, scales = "free", space = "free") + 
        scale_colour_manual(values=Palette1) + 
        theme(axis.text.x = element_text(size=10,angle= 45))  + 
        xlab("Strain_MIC") + ylab("Probablity of resistance") +
		labs(colour='Strain\ngroups') + ggtitle(ltitle) 

    ggsave(plot_file)
}

get_title <- function(features, accuracy, fMap, prefix) {

    
    fVals1 <- fMap[features]
    fVals2 <- substr(fVals1, 1, 50)
    fVals <- paste(as.character(fVals2), collapse = "\n")

	ltitle <- paste0("Confidence of resistance\n", prefix, "\nAccuracy = ", accuracy, "\n", fVals)
    return (ltitle)
}

