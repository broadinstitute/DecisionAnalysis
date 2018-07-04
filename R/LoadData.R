loadData <- function(inpath, sepstr) {

    expdata1 <- read.csv(inpath, header = FALSE, sep = sepstr, row.names = 1, stringsAsFactors = FALSE)

    MIC0 <- expdata1[1,]
    lclass0 <- expdata1[2,]
    strain0 <- expdata1[3,]
    expdata2 <- expdata1[-c(1:3),]

    lclass1 <- as.character(lclass0)
    MIC1 <- as.numeric(MIC0)
    strain1 <- as.character(strain0)
    strain2 <-gsub("_foldInd$", "", strain1) 
    colnames(expdata2) <- strain2

    names(lclass1) <- strain2

    lclass2 <- factor(lclass1)
    names(MIC1) <- strain2
    MIC <- sort(MIC1)
    lclass <- lclass2[names(MIC)] 

    cdata <- data.matrix(expdata2[, names(lclass)])

    probe_names <- rownames(cdata)

    fVals <- lapply(probe_names, getAnnoU)
    fMap <- unlist(fVals)
    names(fMap) <- probe_names

    alldata <- list(lclass = lclass, MIC = MIC, cdata = cdata, fMap = fMap)

    return (alldata)
}

# The annotation of a feature (the species is not provided by the user and hence not checked)
getAnnoU <- function(feature) {
    
    res <- NA
    fval <- entrez_search(db="protein", feature)
    lids <- fval$ids
    if (length(lids) > 0) {
    
        lid <- lids[1]
        fval2 <- entrez_summary(db="protein", id = lid)
        ltitle <- fval2$title
        res1 <- gsub("(.*?)[[:space:]]*\\[.*$", "\\1", fval2$title)
        res <- paste0(feature, ":", res1)
    } else {
        res <- paste0(feature, ":")
    }
    return (res)
}


# The annotation for a feature (checked version), this is a full version when the species is 
# provided by the user and verified.
getAnnoC <- function(feature, species) {
    
    fval <- entrez_search(db="protein", feature)
    for (lid in fval$ids) {
        fval2 <- entrez_summary(db="protein", id = lid)
        ltitle <- fval2$title
        if (grepl(species, ltitle)) {
            res <- gsub("(.*?)[[:space:]]*\\[.*$", "\\1", fval2$title)
            return (paste0(feature_full, ":", res))
        }
    }
    return (paste0(feature, ":"))
}

