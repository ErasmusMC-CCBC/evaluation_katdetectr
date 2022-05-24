# function for calculating performance metrics regarding labeling of variants
calculateConfusionVariants <- function(x, variantData){

    # get de original data on which the model made its prediction
    data <- variantData$genomicVariants[[base::unique(x$sampleNames)]]

    # if a kataegis foci is detected in this sample
    if(base::unique(x$detectedKataegisInSample)){
        # convert the the tibble to a granges
        gr <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
        # subset the original data by the kataegis foci in order to obtain the variants that were classified as kataegis variants
        detectedKataegisVariants <- IRanges::subsetByOverlaps(data, gr)
        # subset the original data by the kataegis foci in order to obtain the variants that were classified as NOT kataegis variants
        detectedNotKataegisVariants <- plyranges::filter_by_non_overlaps(data, gr)
        # calculate the number of true positives, false positives, true negatives and false negatives variants in the sample
        TP <- base::sum(detectedKataegisVariants$kataegis)
        FP <- base::sum(!detectedKataegisVariants$kataegis)
        TN <- base::sum(!detectedNotKataegisVariants$kataegis)
        FN <- base::sum(detectedNotKataegisVariants$kataegis)
        n <- TP + FP + TN + FN
    } else {
        # if no kataegis foci is detected in the sample
        TP <- 0
        FP <- 0
        TN <- base::sum(!data$kataegis)
        FN <- base::sum(data$kataegis)
        n <- TP + FP + TN + FN
    }

    confusionMatrix <- tibble::tibble(TP = TP, FP = FP, TN = TN, FN = FN, n = n, TMB = base::unique(x$TMB), TMBcat = base::unique(x$TMBcat), kataegisInSample = unique(x$kataegisInSample), detectedKataegisInSample = unique(x$detectedKataegisInSample))

    return(confusionMatrix)
}
