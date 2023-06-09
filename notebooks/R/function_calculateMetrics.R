# Calculates performance metrics.
# x should be a dataframe which contains columns that specify the number of true positives, false positives, true negatives and false negatives i.e. all info of a confusion matrix
calculateMetrics <- function(x){

    metrics <- tibble::tibble(
        TP = sum(x$TP),
        FP = sum(x$FP),
        TN = sum(x$TN),
        FN = sum(x$FN),

        accuracy = (TP + TN) / (TP + FP + TN + FN),
        # if the confusion matrix contains no true positives and no false negatives set TPR to 1
        TPR = dplyr::if_else(TP + FN != 0, TP / (TP + FN), 1),
        # if the confucion matrix contains no true negatives and no false positives set TNR to 1
        TNR = dplyr::if_else(TN + FP != 0, TN / (TN + FP), 1),
        # if the confucion matrix contains no false positives and no false negatives set F1 to 1
        F1 =  dplyr::if_else(FP + FN != 0 ,TP / (TP + 0.5 * (FP + FN)), 1),
        # nMCC = 0 when a column or row of the confusion matrix
        nMCC = (mltools::mcc(TP = TP, FP = FP, TN = TN, FN = FN) + 1) / 2,
        meanRuntime = mean(x$runTime, na.rm = TRUE)
    )

    return(metrics)
}
