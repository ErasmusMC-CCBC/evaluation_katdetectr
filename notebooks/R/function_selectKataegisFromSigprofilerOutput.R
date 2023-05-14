. <- NULL

selectKataegisFromSigprofilerOutput <- function(pathToData){

    load(file = pathToData)

    samplesWithKataegis <- allResultsSynthetic |>
        # filter out all rows that do not fit our definition of a kataegis foci (totalVariants >= 6 and mean IMD <= 1000)
        # Keep the rows with totalVariants = NA as these represent no detected kataegis in the sample
        filter(.$totalVariants >= 6 & .$meanIMD <= 1000)

    allSamplesCorrect <- allResultsSynthetic |>
        # keep all rows that DO NOT represent a kataegis foci according to our definition
        filter(! .$sampleNames %in% unique(samplesWithKataegis$sampleNames)) |>
        dplyr::group_by(.$sampleNames) |>
        # give each sample in which no kataegis was detected the correct labeling
        dplyr::summarise(seqnames = NA, start = NA, end = NA, totalVariants = NA, meanIMD = NA, sampleNames = sampleNames, runTime = runTime) |>
        # remove all duplications
        dplyr::distinct() |>
        # combine with the samples that do contain kataegis to end up with the original tibble
        dplyr::bind_rows(samplesWithKataegis) |>
        dplyr::ungroup()

    return(allSamplesCorrect)
}
