runSigProfiler <- function(pathAlexandrov = "data/alexandrov_data_processed.RData", pathSynthetic = "data/synthetic_data.RData", subset = FALSE) {

    futile.logger::flog.info("Loading Alexandrov et al. (2013) data")
    load(pathAlexandrov)

    futile.logger::flog.info("Loading synthetic data")
    load(pathSynthetic)

    if (subset) {
        dataSynthetic$genomicVariants <- dataSynthetic$genomicVariants[1:5]
    }

    # Subset the samples to first 5 samples.
    if (subset) {
        alexandrovData$genomicVariants <- alexandrovData$genomicVariants[1:5]
    }

    # Create the output and tmp folder to run in.
    dir.create("data/SigProfilerClusters/")
    setwd(dir = "data/SigProfilerClusters/")

    futile.logger::flog.info("Running SigProfiler - Alexandrov et al. (2013)")

    allResultsAlexandrov <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariants, runSigProfilerClustered))
    save(object = allResultsAlexandrov, file = "data/allResultsAlexandrovSigProfiler.Rdata")

    futile.logger::flog.info("Running SigProfiler - Synthetic dataset")
    allResultsSynthetic <- dplyr::bind_rows(base::lapply(dataSynthetic$genomicVariants, runSigProfilerClustered))
    save(object = allResultsSynthetic, file = "data/allResultsSyntheticSigProfiler.Rdata")
}


runSigProfilerClustered <- function(data) {

    sampleName <- unique(data@sampleNames)

    # Create a directory for the sample and move into it.
    dir.create(path = paste0("", sampleName), showWarnings = FALSE)
    setwd(dir = paste0("", sampleName))

    # Convert data to a sigprofiler friendly format
    GenomeInfoDb::seqlevelsStyle(data) <- "NCBI"
    vrFormatted <- data |>
        tibble::as_tibble() |>
        dplyr::mutate(
            CancerType = ".",
            referenceGenome = "GRCh37",
            NGS_method = ".",
            mutss = "SOMATIC",
            mutType = "SNV"
        ) |>
        dplyr::select(CancerType, sampleNames, NGS_method, referenceGenome, mutType, seqnames, start, end, ref, alt, mutss)

    # Write to a .txt file.
    write.table(vrFormatted, file = "data_SigProf.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    startTime <- base::proc.time()
    system("python3 ../../../../python/prunSigProfiler.py")
    runTime <- base::proc.time() - startTime

    if (file.exists("./output/vcf_files_corrected/results_SigProfiler_clustered/subclasses/class2/results_SigProfiler_clustered_class2.txt")) {
        if (nrow(readr::read_table("./output/vcf_files_corrected/results_SigProfiler_clustered/subclasses/class2/results_SigProfiler_clustered_class2.txt")) != 0) {
            results <- readr::read_table("./output/vcf_files_corrected/results_SigProfiler_clustered/subclasses/class2/results_SigProfiler_clustered_class2.txt") |>
                dplyr::select(seqnames = chr, sampleNames = samples, start, end, ref, alt, fociID = IMDplot, IMD) |>
                dplyr::group_by(fociID) |>
                dplyr::summarise(
                    sampleNames = as.character(unique(sampleNames)),
                    start = min(start),
                    end = max(end),
                    totalVariants = sum(n()),
                    meanIMD = mean(.data$IMD, na.rm = TRUE),
                    seqnames = paste0("chr", unique(seqnames)),
                    runTime = runTime[3]
                ) |>
                dplyr::select(seqnames, start, end, totalVariants, meanIMD, sampleNames, runTime)
        } else {
            results <- tibble::tibble(seqnames = NA, start = NA, end = NA, totalVariants = NA, meanIMD = NA, sampleNames = sampleName, runTime = runTime[3])
        }
    } else {
        results <- tibble::tibble(seqnames = NA, start = NA, end = NA, totalVariants = NA, meanIMD = NA, sampleNames = sampleName, runTime = runTime[3])
    }

    save(object = results, file = paste0("../../results/", sampleName, ".Rdata"))
    setwd(dir = "../")
    system(paste0("rm -rf *"))

    return(results)
}
