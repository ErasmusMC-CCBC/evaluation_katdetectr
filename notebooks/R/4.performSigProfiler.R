runSigProfiler <- function(dataAlexandrov = "data/alexandrov_data_processed.RData", dataSynthetic = "data/synthetic_data_processed.RData") {
    dir.create("data/SigProfilerClusters/")
    setwd(dir = "data/SigProfilerClusters/")

    futile.logger::flog.info("Loading Alexandrov et al. (2013) data")
    load(dataAlexandrov)

    futile.logger::flog.info("Running SigProfiler - Alexandrov et al. (2013)")

    allResultsAlexandrov <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariants["LUAD-E01014"], runSigProfilerClustered))
    save(object = allResultsAlexandrov, file = "data/allResultsAlexandrovSigProfiler.Rdata")

    futile.logger::flog.info("Loading synthetic data")
    load(dataSynthetic)

    futile.logger::flog.info("Running SigProfiler - Synthetic dataset")
    allResultsSynthetic <- dplyr::bind_rows(base::lapply(dataSynthetic$genomicVariants, runSigProfilerClustered))
    save(object = allResultsSynthetic, file = "data/allResultsSyntheticSigProfiler.Rdata")
}


runSigProfilerClustered <- function() {
    sampleName <- unique(data@sampleNames)

    # create new dir
    dir.create(path = paste0("", sampleName), showWarnings = FALSE)
    setwd(dir = paste0("", sampleName))

    # convert data to a sigprofiler friendly format
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

    # write data to disk as a .txt file
    write.table(vrFormatted, file = "data_SigProf.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    # add the header to the data.txt
    file.copy("data/header.txt", "data_SigProf_header.txt", overwrite = TRUE)
    file.append("data_SigProf_header.txt", "data_SigProf.txt")
    unlink("data_SigProf.txt")


    startTime <- base::proc.time()
    system("python ../../prunSigProfiler.py")
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