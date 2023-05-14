# counter to keep track of progress
i <- 1

# Run all tools on Alexandrov et al. dataset
runTools_Alexandrov <- function(data = "data/alexandrov_data_processed.RData") {
    futile.logger::flog.info("Running all R-based tools on Alexandrov et al. (2013) dataset")

    futile.logger::flog.info("Loading RData object: Alexandrov et al. (2013) dataset")
    alexandrovData <- NULL
    load(file = data, verbose = FALSE)

    futile.logger::flog.info("Running katdetectr")
    resultsAlexKatdetectr <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariants, runKatdetectr, minSizeKataegis = 6, test.stat = "Exponential", penalty = "BIC", pen.value = 0, minseglen = 2))
    save(resultsAlexKatdetectr, file = "./data/resultsAlexKatdetectr.RData")

    futile.logger::flog.info("Running SeqKat")
    resultsAlexSeqkat <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariants, runSeqkat))
    save(resultsAlexSeqkat, file = "data/resultsAlexSeqkat.RData")

    futile.logger::flog.info("Running MafTools")
    resultsAlexMaftools <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariantsMAF, runMaftools))
    save(resultsAlexMaftools, file = "data/resultsAlexMaftools.RData")

    futile.logger::flog.info("Running kataegis")
    resultsAlexKataegis <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariants, runKataegis))
    save(resultsAlexKataegis, file = "data/resultsAlexKataegis.RData")

    futile.logger::flog.info("Running ClusteredMutations")
    resultsAlexClusteredMutations <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariants, runClusteredMutations))
    save(resultsAlexClusteredMutations, file = "data/resultsAlexClusteredMutations.RData")

    futile.logger::flog.info("Combining results and saving to data/resultsAlexandrovAllTools.RData")

    resultsAlexandrovAllTools <- base::list(
        katdetectr = resultsAlexKatdetectr,
        seqkat = resultsAlexSeqkat,
        maftools = resultsAlexMaftools,
        clusteredMutations = resultsAlexClusteredMutations,
        kataegis = resultsAlexKataegis
    )

    save(resultsAlexandrovAllTools, file = "data/resultsAlexandrovAllTools.RData")
}


runTools_Synthetic <- function(data = "data/synthetic_data.RData") {
    futile.logger::flog.info("Running all R-based tools on Alexandrov et al. (2013) dataset")

    futile.logger::flog.info("Loading RData object: synthetic dataset")
    dataSynthetic <- NULL
    load(data)

    futile.logger::flog.info("Running katdetectr")
    # resultsSyntheticKatdetectrAMOC <- dplyr::bind_rows(lapply(dataSynthetic$genomicVariants, runKatdetectr, method = "AMOC"))
    # save(resultsSyntheticKatdetectrAMOC, file = "./data/resultsSyntheticKatdetectrAMOC.RData")

    # resultsSyntheticKatdetectrSegNeigh <- dplyr::bind_rows(lapply(dataSynthetic$genomicVariants, runKatdetectr, method = "SegNeigh"))
    # save(resultsSyntheticKatdetectrSegNeigh, file = "./data/resultsSyntheticKatdetectrSegNeigh.RData")

    # resultsSyntheticKatdetectrBinSeg <- dplyr::bind_rows(lapply(dataSynthetic$genomicVariants, runKatdetectr, method = "BinSeg"))
    # save(resultsSyntheticKatdetectrBinSeg, file = "./data/resultsSyntheticKatdetectrBinSeg.RData")

    resultsSyntheticKatdetectr <- dplyr::bind_rows(lapply(dataSynthetic$genomicVariants, runKatdetectr))
    save(resultsSyntheticKatdetectr, file = "data/resultsSyntheticKatdetectr.RData")

    futile.logger::flog.info("Running SeqKat")
    resultsSyntheticSeqkat <- dplyr::bind_rows(lapply(dataSynthetic$genomicVariants, runSeqkat))
    save(resultsSyntheticSeqkat, file = "data/resultsSyntheticSeqkat.RData")

    futile.logger::flog.info("Running MafTools")
    resultsSyntheticMaftools <- dplyr::bind_rows(lapply(dataSynthetic$genomicVariantsMAF, runMaftools))
    save(resultsSyntheticMaftools, file = "data/resultsSyntheticMaftools.RData")

    futile.logger::flog.info("Running kataegis")
    resultsSyntheticKataegis <- dplyr::bind_rows(lapply(dataSynthetic$genomicVariants, runKataegis))
    save(resultsSyntheticKataegis, file = "data/resultsSyntheticKataegis.RData")

    futile.logger::flog.info("Running ClsuteredMutations")
    resultsSyntheticClusteredMutations <- dplyr::bind_rows(lapply(dataSynthetic$genomicVariants, runClusteredMutations))
    save(resultsSyntheticClusteredMutations, file = "data/resultsSyntheticClusteredMutations.RData")

    futile.logger::flog.info("Combining results and saving to data/resultsSyntheticAllTools.RData")

    resultsSyntheticAllTools <- base::list(
        katdetectr = resultsSyntheticKatdetectr,
        seqkat = resultsSyntheticSeqkat,
        maftools = resultsSyntheticMaftools,
        clusteredMutations = resultsSyntheticClusteredMutations,
        kataegis = resultsSyntheticKataegis
    )

    save(resultsSyntheticAllTools, file = "data/resultsSyntheticAllTools.RData")
}

# Functions to run the various packages. ----

runMaftools <- function(genomicVariants) {
    # Write to temp. file.
    tmpMAF <- tempfile(fileext = ".maf")
    utils::write.table(genomicVariants, file = tmpMAF, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

    # Import MAF to maftools object.
    MAF <- maftools::read.maf(
        maf = tmpMAF,
        clinicalData = NULL,
        rmFlags = FALSE,
        removeDuplicatedVariants = TRUE,
        useAll = TRUE,
        gisticAllLesionsFile = NULL,
        gisticAmpGenesFile = NULL,
        gisticDelGenesFile = NULL,
        gisticScoresFile = NULL,
        cnLevel = "all",
        cnTable = NULL,
        isTCGA = FALSE,
        vc_nonSyn = NULL,
        verbose = FALSE
    )

    startTime <- base::proc.time()

    # maftools saves the result/output to the working directory.
    maftools::rainfallPlot(
        maf = MAF,
        detectChangePoints = TRUE,
        ref.build = "hg19",
        color = NULL,
        savePlot = FALSE,
        width = 6,
        height = 3,
        fontSize = 1.2,
        pointSize = 0.4
    )

    runTime <- base::proc.time() - startTime

    # Check for output (if no kataegis is found, there is no output)
    outFile <- sprintf("%s_Kataegis.tsv", base::unique(genomicVariants$Tumor_Sample_Barcode))

    if (base::file.exists(outFile)) {
        results <- readr::read_delim(sprintf("%s_Kataegis.tsv", base::unique(genomicVariants$Tumor_Sample_Barcode)), delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) |>
            dplyr::mutate(
                seqnames = paste0("chr", Chromosome),
                sampleNames = base::unique(genomicVariants$Tumor_Sample_Barcode),
                runTime = runTime[3]
            ) |>
            dplyr::select(
                seqnames,
                start = Start_Position, end = End_Position, totalVariants = nMuts, meanIMD = Avg_intermutation_dist, sampleNames, runTime
            )
    } else {
        results <- tibble::tibble(seqnames = NA, start = NA, end = NA, totalVariants = NA, meanIMD = NA, sampleNames = base::unique(genomicVariants$Tumor_Sample_Barcode), runTime = runTime[3])
    }

    return(results)
}


runSeqkat <- function(genomicVariants) {

    # Convert to SeqKat-friendly formats.
    data <- GenomicRanges::sort(genomicVariants) |>
        tibble::as_tibble() |>
        dplyr::select(chr = seqnames, location = start, REF = ref, ALT = alt)

    # Write to temp. file.
    tmpFile <- file.path(base::tempdir(), paste0(unique(Biobase::sampleNames(genomicVariants)), ".bed"))
    utils::write.table(data, file = tmpFile, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    if (base::nrow(data) <= 3) {
        return(tibble::tibble(seqnames = NA, start = NA, end = NA, totalVariants = NA, meanIMD = NA, sampleNames = base::unique(genomicVariants@sampleNames), runTime = NA))
    }

    startTime <- base::proc.time()

    SeqKat::seqkat(
        sigcutoff = 5,
        mutdistance = 3.2,
        segnum = 4,
        ref.dir = "./data/perChromosome/",
        bed.file = tmpFile,
        output.dir = file.path(base::tempdir(), paste0(unique(Biobase::sampleNames(genomicVariants)))),
        chromosome = "all",
        chromosome.length.file = NULL,
        trinucleotide.count.file = "misc/trinucleotide_hg19_whole_genome.txt"
    )

    runTime <- base::proc.time() - startTime

    outFile <- list.files(file.path(base::tempdir(), paste0(levels(Biobase::sampleNames(genomicVariants)))), pattern = "_segnum4.txt", all.files = TRUE, recursive = TRUE, full.names = TRUE)

    if (length(outFile) != 0) {
        results <- readr::read_delim(outFile, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) |>
            dplyr::mutate(
                seqnames = paste0("chr", chr),
                sampleNames = base::unique(genomicVariants@sampleNames),
                runTime = runTime[3]
            ) |>
            dplyr::mutate(meanIMD = NA) |>
            dplyr::select(
                seqnames, start, end,
                totalVariants = variants, meanIMD, sampleNames, runTime
            )
    } else {
        results <- tibble::tibble(seqnames = NA, start = NA, end = NA, totalVariants = NA, meanIMD = NA, sampleNames = base::unique(genomicVariants@sampleNames), runTime = runTime[3])
    }
    return(results)
}



runKatdetectr <- function(genomicVariants, minSizeKataegis = 6, test.stat = "Exponential", penalty = "BIC", pen.value = 0, minseglen = 2) {
    futile.logger::flog.trace(unique(genomicVariants@sampleNames))

    startTime <- base::proc.time()

    # Run katdetectr
    kd <- katdetectr::detectKataegis(
        genomicVariants = genomicVariants,
        minSizeKataegis = minSizeKataegis,
        test.stat = test.stat,
        penalty = penalty,
        pen.value = pen.value,
        minseglen = minseglen,
        BPPARAM = BiocParallel::MulticoreParam(workers = 4, progressbar = TRUE)
    )
    # determine the runtime
    runTime <- base::proc.time() - startTime

    # store the results accordingly
    if (base::length(kd@kataegisFoci) != 0) {
        results <- kd@kataegisFoci |>
            tibble::as_tibble() |>
            dplyr::mutate(
                runTime = runTime[3]
            ) |>
            dplyr::select(seqnames, start, end, totalVariants, meanIMD, sampleNames, runTime)
    } else {
        # if no kataegis is detected return a tibble filled with NA and column with parameter settings
        results <- tibble::tibble(
            seqnames = NA,
            start = NA,
            end = NA,
            totalVariants = NA,
            meanIMD = NA,
            sampleNames = base::unique(kd@genomicVariants@sampleNames),
            runTime = runTime[3]
        )
    }
    return(results)
}


runClusteredMutations <- function(genomicVariants) {
    print(i)
    i <<- i + 1

    startTime <- base::proc.time()

    print(unique(genomicVariants@sampleNames))

    # convert to format that ClusteredMutations can use
    data <- tibble::as_tibble(genomicVariants) |>
        dplyr::mutate(
            Chr = substr(seqnames, start = 4, stop = 6)
        ) |>
        dplyr::select(Chr, Position = start, Ref_base = ref, Mutant_base = alt, Sample_id = sampleNames)

    results <- ClusteredMutations::showers(data = data, chr = Chr, position = Position)

    runTime <- base::proc.time() - startTime

    if (base::nrow(results) != 0) {
        results <- results |>
            dplyr::mutate(
                seqnames = base::ifelse(!base::is.na(chr), base::paste0("chr", chr), chr),
                runTime = runTime[3],
                sampleNames = unique(genomicVariants@sampleNames),
                meanIMD = NA
            ) |>
            dplyr::select(seqnames, start = pstart, end = pend, totalVariants = number, sampleNames, runTime)
    } else {
        # if no kataegis is detected return a tibble filled with NA and column with parameter settings
        results <- tibble::tibble(seqnames = NA, start = NA, end = NA, totalVariants = NA, meanIMD = NA, sampleNames = unique(genomicVariants@sampleNames), runTime = runTime[3])
    }
    return(results)
}


runKataegis <- function(genomicVariants) {
    print(unique(genomicVariants@sampleNames))

    # Convert data to correct format
    data <- tibble::as_tibble(genomicVariants) |>
        dplyr::mutate(
            "#CHROM" = base::as.character(seqnames),
            POS = base::trimws(base::as.character(start)),
            FILTER = "PASS"
        ) |>
        dplyr::select(
            "#CHROM",
            POS,
            REF = ref,
            ALT = alt,
            FILTER
        ) |>
        base::as.matrix()

    startTime <- base::proc.time()

    # kataegis needs more than 1 mutation or throws an error
    if (base::nrow(data) < 2) {
        seg <- tibble::tibble(chr = NA, arm = NA, start.pos = NA, end.pos = NA, n.mut = NA, dist.mean = NA)
    } else {
        # Calculate IMD.
        IMD <- kataegis::interMutDist(data)

        # Detect kataegis; try-catch is required as it gives an error if no kataegis is detected!?
        seg <- tryCatch(
            {
                kataegis::kata(IMD, kmin = 2, gamma = 25, assembly = "hg19", len = 1000, nmut = 6, verbose = TRUE)
            },
            error = function(e) {
                tibble::tibble(chr = NA, arm = NA, start.pos = NA, end.pos = NA, n.mut = NA, dist.mean = NA)
            }
        )
    }

    runTime <- base::proc.time() - startTime

    # Clean-up results.
    results <- tibble::as_tibble(seg) |>
        dplyr::mutate(
            # Add chr-prefix.
            chr = base::ifelse(!base::is.na(chr), base::paste0("chr", chr), chr),
            sampleNames = base::unique(tibble::as_tibble(genomicVariants)$sampleNames),
            runTime = runTime[3]
        ) |>
        # remove arm column
        dplyr::select(seqnames = chr, start = start.pos, end = end.pos, totalVariants = n.mut, meanIMD = dist.mean, sampleNames, runTime)

    return(results)
}