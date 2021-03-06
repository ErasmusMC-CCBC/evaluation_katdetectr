runMaftools <- function(genomicVariants){

    # Write to temp. file.
    tmpMAF <- tempfile(fileext = '.maf')
    utils::write.table(genomicVariants, file = tmpMAF, append = F, quote = F, sep = '\t', row.names = F)

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
        cnLevel = 'all',
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
        ref.build = 'hg19',
        color = NULL,
        savePlot = FALSE,
        width = 6,
        height = 3,
        fontSize = 1.2,
        pointSize = 0.4
    )

    runTime <- base::proc.time() - startTime

    # Check for output (if no kataegis is found, there is no output)
    outFile <- sprintf('%s_Kataegis.tsv', base::unique(genomicVariants$Tumor_Sample_Barcode))

    if(base::file.exists(outFile)){
        results <- readr::read_delim(sprintf('%s_Kataegis.tsv', base::unique(genomicVariants$Tumor_Sample_Barcode)), delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) |>
            dplyr::mutate(
                seqnames = paste0('chr', Chromosome),
                sampleNames = base::unique(genomicVariants$Tumor_Sample_Barcode),
                runTime = runTime[3]
            ) |>
            dplyr::select(
                seqnames, start = Start_Position, end = End_Position, totalVariants = nMuts, meanIMD = Avg_intermutation_dist, sampleNames, runTime
            )
    }else{
        results <- tibble::tibble(seqnames = NA, start = NA, end = NA, totalVariants = NA, meanIMD = NA, sampleNames = base::unique(genomicVariants$Tumor_Sample_Barcode), runTime = runTime[3])
    }

    return(results)
}


runSeqkat <- function(genomicVariants){

    # Convert to SeqKat-friendly formats.
    data <- GenomicRanges::sort(genomicVariants) |>
        tibble::as_tibble(genomicVariants) |>
        dplyr::select(chr = seqnames, location = start, REF = ref, ALT = alt)

    # Write to temp. file.
    tmpFile = file.path(base::tempdir(), paste0(unique(Biobase::sampleNames(genomicVariants)), '.bed'))
    utils::write.table(data, file = tmpFile, append = F, quote = F, sep = '\t', row.names = F, col.names = F)

    if(base::nrow(data) <= 3){
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
        chromosome = 'all',
        chromosome.length.file = NULL,
        trinucleotide.count.file = "./data/trinucleotide_hg19_whole_genome.txt"
    )

    runTime <- base::proc.time() - startTime

    outFile <- list.files(file.path(base::tempdir(), paste0(levels(Biobase::sampleNames(genomicVariants)))), pattern = '_segnum4.txt', all.files = T, recursive = T, full.names = T)

    if(length(outFile) != 0){
        results <- readr::read_delim(outFile, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) |>
            dplyr::mutate(
                seqnames = paste0('chr', .data$chr),
                sampleNames = base::unique(genomicVariants@sampleNames),
                runTime = runTime[3]
            ) |>
            dplyr::mutate(meanIMD = NA) |>
            dplyr::select(
                seqnames, start, end, totalVariants = variants, meanIMD, sampleNames, runTime
            )
    }else{
        results <- tibble::tibble(seqnames = NA, start = NA, end = NA, totalVariants = NA, meanIMD = NA, sampleNames = base::unique(genomicVariants@sampleNames), runTime = runTime[3])
    }
    return(results)
}



runKatdetectr <- function(genomicVariants, minSizeKataegis = 6, maxMeanIMD = 1000, method = "PCF", test.stat = "Exponential", penalty = "BIC", pen.value = 0, minseglen = 2){

    print(unique(genomicVariants@sampleNames))

    startTime <- base::proc.time()

    # Run katdetectr
    kd <- katdetectr::detectKataegis(
        genomicVariants = genomicVariants,
        minSizeKataegis = minSizeKataegis,
        maxMeanIMD = maxMeanIMD,
        test.stat = test.stat,
        penalty = penalty,
        pen.value = pen.value,
        minseglen = minseglen
    )
    # determine the runtime
    runTime <- base::proc.time() - startTime

    # store the results accordingly
    if(base::length(kd@kataegisFoci) != 0){
        results <- kd@kataegisFoci  |>
            tibble::as_tibble() |>
            dplyr::mutate(
                runTime = runTime[3]
            ) |>
            dplyr::select(.data$seqnames, .data$start, .data$end, .data$totalVariants, .data$meanIMD, .data$sampleNames, .data$runTime)
    } else {
        # if no kataegis is detected return a tibble filled with NA and column with parameter settings
        results <- tibble::tibble(seqnames = NA,
                                  start = NA,
                                  end = NA,
                                  totalVariants = NA,
                                  meanIMD = NA,
                                  sampleNames = base::unique(kd@genomicVariants@sampleNames),
                                  runTime = runTime[3])
    }
    return(results)
}


runClusteredMutations <- function(genomicVariants){

    print(i)
    i <<- i+1

    startTime <- base::proc.time()

    print(unique(genomicVariants@sampleNames))

    # convert to format that ClusteredMutations can use
    data <- tibble::as_tibble(genomicVariants) |>
        dplyr::mutate(
            Chr = substr(.data$seqnames, start = 4, stop = 6)
        ) |>
        dplyr::select(Chr, Position = start, Ref_base = ref, Mutant_base = alt, Sample_id = sampleNames)

    results <- ClusteredMutations::showers(data = data, chr = Chr, position = Position)

    runTime <- base::proc.time() - startTime

    if(base::nrow(results) != 0){
        results <- results |>
            dplyr::mutate(
                seqnames = base::ifelse(!base::is.na(chr), base::paste0('chr', chr), chr),
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


runKataegis <- function(genomicVariants){

    print(unique(genomicVariants@sampleNames))

    # Convert data to correct format
    data <- tibble::as_tibble(genomicVariants) |>
        dplyr::mutate(
            '#CHROM' = base::as.character(.data$seqnames),
            POS = base::trimws(base::as.character(.data$start)),
            FILTER = 'PASS'
        ) |>
        dplyr::select(
            '#CHROM',
            POS,
            REF = .data$ref,
            ALT = .data$alt,
            FILTER
        ) |>
        base::as.matrix()

    startTime <- base::proc.time()

    # kataegis needs more than 1 mutation or throws an error
    if(base::nrow(data) < 2){
        seg <- tibble::tibble(chr = NA, arm = NA, start.pos = NA, end.pos = NA, n.mut = NA, dist.mean = NA)
    } else{
        # Calculate IMD.
        IMD <- kataegis::interMutDist(data)

        # Detect kataegis; try-catch is required as it gives an error if no kataegis is detected!?
        seg <- tryCatch({
            kataegis::kata(IMD, kmin = 2, gamma = 25, assembly = 'hg19', len = 1000, nmut = 6, verbose = TRUE)
        }, error = function(e) {
            tibble::tibble(chr = NA, arm = NA, start.pos = NA, end.pos = NA, n.mut = NA, dist.mean = NA)
        })
    }

    runTime <- base::proc.time() - startTime

    # Clean-up results.
    results <- tibble::as_tibble(seg) |>
        dplyr::mutate(
            # Add chr-prefix.
            chr = base::ifelse(!base::is.na(chr), base::paste0('chr', chr), chr),
            sampleNames = base::unique(tibble::as_tibble(genomicVariants)$sampleNames),
            runTime = runTime[3]
        ) |>
        # remove arm column
        dplyr::select(seqnames = chr, start = start.pos, end = end.pos, totalVariants = n.mut, meanIMD = dist.mean, sampleNames, runTime)

    return(results)
}
