library("VariantAnnotation")
library("dplyr")

setwd(dir = "~/sig/")
load(file = "data/synthetic_data.RData")
load(file = "data/alexandrov_data_processed.RData")

setwd(dir = "sigProfilerResults")

runSigProfilerClustered <- function(data){

    sampleName <- unique(data@sampleNames)

    # create new dir
    dir.create(path = paste0("", sampleName), showWarnings = FALSE)
    setwd(dir = paste0("", sampleName))

    # convert data to a sigprofiler friendly format
    GenomeInfoDb::seqlevelsStyle(data) <- 'NCBI'
    vrFormatted <- data |> tibble::as_tibble() |>
        mutate(
            CancerType = ".",
            referenceGenome = "GRCh37",
            NGS_method = ".",
            mutss = "SOMATIC",
            mutType = "SNV"
        ) |>
        select(CancerType, sampleNames, NGS_method, referenceGenome, mutType, seqnames, start, end, ref, alt, mutss)

    # write data to disk as a .txt file
    write.table(vrFormatted, file = "data_SigProf.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    # add the header to the data.txt
    file.copy("data/header.txt", "data_SigProf_header.txt", overwrite = TRUE)
    file.append("data_SigProf_header.txt", "data_SigProf.txt")
    unlink("data_SigProf.txt")


    startTime <- base::proc.time()
    system("python ../../prunSigProfiler.py")
    runTime <- base::proc.time() - startTime

    if(file.exists("./output/vcf_files_corrected/results_SigProfiler_clustered/subclasses/class2/results_SigProfiler_clustered_class2.txt")){
        if(nrow(readr::read_table("./output/vcf_files_corrected/results_SigProfiler_clustered/subclasses/class2/results_SigProfiler_clustered_class2.txt")) != 0){
            results <- readr::read_table("./output/vcf_files_corrected/results_SigProfiler_clustered/subclasses/class2/results_SigProfiler_clustered_class2.txt") |>
                select(seqnames = chr, sampleNames = samples, start, end, ref, alt, fociID = IMDplot, IMD) |>
                group_by(fociID) |>
                summarise(
                    sampleNames = as.character(unique(sampleNames)),
                    start = min(start),
                    end = max(end),
                    totalVariants = sum(n()),
                    meanIMD = mean(.data$IMD, na.rm = TRUE),
                    seqnames = paste0("chr", unique(seqnames)),
                    runTime = runTime[3]
                ) |>
                select(seqnames, start, end, totalVariants, meanIMD, sampleNames, runTime)
        }else {
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

allResultsAlexandrov <- bind_rows(lapply(alexandrovData$genomicVariants["LUAD-E01014"], runSigProfilerClustered))
save(object = allResultsAlexandrov, file = "data/allResultsAlexandrovSigProfiler.Rdata")

allResultsSynthetic <- bind_rows(lapply(dataSynthetic$genomicVariants, runSigProfilerClustered))
save(object = allResultsSynthetic, file = "data/allResultsSyntheticSigProfiler.Rdata")
