library("VariantAnnotation")
library("dplyr")
library("SeqKat")
source(file = "./R/function_package_workflows.R")

# counter to keep track of progress
i <- 1

# Run all tools on Alexandrov et al. dataset

load(file = "./data/alexandrov_data_processed.RData")

resultsAlexKatdetectr <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariants, runKatdetectr, minSizeKataegis = 6, IMDcutoff = 1000, method = "PCF", test.stat = "Exponential", penalty = "BIC", pen.value = 0, minseglen = 2))
save(resultsAlexKatdetectr, file = "./results/resultsAlexKatdetectr.RData")

resultsAlexSeqkat <- dplyr::bind_rows(pbapply::pblapply(alexandrovData$genomicVariants, runSeqkat, cl = 15))
save(resultsAlexSeqkat, file = "./results/resultsAlexSeqkat.RData")

resultsAlexMaftools <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariantsMAF, runMaftools))
save(resultsAlexMaftools, file = "./results/resultsAlexMaftools.RData")

resultsAlexKataegis <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariants, runKataegis))
save(resultsAlexKataegis, file = "./results/resultsAlexKataegis.RData")

resultsAlexClusteredMutations <- dplyr::bind_rows(base::lapply(alexandrovData$genomicVariants, runClusteredMutations))
save(resultsAlexClusteredMutations, file = "./results/resultsAlexClusteredMutations.RData")

resultsAlexandrovAllTools <- base::list(
    katdetectr = resultsAlexKatdetectr,
    seqkat = resultsAlexSeqkat,
    maftools = resultsAlexMaftools,
    clusteredMutations = resultsAlexClusteredMutations,
    kataegis = resultsAlexKataegis
)

save(resultsAlexandrovAllTools, file = "./results/resultsAlexandrovAllTools.RData")



# Run all tools on the synthetic dataset dataset

load(file = "./data/synthetic_data.RData")

resultsSyntheticKatdetectrAMOC <- dplyr::bind_rows(pbapply::pblapply(dataSynthetic$genomicVariants, runKatdetectr, method = "AMOC", cl = 6))
save(resultsSyntheticKatdetectrAMOC, file = "./results/resultsSyntheticKatdetectrAMOC.RData")

resultsSyntheticKatdetectrSegNeigh <- dplyr::bind_rows(pbapply::pblapply(dataSynthetic$genomicVariants, runKatdetectr, method = "SegNeigh", cl = 6))
save(resultsSyntheticKatdetectrSegNeigh, file = "./results/resultsSyntheticKatdetectrSegNeigh.RData")

resultsSyntheticKatdetectrBinSeg <- dplyr::bind_rows(pbapply::pblapply(dataSynthetic$genomicVariants, runKatdetectr, method = "BinSeg", cl = 6))
save(resultsSyntheticKatdetectrBinSeg, file = "./results/resultsSyntheticKatdetectrBinSeg.RData")

resultsSyntheticKatdetectr <- dplyr::bind_rows(pbapply::pblapply(dataSynthetic$genomicVariants, runKatdetectr, cl = 15))
save(resultsSyntheticKatdetectr, file = "./results/resultsSyntheticKatdetectr.RData")

resultsSyntheticSeqkat <- dplyr::bind_rows(pbapply::pblapply(dataSynthetic$genomicVariants, runSeqkat, cl = 15))
save(resultsSyntheticSeqkat, file = "./results/resultsSyntheticSeqkat.RData")

resultsSyntheticMaftools <- dplyr::bind_rows(pbapply::pblapply(dataSynthetic$genomicVariantsMAF, runMaftools, cl = 15))
save(resultsSyntheticMaftools, file = "./results/resultsSyntheticMaftools.RData")

resultsSyntheticKataegis <- dplyr::bind_rows(pbapply::pblapply(dataSynthetic$genomicVariants, runKataegis, cl = 15))
save(resultsSyntheticKataegis, file = "./results/resultsSyntheticKataegis.RData")

resultsSyntheticClusteredMutations <- dplyr::bind_rows(pbapply::pblapply(dataSynthetic$genomicVariants, runClusteredMutations, cl = 15))
save(resultsSyntheticClusteredMutations, file = "./results/resultsSyntheticClusteredMutations.RData")

resultsSyntheticAllTools <- base::list(
    katdetectr = resultsSyntheticKatdetectr,
    seqkat = resultsSyntheticSeqkat,
    maftools = resultsSyntheticMaftools,
    clusteredMutations = resultsSyntheticClusteredMutations,
    kataegis = resultsSyntheticKataegis
)

save(resultsSyntheticAllTools, file = "./results/resultsSyntheticAllTools.RData")

