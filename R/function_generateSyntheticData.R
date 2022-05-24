# Author:   Daan Hazelaar / Job van Riet
# Date:     05-05-2022
# Function:

# Import libraries ----

#' @title Generate synthetic dataset
#' @description Generate a synthetic dataset of samples that harbor mutations that resamble kataegis for testing and evaluation purposes.
#'
#' @return list that contains three elements: genomicVariants (all variants in a VRangesList object), genomicVariantsMAF (all variants in a tibble in MAF format) and
#' reportedKataegisFoci (a tibble that contains all kataegis foci)
#'
#' @author Daan Hazelaar

library(katdetectr)
library(dplyr)


generateSyntheticData <- function(){

    syntheticData <- setParameters() |>
        generateData()

    kataegisInformation <- getKataegisInformation(syntheticData)

    syntheticDataMAF <- convertoMafSynth(syntheticData)

    saveData(syntheticData, syntheticDataMAF, kataegisInformation)
}

setParameters <- function(){

    nMutations <- base::round(c(0.1, 0.5, 1, 5, 10, 50, 100, 500) * base::length(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19$chr1) / 10^6)

    # Generate a framework with all possible combinations of parameters.
    parameters <- dplyr::bind_rows(
        tibble::tibble(nBackgroundMuts = base::rep(nMutations, 64),
                       nKataegisFoci = 0,
                       nKataegisMuts = 0,
                       expectedIMD = 0,
                       iteration = base::rep(1:64, each = 8)
        ),
        base::expand.grid(nBackgroundMuts = nMutations,
                          nKataegisFoci = c(1, 2, 3, 5),
                          nKataegisMuts = c(6, 10, 25, 50),
                          expectedIMD = c(100, 250, 500, 750)
        )
    ) %>%
        dplyr::mutate(
            parmNames = sprintf('%s_%s_%s_%s%s', .data$nBackgroundMuts, .data$nKataegisFoci, .data$nKataegisMuts, .data$expectedIMD, base::ifelse(base::is.na(.data$iteration), '', paste0('_(', .data$iteration, ')')))
        )

    return(parameters)
}

generateData <- function(parameters){

    base::set.seed(1)

    allSyntheticDataList <- base::lapply(base::seq_len(base::nrow(parameters)), function(i){

        print(i)

        # Generate synthetic data (SNV only) ----
        syntheticData <- katdetectr::generateSyntheticData(
            genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
            nBackgroundVariants = parameters$nBackgroundMuts[i],
            seqnames = 'chr1',
            probMutationType = c(1,0,0),
            nKataegisFoci = parameters$nKataegisFoci[i],
            nKataegisVariants = parameters$nKataegisMuts[i],
            expectedIMD = parameters$expectedIMD[i],
            sampleName = parameters$parmNames[i],
            removeValidationColumns = FALSE
        )

        return(syntheticData)
    })

    base::set.seed(NULL)

    # Convert to VRanges.
    allSyntheticDataList <- VariantAnnotation::VRangesList(allSyntheticDataList)

    # Add sample name.
    sampleNames <- BiocGenerics::sapply(allSyntheticDataList, function(x){unique(x@sampleNames)})
    base::names(allSyntheticDataList) <- sampleNames

    return(allSyntheticDataList)
}


getKataegisInformation <- function(syntheticData){

    allSyntheticKataegisFoci <- dplyr::bind_rows(BiocGenerics::lapply(syntheticData, function(x){

        # Add info into a GRanges.
        if(any(x$kataegis)){

            syntheticKataegisFoci <- GenomicRanges::GRanges(
                seqnames = base::unique(x@seqnames),
                ranges = IRanges::IRanges(start = unique(na.omit(x$startKataegisFoci)), end = unique(na.omit(x$endKataegisFoci))),
                sampleNames = unique(x@sampleNames)
            )

            # Add column specifying the nr. of variants in kataegis foci.
            syntheticKataegisFoci$totalVariants <- IRanges::countOverlaps(syntheticKataegisFoci, x)

            # Add column specifying the presence of kataegis.
            syntheticKataegisFoci$kataegis <- TRUE

            syntheticKataegisFoci <- tibble::as_tibble(syntheticKataegisFoci)

        } else {
            # If no kataegis was added, return NAs.
            syntheticKataegisFoci <- tibble::tibble(
                seqnames = NA, start = NA, end = NA, width = NA, strand = NA, totalVariants = NA,
                sampleNames = base::unique(x@sampleNames), kataegis = FALSE
            )
        }

        return(syntheticKataegisFoci)

    }))

    return(allSyntheticKataegisFoci)
}


convertoMafSynth <- function(syntheticData){

    # Convert to MAF format (for maftools) ----
    syntheticDataMAF <- tibble::as_tibble(base::unlist(syntheticData)) %>%
        dplyr::mutate(
            Hugo_Symbol = 'Unknown',
            Variant_Classification = 'Missense_Mutation',
            Variant_Type = 'SNP',
        ) %>%
        dplyr::select(
            Hugo_Symbol,
            Chromosome = .data$seqnames,
            Start_Position = .data$start,
            End_Position = .data$end,
            Reference_Allele = .data$ref,
            Tumor_Seq_Allele2 = .data$alt,
            Variant_Classification,
            Variant_Type,
            Tumor_Sample_Barcode = .data$sampleNames
        ) |>
        dplyr::group_by(Tumor_Sample_Barcode) |>
        dplyr::group_split()


    return(syntheticDataMAF)
}


saveData <- function(syntheticData, syntheticDataMAF, kataegisInformation){

    dataSynthetic <- list(
        genomicVariants = syntheticData,
        genomicVariantsMAF = syntheticDataMAF,
        reportedKataegisFoci = kataegisInformation
    )

    save(dataSynthetic, file = "~/local_data/katdetectr/data/synthetic_data.RData")
}
