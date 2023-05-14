#' @title Import Alexandrov dataset
#' @description Download, process and annotate Alexandrov dataset for further use in performance evaluation of kataegis detection tools.
#' Signatures of mutational processes in human cancer (Nature, 2013). We remove the indels from the data as some tools can only consider SNVs.
#' The data is converted to MAF format for use in maftools. The resulting datasset is also saved to inst/extdata/alexandrovData.RData
#'
#' @return list that contains three elements: genomicVariants (all variants in a VRangesList object), genomicVariantsMAF (all variants in a tibble in MAF format) and
#' reportedKataegisFoci (a tibble that contains all kataegis foci reported by Alexandrov et al.)
#'
#' @author Daan Hazelaar
#' @author Job van Riet

. <- NULL
.data <- NULL

importAlexandrovData <- function(path = "notebooks/data/") {
    futile.logger::flog.info("Downloading somatic variants - Alexandrov et al. (2013)")

    AlexandrovDataProcessed <- getAlexandrovLinks() |>
        downloadAlexandrovData() |>
        processAlexandrovData()

    futile.logger::flog.info("Downloading and importing kataegis calls - Alexandrov et al. (2013)")

    reportedKataegisFoci <- importReportedKataegisFoci()

    futile.logger::flog.info("Annotating samples with kataegis calls - Alexandrov et al. (2013)")

    AlexandrovDataAnnotated <- addSequencingTechnique(AlexandrovDataProcessed) |>
        dplyr::filter(.data$sequencingType == "Whole_genome") |>
        annotateAlexandrovData(reportedKataegisFoci = reportedKataegisFoci)

    futile.logger::flog.info("Converting to MAF files - Alexandrov et al. (2013)")

    AlexandrovDataMaf <- convertToMaf(AlexandrovDataAnnotated)

    alexandrovData <- base::list(
        genomicVariants = AlexandrovDataAnnotated,
        genomicVariantsTibble = tibble::as_tibble(AlexandrovDataAnnotated),
        genomicVariantsMAF = AlexandrovDataMaf,
        reportedKataegisFoci = reportedKataegisFoci
    )

    futile.logger::flog.info("Saving processed Alexandrov et al. (2013) data under: %s/alexandrov_data_processed.RData", path)

    base::save(alexandrovData, file = file.path(path, "alexandrov_data_processed.RData"))

    futile.logger::flog.info("Done")
}


getAlexandrovLinks <- function() {
    linksAlexandrov <- c(
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/ALL/ALL_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/AML/AML_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Bladder/Bladder_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Breast/Breast_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Cervix/Cervix_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/CLL/CLL_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Colorectum/Colorectum_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Esophageal/Esophageal_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Glioblastoma/Glioblastoma_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Glioma%20Low%20Grade/Glioma%20Low%20Grade_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Head%20and%20Neck/Head%20and%20Neck_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Kidney%20Chromophobe/Kidney%20Chromophobe_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Kidney%20Clear%20Cell/Kidney%20Clear%20Cell_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Kidney%20Papillary/Kidney%20Papillary_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Liver/Liver_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lung%20Adeno/Lung%20Adeno_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lung%20Small%20Cell/Lung%20Small%20Cell_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lung%20Squamous/Lung%20Squamous_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma%20B-cell/Lymphoma%20B-cell_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Medulloblastoma/Medulloblastoma_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Melanoma/Melanoma_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Myeloma/Myeloma_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Neuroblastoma/Neuroblastoma_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Ovary/Ovary_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Pancreas/Pancreas_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Pilocytic%20Astrocytoma/Pilocytic%20Astrocytoma_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Prostate/Prostate_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Stomach/Stomach_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Thyroid/Thyroid_clean_somatic_mutations_for_signature_analysis.txt",
        "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Uterus/Uterus_clean_somatic_mutations_for_signature_analysis.txt"
    )

    return(linksAlexandrov)
}

downloadAlexandrovData <- function(linksAlexandrovData) {

    AlexandrovDataRaw <- base::lapply(linksAlexandrovData, function(link) {
        data <- readr::read_delim(
            file = as.character(link),
            col_names = c("sampleNames", "type", "chr", "start", "end", "REF", "ALT", "origin"),
            trim_ws = TRUE,
            col_types = "cccnnccc",
            delim = "\t"
        ) |>
            # add tissue type (this is specified in the url) the cancer originates from
            dplyr::mutate(
                tissue = base::gsub("%20", " ", base::gsub("_.*", "", base::gsub(".*/", "", link)))
            )

        return(data)
    }) |>
        dplyr::bind_rows()

    return(AlexandrovDataRaw)
}

processAlexandrovData <- function(AlexandrovDataRaw) {
    futile.logger::flog.info("Processing somatic variants - Alexandrov et al. (2013)")

    AlexandrovDataProcessed <- AlexandrovDataRaw |>
        # The following three sample names are duplicated for different samples
        # change the names of these samples
        dplyr::filter(.data$sampleNames == 266 | .data$sampleNames == 325 | .data$sampleNames == 91) |>
        dplyr::group_by(.$tissue) |>
        dplyr::mutate(sampleNames = base::paste0(.data$sampleNames, "_", .data$tissue)) |>
        # add the renamed samples back to the complete dataset
        dplyr::bind_rows(AlexandrovDataRaw[!AlexandrovDataRaw$sampleNames %in% c(266, 325, 91), ]) |>
        # Subset on Single Nucleotide Variants only.
        dplyr::filter(.data$type == "subs")

    return(AlexandrovDataProcessed)
}

importReportedKataegisFoci <- function() {
    # Download the detected kataegis regions from Supplementary Table 3 into a tmp. file.
    link <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature12477/MediaObjects/41586_2013_BFnature12477_MOESM86_ESM.xls"
    tmpFile <- base::tempfile(fileext = ".xls")
    curl::curl_download(link, tmpFile)

    reportedKataegisFoci <- readxl::read_xls(tmpFile, range = "B4:G877") |>
        dplyr::select(sampleNames = .data$`Sample Name`, totalVariants = .data$`No of variants`, dplyr::everything()) |>
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    # Change to UCSC nomenclature.
    GenomeInfoDb::seqlevelsStyle(reportedKataegisFoci) <- "UCSC"

    # add all seqlevels i.e. chromosomes
    GenomeInfoDb::seqlevels(reportedKataegisFoci) <- c(base::paste0(base::rep("chr", 22), base::seq(1:22)), "chrX", "chrY", "chrM")

    return(reportedKataegisFoci)
}

addSequencingTechnique <- function(AlexandrovDataProcessed) {
    sampleOrigin <- readr::read_tsv("misc/sampleOrigin.txt", col_types = 'c')
    AlexandrovDataAnnotated <- dplyr::left_join(AlexandrovDataProcessed, sampleOrigin, by = "sampleNames")

    return(AlexandrovDataAnnotated)
}

annotateAlexandrovData <- function(AlexandrovDataProcessed, reportedKataegisFoci) {

    # convert data to VRanges and add a column which specifies if a variant lies within an kataegis foci
    AlexandrovDataAnnotated <- AlexandrovDataProcessed |>
        # dplyr::filter(sampleNames %in% unique(reportedKataegisFoci$sampleNames)) |>
        dplyr::group_by(.data$sampleNames) |>
        dplyr::group_map(~ annotatePerSample(sampleData = .x, reportedKataegisFoci = reportedKataegisFoci), .keep = TRUE)

    names(AlexandrovDataAnnotated) <- base::unlist(base::lapply(AlexandrovDataAnnotated, function(x) {
        base::unique(x@sampleNames)
    }))
    AlexandrovDataAnnotated <- VariantAnnotation::VRangesList(AlexandrovDataAnnotated)

    return(AlexandrovDataAnnotated)
}

annotatePerSample <- function(sampleData, reportedKataegisFoci) {
    # convert to GRanges object
    sampleDataGr <- GenomicRanges::makeGRangesFromDataFrame(
        df = sampleData,
        seqnames = "chr",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = TRUE
    )
    # Change styling
    GenomeInfoDb::seqlevelsStyle(sampleDataGr) <- "UCSC"
    # convert to VRanges object
    sampleDataVr <- VariantAnnotation::makeVRangesFromGRanges(sampleDataGr)

    # select the kataegis foci that were detected by Alexandrov in this specific sample
    kataegisFociSample <- reportedKataegisFoci[reportedKataegisFoci$sampleNames == base::unique(sampleDataVr@sampleNames)]

    # use overlap functions to find which vairants lie withing a kataegis foci and which variants do not lie withing a kataegis foci
    kataegisVariants <- IRanges::subsetByOverlaps(sampleDataVr, kataegisFociSample)
    NonKataegisVariants <- plyranges::filter_by_non_overlaps(sampleDataVr, kataegisFociSample)

    # add a column which specifies if a variants lies within a kataegis foci (or not)
    if (base::length(kataegisVariants) != 0) {
        kataegisVariants$kataegis <- TRUE
    }
    NonKataegisVariants$kataegis <- FALSE

    # combine the previous two vranges into a single vranges
    sampleDataAnnotated <- c(kataegisVariants, NonKataegisVariants) |>
        GenomicRanges::sort()

    return(sampleDataAnnotated)
}


convertToMaf <- function(AlexandrovDataAnnotated) {
    AlexandrovDataMaf <- AlexandrovDataAnnotated |>
        tibble::as_tibble() |>
        base::suppressWarnings() |>
        dplyr::mutate(
            Hugo_Symbol = "Unknown",
            Variant_Classification = "Missense_Mutation",
            Variant_Type = "SNP"
        ) |>
        dplyr::select(
            .$Hugo_Symbol,
            Chromosome = .$seqnames,
            Start_Position = .$start,
            End_Position = .$end,
            Reference_Allele = .$ref,
            Tumor_Seq_Allele2 = .$alt,
            .$Variant_Classification,
            .$Variant_Type,
            Tumor_Sample_Barcode = .$sampleNames
        ) |>
        dplyr::group_by(.$Tumor_Sample_Barcode) |>
        dplyr::group_split()

    return(AlexandrovDataMaf)
}
