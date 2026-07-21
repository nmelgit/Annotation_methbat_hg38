###############################################################
#
# annotate_methbat_dmrs.R
#
# Annotate MethBat DMR regions (hg38/GRCh38)
#
# Input:
#   MethBat signature output:
#
#   #chrom start end summary_label
#
# Output:
#   Annotated DMR table
#
###############################################################


###############################################################
## USER SETTINGS
###############################################################

## MethBat output file

input_file <- "/scratch/USERS/nmele/nmele/PacBio/MethBat_old_v0.17.0/methbat-v0.17.0-x86_64-unknown-linux-gnu/MS_1_batch_signature_Smoking_no_env.signature_regions.bed"


## GENCODE GTF matching your PacBio alignment

gtf_file <- "/scratch/USERS/nmele/nmele/rres_dmr_analysis/annotation/gencode.v50.annotation.gtf.gz"


## Output

output_file <- "methbat_DMR_annotation.csv"



###############################################################
## Packages
###############################################################

suppressPackageStartupMessages({

    library(GenomicRanges)
    library(txdbmaker)
    library(GenomeInfoDb)
    library(GenomicFeatures)
    library(rtracklayer)
    library(AnnotationHub)
    library(annotatr)
    library(dplyr)
    library(stringr)

})


###############################################################
## Helper function
###############################################################

message("Loading MethBat DMR file...")


read_methbat <- function(file){

    df <- read.delim(
        file,
        comment.char = "#",
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )


    ## Standardize column names

    colnames(df)[1:4] <-
        c(
            "chrom",
            "start",
            "end",
            "summary_label"
        )


    return(df)

}



###############################################################
## Read MethBat output
###############################################################

dmr_df <- read_methbat(input_file)


message(
    "Loaded ",
    nrow(dmr_df),
    " DMRs"
)



###############################################################
## Convert to GRanges
###############################################################

dmr_gr <- GRanges(

    seqnames =
        dmr_df$chrom,

    ranges =
        IRanges(
            start = dmr_df$start,
            end   = dmr_df$end
        ),

    summary_label =
        dmr_df$summary_label

)


## Width

dmr_df$DMR_width <-
    width(dmr_gr)



###############################################################
## Import GENCODE annotation
###############################################################

message("Importing GENCODE annotation...")


gtf <- import(gtf_file)



###############################################################
## Create TxDb object
###############################################################

message("Building transcript database...")


txdb <-
    txdbmaker::makeTxDbFromGRanges(
    gtf,
    drop.stop.codons = FALSE
)



###############################################################
## Extract genomic features
###############################################################

message("Preparing gene models...")


## Genes

genes_gr <-
    gtf[
        gtf$type == "gene"
    ]


## Gene metadata

gene_symbols <-
    mcols(genes_gr)$gene_name


gene_ids <-
    mcols(genes_gr)$gene_id


gene_types <-
    mcols(genes_gr)$gene_type



## Exons

exons_gr <-
    gtf[
        gtf$type == "exon"
    ]



## Promoters
##
## Definition:
##   -2000 bp upstream
##   +500 bp downstream
## relative to TSS

promoters_gr <-
    promoters(
        genes(txdb),
        upstream = 2000,
        downstream = 500
    )




###############################################################
## Introns
###############################################################

message("Building introns...")


introns_gr <-
    unlist(
        intronsByTranscript(txdb)
    )

###############################################################
## Keep annotation objects on standard GRCh38 chromosomes
###############################################################

standard_chromosomes <- paste0(
    "chr",
    c(1:22, "X", "Y", "M")
)


## Remove random/alternative contigs from annotation objects only

genes_gr <- keepSeqlevels(
    genes_gr,
    intersect(seqlevels(genes_gr), standard_chromosomes),
    pruning.mode="coarse"
)


exons_gr <- keepSeqlevels(
    exons_gr,
    intersect(seqlevels(exons_gr), standard_chromosomes),
    pruning.mode="coarse"
)


promoters_gr <- keepSeqlevels(
    promoters_gr,
    intersect(seqlevels(promoters_gr), standard_chromosomes),
    pruning.mode="coarse"
)


introns_gr <- keepSeqlevels(
    introns_gr,
    intersect(seqlevels(introns_gr), standard_chromosomes),
    pruning.mode="coarse"
)

#dmr_gr <- keepSeqlevels(
#    dmr_gr,
#    intersect(seqlevels(dmr_gr), standard_chromosomes),
#    pruning.mode="coarse"
#)

message(
    "Checking DMR count: ",
    length(dmr_gr),
    " regions"
)

###############################################################
## Save objects temporarily
###############################################################

message("Genome annotation objects created")

###############################################################
## Part 2/3
##
## DMR annotation
##
###############################################################


###############################################################
## Helper function:
## collapse multiple hits into one string
###############################################################

collapse_hits <- function(x){

    x <- unique(x)

    if(length(x)==0){
        return(NA)
    }

    paste(x, collapse=";")

}



###############################################################
## Nearest gene annotation (NA-safe)
###############################################################

message("Finding nearest genes...")


## initialize columns

dmr_df$nearest_gene <- NA
dmr_df$nearest_gene_id <- NA
dmr_df$nearest_gene_type <- NA
dmr_df$distance_to_gene <- NA



nearest_hits <-
    distanceToNearest(
        dmr_gr,
        genes_gr
    )


q <- queryHits(nearest_hits)
s <- subjectHits(nearest_hits)


dmr_df$nearest_gene[q] <-
    gene_symbols[s]


dmr_df$nearest_gene_id[q] <-
    gene_ids[s]


dmr_df$nearest_gene_type[q] <-
    gene_types[s]


dmr_df$distance_to_gene[q] <-
    mcols(nearest_hits)$distance



###############################################################
## Distance to nearest TSS
###############################################################


tss_gr <-
    promoters(
        genes(txdb),
        upstream = 0,
        downstream = 1
    )


###############################################################
## Distance to nearest TSS (NA-safe)
###############################################################

message("Calculating distance to TSS...")


dmr_df$distance_to_TSS <- NA


tss_hits <-
    distanceToNearest(
        dmr_gr,
        tss_gr
    )


q <- queryHits(tss_hits)


dmr_df$distance_to_TSS[q] <-
    mcols(tss_hits)$distance



###############################################################
## Overlapping genes
###############################################################

message("Finding overlapping genes...")


gene_overlap <-
    findOverlaps(
        dmr_gr,
        genes_gr
    )


overlap_gene_list <-
    split(
        gene_symbols[
            subjectHits(gene_overlap)
        ],
        queryHits(gene_overlap)
    )


dmr_df$overlapping_genes <- NA


dmr_df$overlapping_genes[
    as.numeric(names(overlap_gene_list))
] <-
    sapply(
        overlap_gene_list,
        collapse_hits
    )



###############################################################
## Genomic feature annotation
###############################################################

message("Annotating genomic regions...")


dmr_df$genomic_feature <-
    "Intergenic"



## Exons

hits <-
    findOverlaps(
        dmr_gr,
        exons_gr
    )


dmr_df$genomic_feature[
    unique(queryHits(hits))
] <- "Exon"



## Introns

hits <-
    findOverlaps(
        dmr_gr,
        introns_gr
    )


idx <-
    unique(queryHits(hits))


dmr_df$genomic_feature[idx] <-
    ifelse(
        dmr_df$genomic_feature[idx]=="Intergenic",
        "Intron",
        dmr_df$genomic_feature[idx]
    )



## Promoters override other categories

hits <-
    findOverlaps(
        dmr_gr,
        promoters_gr
    )


dmr_df$genomic_feature[
    unique(queryHits(hits))
] <- "Promoter"



###############################################################
## CpG context annotation using annotatr
###############################################################

message("Loading CpG context annotation...")


cpg_annotations <-
    build_annotations(
        genome = "hg38",
        annotations = c(
            "hg38_cpg_islands",
            "hg38_cpg_shores",
            "hg38_cpg_shelves"
        )
    )


###############################################################
## Default = OpenSea
###############################################################

dmr_df$CpG_context <- "OpenSea"



###############################################################
## Island
###############################################################

hits <-
    findOverlaps(
        dmr_gr,
        cpg_annotations[
            cpg_annotations$type ==
                "hg38_cpg_islands"
        ]
    )


dmr_df$CpG_context[
    unique(queryHits(hits))
] <- "Island"



###############################################################
## Shores
###############################################################

hits <-
    findOverlaps(
        dmr_gr,
        cpg_annotations[
            cpg_annotations$type ==
                "hg38_cpg_shores"
        ]
    )


idx <- unique(queryHits(hits))


dmr_df$CpG_context[idx] <-
    ifelse(
        dmr_df$CpG_context[idx]=="OpenSea",
        "Shore",
        dmr_df$CpG_context[idx]
    )



###############################################################
## Shelves
###############################################################

hits <-
    findOverlaps(
        dmr_gr,
        cpg_annotations[
            cpg_annotations$type ==
                "hg38_cpg_shelves"
        ]
    )


idx <- unique(queryHits(hits))


dmr_df$CpG_context[idx] <-
    ifelse(
        dmr_df$CpG_context[idx]=="OpenSea",
        "Shelf",
        dmr_df$CpG_context[idx]
    )


message("CpG annotation completed")

###############################################################
## Part 3/3
##
## Finalize annotations
## Export results
##
###############################################################



###############################################################
## Calculate overlap fraction (NA-safe)
###############################################################

message("Calculating overlap fractions...")

calculate_overlap_fraction <- function(query, subject){

    n <- length(query)

    fraction <- rep(0, n)


    hits <- findOverlaps(
        query,
        subject
    )


    if(length(hits) == 0){
        return(fraction)
    }


    overlap_width <- width(
        pintersect(
            query[queryHits(hits)],
            subject[subjectHits(hits)]
        )
    )


    overlap_sum <- tapply(
        overlap_width,
        queryHits(hits),
        sum
    )


    idx <- as.integer(names(overlap_sum))


    fraction[idx] <-
        as.numeric(overlap_sum) /
        width(query)[idx]


    ## safety check
    if(length(fraction) != n){
        stop(
            "Overlap fraction length mismatch: ",
            length(fraction),
            " instead of ",
            n
        )
    }


    return(fraction)

}



###############################################################
## Feature overlap fractions
###############################################################

print(length(dmr_gr))
print(length(calculate_overlap_fraction(dmr_gr, promoters_gr)))

dmr_df$promoter_fraction <-
    calculate_overlap_fraction(
        dmr_gr,
        promoters_gr
    )


dmr_df$exon_fraction <-
    calculate_overlap_fraction(
        dmr_gr,
        exons_gr
    )


dmr_df$intron_fraction <-
    calculate_overlap_fraction(
        dmr_gr,
        introns_gr
    )


###############################################################
## CpG island overlap fraction
###############################################################

cpg_islands_gr <-
    cpg_annotations[
        cpg_annotations$type ==
            "hg38_cpg_islands"
    ]


dmr_df$cpg_island_fraction <-
    calculate_overlap_fraction(
        dmr_gr,
        cpg_islands_gr
    )



## Replace missing fractions

fraction_columns <-
    grep(
        "_fraction$",
        colnames(dmr_df),
        value=TRUE
    )


dmr_df[fraction_columns] <-
    lapply(
        dmr_df[fraction_columns],
        function(x){
            x[is.na(x)] <- 0
            x
        }
    )



###############################################################
## Add useful QC columns
###############################################################

dmr_df$DMR_center <-
    round(
        (dmr_df$start + dmr_df$end) / 2
    )


dmr_df$DMR_width <-
    dmr_df$end -
    dmr_df$start + 1



###############################################################
## Arrange columns
###############################################################

final_columns <- c(

    "chrom",
    "start",
    "end",
    "DMR_center",
    "DMR_width",

    "summary_label",

    "nearest_gene",
    "nearest_gene_id",
    "nearest_gene_type",

    "distance_to_gene",
    "distance_to_TSS",

    "overlapping_genes",

    "genomic_feature",

    "CpG_context",

    "promoter_fraction",
    "exon_fraction",
    "intron_fraction",
    "cpg_island_fraction"

)



## Keep only existing columns

final_columns <-
    final_columns[
        final_columns %in% colnames(dmr_df)
    ]


dmr_annotation <-
    dmr_df[,final_columns]



###############################################################
## Save final annotation
###############################################################

message(
    "Writing annotation table..."
)


write.csv(
    dmr_annotation,
    output_file,
    row.names = FALSE
)



###############################################################
## Summary reports
###############################################################

message("Generating summary tables...")


feature_summary <-
    dmr_annotation %>%
    count(
        genomic_feature
    ) %>%
    arrange(
        desc(n)
    )


cpg_summary <-
    dmr_annotation %>%
    count(
        CpG_context
    ) %>%
    arrange(
        desc(n)
    )


methylation_summary <-
    dmr_annotation %>%
    count(
        summary_label
    )



write.csv(
    feature_summary,
    "DMR_feature_summary.csv",
    row.names=FALSE
)


write.csv(
    cpg_summary,
    "DMR_CpG_context_summary.csv",
    row.names=FALSE
)


write.csv(
    methylation_summary,
    "DMR_direction_summary.csv",
    row.names=FALSE
)



###############################################################
## Console report
###############################################################

cat("\n")
cat("====================================\n")
cat(" MethBat DMR annotation completed\n")
cat("====================================\n\n")


cat(
    "Number of DMRs:",
    nrow(dmr_annotation),
    "\n\n"
)


cat("Genomic feature distribution:\n")
print(feature_summary)


cat("\nCpG context distribution:\n")
print(cpg_summary)


cat("\nMethylation direction:\n")
print(methylation_summary)


cat("\nOutput files:\n")
cat(
    "- ",
    output_file,
    "\n"
)

cat(
    "- DMR_feature_summary.csv\n"
)

cat(
    "- DMR_CpG_context_summary.csv\n"
)

cat(
    "- DMR_direction_summary.csv\n"
)

cat("\nDone!\n")

