#' Bedgraph to TranscriptDb object
#'
#' Function to transform BED format file to Bioconductor TranscriptDb object
#'
#' @param input_bed_path Path to BED file. If 12 column BED provided, function
#' is splice aware. If 6 column BED provided, function assumes no splicing.
#' @return TranscriptDb object
#' @author Lukasz Jan Kielpinski, Nikos Sidiropoulos
#' @examples
#'
#' write(strwrap("chr1\t134212702\t134229870\tENSMUST00000072177\t0\t+
#'              \t134212806\t134228958\t0\t8\t347,121,24,152,66,120,133,1973,
#'              \t0,8827,10080,11571,12005,13832,14433,15195,", width = 300),
#'       file="dummy.bed")
#' BED2txDb("dummy.bed")
#'
#' @import GenomicRanges
#' @export BED2txDb
BED2txDb <- function(input_bed_path)
{

    ##Read in BED file and add unique IDs:
    input_bed <- read.table(input_bed_path)

    #Extend BED to 12 columns if only 6 provided:
    if(ncol(input_bed)==6){

        colnames(input_bed) <- c("chrom", "chromStart", "chromEnd", "name",
                                 "score", "strand")

        input_bed$thickStart <- input_bed$chromStart
        input_bed$thickEnd <- input_bed$chromStart
        input_bed$itemRgb <- 0
        input_bed$blockCount <- 1
        input_bed$blockSizes <- input_bed$chromEnd - input_bed$chromStart
        input_bed$blockStarts <- 0
    } else
        colnames(input_bed) <- c("chrom", "chromStart", "chromEnd", "name",
                                 "score", "strand", "thickStart", "thickEnd",
                                 "itemRgb", "blockCount", "blockSizes",
                                 "blockStarts")

    # This condition should never be fulfilled, just in case if someone tampers
    # with BED file.
    if(is.element("tx_id", names(input_bed)))
        stop("ERROR: BED file not allowed to contain column named tx_id")

    #Add internal transcript ID (required for making TranscriptDb)
    input_bed$tx_id <- (1:nrow(input_bed))

    ##create first data frame required for makeTranscriptDb function
    #(transcripts):
    transcripts <- data.frame(tx_id=input_bed$tx_id, tx_name=input_bed$name,
                              tx_chrom=input_bed$chrom,
                              tx_strand=input_bed$strand,
                              tx_start=input_bed$chromStart + 1,
                              tx_end=input_bed$chromEnd)

    ##create second data frame required for makeTranscriptDb function
    #(splicing):
    #column 11 of BED file - exon sizes;  change factor to character
    input_bed$blockSizes <- as.character(input_bed$blockSizes)
    #column 12 of BED file - exon starts; change factor to character
    input_bed$blockStarts <- as.character(input_bed$blockStarts)

    #Split strings into list of integers
    exon_size_tx_list <- lapply(strsplit(input_bed$blockSizes,","), as.integer)
    #Split strings into list of integers
    exon_start_tx_list <- lapply(strsplit(input_bed$blockStarts,","),
                                 as.integer)

    #Calculate exon starts in genome coordinates, and change to 1-counted
    #notation.
    exon_start_genome_list <- mapply(FUN=function(x,y){x + y + 1},
                                     exon_start_tx_list, input_bed$chromStart,
                                     SIMPLIFY = FALSE)
    #Calculate exon ends in genome coordinates, and change to 1-counted
    #notation.
    exon_end_genome_list <- mapply(FUN=function(x,y){x + y - 1},
                                   exon_start_genome_list, exon_size_tx_list,
                                   SIMPLIFY = FALSE)

    #Calcualte exon rank: ordering within transcript
    exon_rank_list <- mapply(FUN=.rank_fun, input_bed$blockCount,
                             input_bed$strand, SIMPLIFY = FALSE)

    #transcript ID
    tx_id_list <- mapply(FUN=rep, input_bed$tx_id, input_bed$blockCount,
                         SIMPLIFY = FALSE)

    splicings <- data.frame(tx_id=unlist(tx_id_list),
                            exon_rank=unlist(exon_rank_list),
                            exon_start=unlist(exon_start_genome_list),
                            exon_end=unlist(exon_end_genome_list))

    ##Make TranscriptDb object:
    txDb_from_BED <- suppressWarnings(makeTranscriptDb(transcripts=transcripts,
                                      splicings=splicings))
    txDb_from_BED
}

###Auxiliary functions

.rank_fun <- function(exon_count, strandness){
    ranks <- 1:exon_count
    if(strandness=="-")
        ranks <- rev(ranks)

    ranks
}