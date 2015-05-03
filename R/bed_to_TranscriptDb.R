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
    input_bed <- rtracklayer::import(input_bed_path, format="bed")

    ##Check if blocks are defined. If not - assume no splicing.
    if(!is.element("blocks", colnames(mcols(input_bed)))){
        input_bed$blocks <- split(IRanges(start=1, end=width(input_bed)), 1:length(input_bed))
        warning("BED file do not contain splicing information. No splicing assumed.")
    }

    # This condition should never be fulfilled, just in case if someone tampers
    # with BED file.
    if(is.element("tx_id", colnames(mcols(input_bed))))
        stop("ERROR: BED file not allowed to contain column named tx_id")

    #Add internal transcript ID (required for making TranscriptDb)
    input_bed$tx_id <- 1:length(input_bed)

    ##create first data frame required for makeTxDb function
    #(transcripts):
    transcripts <- data.frame(tx_id=input_bed$tx_id, tx_name=input_bed$name,
                              tx_chrom=seqnames(input_bed),
                              tx_strand=strand(input_bed),
                              tx_start=start(input_bed),
                              tx_end=end(input_bed))

    ##create second data frame required for makeTxDb function
    #(splicing):

    #Split strings into list of integers
    exon_size_tx_list <- width(input_bed$blocks)
    #Split strings into list of integers
    exon_start_tx_list <- start(input_bed$blocks)
    exon_end_tx_list <- end(input_bed$blocks)

    #Calculate exon starts in genome coordinates, and change to 1-counted
    #notation.
    exon_start_genome_list <- mapply(FUN=function(x,y){x + y - 1},
                                     exon_start_tx_list, start(input_bed),
                                     SIMPLIFY = FALSE)
    #Calculate exon ends in genome coordinates, and change to 1-counted
    #notation.
    exon_end_genome_list <- mapply(FUN=function(x,y){x + y - 1},
                                   exon_end_tx_list, start(input_bed),
                                   SIMPLIFY = FALSE)

    #Calculate exon rank: ordering within transcript
    blockCount <- unlist(lapply(blocks(input_bed), length))
    exon_rank_list <- mapply(FUN=.rank_fun, exon_count = blockCount,
                             strandness = as.character(strand(input_bed)), SIMPLIFY = FALSE)

    #transcript ID
    tx_id_list <- mapply(FUN=rep, input_bed$tx_id, blockCount,
                         SIMPLIFY = FALSE)

    splicings <- data.frame(tx_id=unlist(tx_id_list),
                            exon_rank=unlist(exon_rank_list),
                            exon_start=unlist(exon_start_genome_list),
                            exon_end=unlist(exon_end_genome_list))

    ##Make TranscriptDb object:
    suppressWarnings(makeTxDb(transcripts=transcripts, splicings=splicings))
}

###Auxiliary functions

.rank_fun <- function(exon_count, strandness){

    if(strandness=="-")
        exon_count:1
    else
        1:exon_count
}
