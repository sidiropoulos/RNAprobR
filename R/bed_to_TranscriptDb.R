###Function to transform BED format file to Bioconductor TranscriptDb object
# Description of options:
# input_bed_path - path to BED file (12 column BED needed)
# Value: TranscriptDb object












#' Bedgraph to TranscriptDb object
#' 
#' Function to transform BED format file to Bioconductor TranscriptDb object
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param input_bed_path Path to 12 column BED file
#' @return \item{comp1 }{Description of 'comp1'} \item{comp2 }{Description of
#' 'comp2'}
#' @note %% ~~further notes~~
#' @author Lukasz Jan Kielpinski
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' write("chr1\t134212702\t134229870\tENSMUST00000072177\t0\t+\t134212806\t134228958\t0\t8\t347,121,24,152,66,120,133,1973,\t0,8827,10080,11571,12005,13832,14433,15195,", file="dummy.bed")
#' BED2txDb("dummy.bed")
#' 
#' ## The function is currently defined as
#' function (input_bed_path)
#' {
#'     rank_fun <- function(exon_count, strandness) {
#'         ranks <- 1:exon_count
#'         if (strandness == "-") {
#'             ranks <- rev(ranks)
#'         }
#'         ranks
#'     }
#'     input_bed <- read.table(input_bed_path)
#'     if (is.element("tx_id", names(input_bed))) {
#'         print("ERROR: BED file not allowed to contain column named tx_id")
#'         stop()
#'     }
#'     input_bed$tx_id <- (1:nrow(input_bed))
#'     transcripts <- data.frame(tx_id = input_bed$tx_id, tx_name = input_bed[,
#'         4], tx_chrom = input_bed[, 1], tx_strand = input_bed[,
#'         6], tx_start = input_bed[, 2] + 1, tx_end = input_bed[,
#'         3])
#'     input_bed[, 11] <- as.character(input_bed[, 11])
#'     input_bed[, 12] <- as.character(input_bed[, 12])
#'     exon_size_tx_list <- lapply(strsplit(input_bed[, 11], ","),
#'         as.integer)
#'     exon_start_tx_list <- lapply(strsplit(input_bed[, 12], ","),
#'         as.integer)
#'     exon_start_genome_list <- mapply(FUN = function(x, y) {
#'         x + y + 1
#'     }, exon_start_tx_list, input_bed[, 2], SIMPLIFY = FALSE)
#'     exon_end_genome_list <- mapply(FUN = function(x, y) {
#'         x + y - 1
#'     }, exon_start_genome_list, exon_size_tx_list, SIMPLIFY = FALSE)
#'     exon_rank_list <- mapply(FUN = rank_fun, input_bed[, 10],
#'         input_bed[, 6], SIMPLIFY = FALSE)
#'     tx_id_list <- mapply(FUN = rep, input_bed$tx_id, input_bed[,
#'         10], SIMPLIFY = FALSE)
#'     splicings <- data.frame(tx_id = unlist(tx_id_list), exon_rank = unlist(exon_rank_list),
#'         exon_start = unlist(exon_start_genome_list), exon_end = unlist(exon_end_genome_list))
#'     txDb_from_BED <- makeTranscriptDb(transcripts = transcripts,
#'         splicings = splicings)
#'     txDb_from_BED
#'   }
#' 
#' @import GenomicRanges
#' @export BED2txDb
BED2txDb <- function(input_bed_path)
{

    #Define functions:
    rank_fun <- function(exon_count, strandness){
        ranks <- 1:exon_count
        if(strandness=="-"){ranks <- rev(ranks)}
        ranks
    }
    #End of define functions

    ##Read in BED file and add unique IDs:
    input_bed <- read.table(input_bed_path)

    # This condition should never be fulfilled, just in case if someone tampers
    # with BED file.
    if(is.element("tx_id", names(input_bed)))
    {
        print("ERROR: BED file not allowed to contain column named tx_id")
        stop()
    } 
    
    #Add internal transcript ID (required for making TranscriptDb)
    input_bed$tx_id <- (1:nrow(input_bed))

    ##create first data frame required for makeTranscriptDb function 
    #(transcripts):
    transcripts <- data.frame(tx_id=input_bed$tx_id, tx_name=input_bed[,4], 
        tx_chrom=input_bed[,1], tx_strand=input_bed[,6], 
        tx_start=input_bed[,2]+1, tx_end=input_bed[,3])

    ##create second data frame required for makeTranscriptDb function 
    #(splicing):
    #column 11 of BED file - exon sizes;  change factor to character
    input_bed[,11] <- as.character(input_bed[,11])
    #column 12 of BED file - exon starts; change factor to character
    input_bed[,12] <- as.character(input_bed[,12])

    #Split strings into list of integers
    exon_size_tx_list <- lapply(strsplit(input_bed[,11],","), as.integer)
    #Split strings into list of integers  
    exon_start_tx_list <- lapply(strsplit(input_bed[,12],","), as.integer)

    #Calculate exon starts in genome coordinates, and change to 1-counted 
    #notation.
    exon_start_genome_list <- mapply(FUN=function(x,y){x+y+1}, 
        exon_start_tx_list, input_bed[,2], SIMPLIFY = FALSE)
    #Calculate exon ends in genome coordinates, and change to 1-counted 
    #notation.
    exon_end_genome_list <- mapply(FUN=function(x,y){x+y-1}, 
        exon_start_genome_list, exon_size_tx_list, SIMPLIFY = FALSE)

    #Calcualte exon rank: ordering within transcript
    exon_rank_list <- mapply(FUN=rank_fun, input_bed[,10], input_bed[,6], 
        SIMPLIFY = FALSE)

    #transcript ID
    tx_id_list <- mapply(FUN=rep, input_bed$tx_id, input_bed[,10], 
        SIMPLIFY = FALSE)

    splicings <- data.frame(tx_id=unlist(tx_id_list), 
        exon_rank=unlist(exon_rank_list), 
        exon_start=unlist(exon_start_genome_list), 
        exon_end=unlist(exon_end_genome_list))

    ##Make TranscriptDb object:
    txDb_from_BED <- makeTranscriptDb(transcripts=transcripts, 
        splicings=splicings)
    txDb_from_BED
}

