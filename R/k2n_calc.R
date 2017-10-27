#' Calculate number of Estimated Unique Counts (EUC's)
#' corresponding to given number of observed unique barcodes.
#'
#' Function calculates EUC's for each number of observed barcodes accounting
#' for differential ligation probability of different barcodes. Function
#' k2n_calc() writes file with a vector in which an i-th element is an
#' estimated unique count given observing i unique barcodes.
#'
#' @param merged_file path to merged_temp file containing 4 column: 1) RNAid,
#' 2) Start, 3) End, 4) Barcode sequence (required)
#' @param unique_barcode_file character with path to unique_barcode file
#' (required)
#' @param output_file name of a file to be generated (if specified
#' [recommended] function will write a file, if not - function will return a
#' vector)
#' @return If output_file specified function writes a file, if not - returns a
#' vector.
#' @author Lukasz Jan Kielpinski, Nikos Sidiropoulos
#' @seealso \code{\link{readsamples}}
#' @references Kielpinski, L.J., and Vinther, J. (2014). Massive
#' parallel-sequencing-based hydroxyl radical probing of RNA accessibility.
#' Nucleic Acids Res.
#' @examples
#'
#' write(c("DummyRNA\t1\t1\tA", "DummyRNA\t1\t1\tC", "DummyRNA\t2\t2\tG",
#'         "DummyRNA\t2\t2\tT"),file="dummy_merged_file")
#' write(c("DummyRNA\t1\t1\t2", "DummyRNA\t2\t2\t2"),
#'         file="dummy_unique_barcode")
#' k2n_calc(merged_file = "dummy_merged_file",
#'         unique_barcode_file = "dummy_unique_barcode")
#'
#' @import plyr
#' @importFrom stats quantile
#' @importFrom utils read.table
#' @export k2n_calc
k2n_calc <- function(merged_file, unique_barcode_file, output_file){

    Barcodes <- NULL

    # Read in inputs:
    merged <- read.table( merged_file )
    colnames(merged) <- c("RNAid", "Start", "End", "Barcodes")
    barcode_length <- max(nchar(as.character(merged$Barcodes)))
    read_counts <- count(subset(merged, select = - Barcodes))
    max_observed <- max(read.table(unique_barcode_file)[,4])

    #Check if max_observed equals max possible complexity. If yes - subtract 1
    #(otherwise 'while' loop below never ends)
    if( max_observed == 4**barcode_length ) {
        max_observed <- max_observed - 1
        barcodes_saturated <- TRUE
    } else
        barcodes_saturated <- FALSE

    # Remove top-quartile of reads:
    barcodes_nt <- merged[ do.call(
        paste, as.list(subset(merged, select = - Barcodes))) %in%
            do.call(paste, as.list(read_counts[(
                read_counts$freq ) <= quantile(read_counts$freq, 0.75),
                c("RNAid", "Start", "End") ] ) ) , "Barcodes" ]

    # Convert to array
    barcodes_nt <- matrix(laply(seq(1,barcode_length),
                                function(i) substr(barcodes_nt, i, i)),
                          ncol = barcode_length, byrow = TRUE)

    # make the matrix with the nucleotide freqs per position:
    nt_counts <- apply( barcodes_nt, 2, .count_nt)

    # Calculate frequencies
    nt_freqs <- nt_counts / colSums( nt_counts )

    nt_values <- split(nt_freqs, rep(1:ncol(nt_freqs), each = nrow(nt_freqs)))

    all_posible_comb <- expand.grid( nt_values )

    probs <- apply( all_posible_comb, 1, prod )

    ###Create Mf_to_Sf:
    results <- c()

    i <- 1
    results[ i ] <- sum( 1 - ( 1 - probs )**i )
    j <- 1
    while( results[ i ] <= max_observed ) {
        i <- i + 1
        results[ i ] <- sum( 1 - ( 1 - probs )**i )   #Mf to Sf
        if(results[ i ] == results[ i - 1]){
            if(j > 2)
                stop("Too low complexity to estimate k2n distribution")
            else
                j <- j + 1
        }else
            j <- 1
    }

    #assign molecules number to unique barcode number:
    k2n <- laply( 1:floor(max(results)),
                  function (i) .absolute_difference(results,i))

    #Assign +Inf for fragments which have saturated the barcodes
    #(see correct_oversaturation function):
    if( barcodes_saturated )
        k2n[max_observed + 1] <- Inf

    if(!missing(output_file)) {
        write(k2n, file = output_file )
    } else
        k2n
}

#Count nucleotide occurencies
.count_nt <- function(nt)
{
    nA <- sum(nt=="A")
    nC <- sum(nt=="C")
    nG <- sum(nt=="G")
    nT <- sum(nt=="T")

    c(nA, nC, nG, nT)
}

.absolute_difference <- function( results, i)
{
    difference <- abs( results - i )
    which( difference == min( difference ) )
}