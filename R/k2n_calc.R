#' Function to calculate number of Estimated Unique Counts (EUC's)
#' corresponding to given number of observed unique barcodes.
#' 
#' Function calculates EUC's for each number of observed barcodes accounting
#' for differential ligation probability of different barcodes. Function
#' k2n_calc() writes file with a vector in which an i-th element is an
#' estimated unique count given observing i unique barcodes.
#' 
#' %% ~~ If necessary, more details than the description above ~~
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
#' @note %% ~~further notes~~
#' @author Lukasz Jan Kielpinski, Nikos Sidiropoulos
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Kielpinski, L.J., and Vinther, J. (2014). Massive
#' parallel-sequencing-based hydroxyl radical probing of RNA accessibility.
#' Nucleic Acids Res.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##-- or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function (merged_file, unique_barcode_file, output_file)
#' {
#'     merged <- read.table(merged_file)
#'     barcode_length <- max(nchar(as.character(merged[, 4])))
#'     read_counts <- count(merged[1:3])
#'     max_observed <- max(read.table(unique_barcode_file)[, 4])
#'     barcodes_nt <- merged[do.call(paste, as.list(merged[, 1:3])) %in%
#'         do.call(paste, as.list(read_counts[(read_counts[, 4]) <=
#'             summary(read_counts[, 4])[5], 1:3])), 4]
#'     nt_counts <- matrix(nrow = 4, ncol = barcode_length)
#'     for (h in 1:ncol(nt_counts)) {
#'         j <- 1
#'         for (nt_local in c("A", "C", "G", "T")) {
#'             nt_counts[j, h] <- sum(substr(as.character(barcodes_nt),
#'                 h, h) == nt_local)
#'             j <- j + 1
#'         }
#'     }
#'     nt_freqs <- nt_counts/colSums(nt_counts)
#'     nt_values <- list()
#'     for (i in 1:ncol(nt_freqs)) {
#'         nt_values[[i]] <- nt_freqs[, i]
#'     }
#'     all_posible_comb <- expand.grid(nt_values)
#'     probs <- apply(all_posible_comb, 1, prod)
#'     results <- c()
#'     i <- 1
#'     results[i] <- sum(1 - (1 - probs)^i)
#'     while (results[i] <= max_observed) {
#'         i <- i + 1
#'         results[i] <- sum(1 - (1 - probs)^i)
#'     }
#'     k2n <- c()
#'     for (i in 1:floor(max(results))) {
#'         difference <- abs(results - i)
#'         k2n[i] <- which(difference == min(difference))
#'     }
#'     if (!missing(output_file)) {
#'         write(k2n, file = output_file)
#'     }
#'     else {
#'         return(k2n)
#'     }
#'   }
#' 
#' @import plyr
#' @export k2n_calc
k2n_calc <- function(merged_file, unique_barcode_file, output_file){


    # Read in inputs:
    merged <- read.table( merged_file )
    barcode_length <- max(nchar(as.character(merged[,4])))
    read_counts <- count(merged[1:3])
    max_observed <- max(read.table(unique_barcode_file)[,4])

    # Remove top-quartile of reads:
    barcodes_nt <- merged[ do.call( paste, as.list( merged[ ,1:3 ] ) ) %in% do.call( paste, as.list( read_counts[ ( read_counts[ , 4 ] ) <= summary( read_counts[ , 4 ] )[ 5 ], 1:3 ] ) ) , 4 ]

    # make the matrix with the nucleotide freqs per position:
    nt_counts <- matrix( nrow = 4, ncol = barcode_length )
    for( h in 1:ncol( nt_counts ) ){
        j <- 1
        for( nt_local in c( "A","C","G","T" ) ) {
            nt_counts[ j, h ] <- sum( substr( as.character( barcodes_nt ), h, h) == nt_local )
        j <- j+1
        }
    }

    # Calculate frequencies
    nt_freqs <- nt_counts / colSums( nt_counts )

    nt_values <- list()
    for( i in 1:ncol( nt_freqs ) ) {
        nt_values[[ i ]] <- nt_freqs[ , i ]
    }

    all_posible_comb <- expand.grid( nt_values )

    probs <- apply( all_posible_comb, 1, prod )

    ###Create Mf_to_Sf:
    results <- c()

    i <- 1
    results[ i ] <- sum( 1 - ( 1 - probs )**i )

    while( results[ i ] <= max_observed ) {
        i <- i + 1
        results[ i ] <- sum( 1 - ( 1 - probs )**i )   #Mf to Sf
    }

    #assign molecules number to unique barcode number:
    k2n <- c()
    for( i in 1:floor( max( results ) ) ) {
        abs( results - i ) -> difference
        k2n[ i ] <- which( difference == min( difference ) ) #if you want to know how many molecules underlie n unique barcodes ask k2n[n]
    }

    if(!missing(output_file)){write(k2n, file = output_file )}else{return(k2n)}
}
