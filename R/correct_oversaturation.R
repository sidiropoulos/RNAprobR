#' Correcting EUC of oversaturated fragments.
#'
#' If for a given fragment the number of observed unique barcodes is equal to
#' the total barcode complexity (all combinations of barcodes are associated
#' with a given fragment), then the readsamples function assignes infinite EUC.
#' This can be corrected by the function correct_oversaturation(). By comparing
#' observed read counts with EUCs for other fragments it calculates the
#' correction factor.
#' Then, for the oversaturated fragments it multiplies the observed read counts
#' by the correction factor to estimate EUC. The assumption behind this
#' correction is that fragments have similar rate of PCR duplicates production.
#'
#' @param euc_GR GRanges produced by readsamples() function
#' @param read_counts_file path to a file with observed read counts.
#' @return euc_GR GRanges analogous to the readsamples() function output,
#' but with finite EUCs where infinity was present.
#' @examples
#' write(c("DummyRNA\t1\t2\t1000", "DummyRNA\t3\t4\t1024"),
#'      file="dummy_unique_barcode")
#' write(c("DummyRNA\t1\t2\t5000", "DummyRNA\t3\t4\t10000"),
#'       file="dummy_read_counts")
#' my_EUCs <- readsamples(samples = "dummy_unique_barcode", euc = "Fu", m=1024)
#' correct_oversaturation(euc_GR = my_EUCs,
#'                        read_counts_file = "dummy_read_counts")
#'
#' @importFrom stats lm
#' @export

correct_oversaturation <- function(euc_GR, read_counts_file){

    #Read in read counts file:
    read_counts_euc <- readsamples(read_counts_file, euc="counts")
    read_counts <- read_counts_euc[match(euc_GR, read_counts_euc)]$EUC

    #Calculate correction factor using linear model:
    scores_matrix <- matrix(c(euc_GR$EUC, read_counts), ncol=2)
    #(do not use points with Inf or "sporadic" fragments)
    scores_matrix[scores_matrix==Inf] <- NA
    scores_matrix[scores_matrix==1] <- NA
    cor_fact <- lm(scores_matrix[,1] ~ scores_matrix[,2] + 0)$coefficients

    #Which measurements need to be corrected:
    to_correct <- which(euc_GR$EUC == Inf)

    #Correct and return corrected:
    euc_GR$EUC[to_correct] <- round(read_counts[to_correct] * cor_fact)
    euc_GR
}
