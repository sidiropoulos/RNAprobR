#' Function creating norm_GR GRanges from Comp_GR GRanges which contains
#' termination count (TC), termination-coverage ratio (TCR), coverage (Cover)
#' and priming count (PC) metadata
#' 
#' Main use is to add metadata present in GRanges made by comp() function to
#' GRanges made by normalizing functions (dtcr(), slograt(), swinsor(),
#' compdata())
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param Comp_GR GRanges object made by comp() function.
#' @param nt_offset how many nucleotides before modification the reverse
#' transcription terminates (default: 1)
#' @param add_to normalized data frame with already performed normalization of another kind. Results will be merged
#' @return 
#' \item{comp1 }{Description of 'comp1'}
#' \item{comp2 }{Description of comp2'}
#' @note %% ~~further notes~~
#' @author Lukasz Jan Kielpinski
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' dummy_euc_GR_treated <- GRanges(seqnames="DummyRNA", IRanges(start=round(runif(100)*100), width=round(runif(100)*100+1)), strand="+", EUC=round(runif(100)*100))
#' dummy_comp_GR_treated <- comp(dummy_euc_GR_treated)
#' dummy_swinsor <- swinsor(dummy_comp_GR_treated)
#' dummy_swinsor <- compdata(Comp_GR=dummy_comp_GR_treated, add_to=dummy_swinsor)
#' dummy_swinsor
#' 
#' ## The function is currently defined as
#' function (Comp_GR, nt_offset = 1, add_to) 
#' {
#'     if (nt_offset < 0) {
#'         print("error: nt_offset must be >= 0")
#'         stop()
#'     }
#'     offset_oneRNA <- function(oneRNA_comp_df) {
#'         if (prod(diff(oneRNA_comp_df$Pos) == rep(1, nrow(oneRNA_comp_df) - 
#'             1)) == 1) {
#'             oneRNA_comp_df[, columns_to_extract] <- oneRNA_comp_df[(1 + 
#'                 nt_offset):(nrow(oneRNA_comp_df) + nt_offset), 
#'                 columns_to_extract]
#'             return(oneRNA_comp_df)
#'         }
#'         else {
#'             stop(paste("Check if data was properly sorted by comp() function. Problem with", 
#'                 oneRNA_comp_df$RNAid[1]))
#'         }
#'     }
#'     Comp_df <- GR2norm_df(Comp_GR)
#'     columns_to_extract <- names(Comp_df) %in% c("TC", "TCR", 
#'         "Cover", "PC")
#'     Comp_df[is.na(Comp_df)] <- 0
#'     Comp_df <- Comp_df[order(Comp_df$RNAid, Comp_df$Pos), ]
#'     Comp_by_RNA <- split(Comp_df, f = Comp_df$RNAid, drop = T)
#'     normalized <- do.call(rbind, lapply(Comp_by_RNA, FUN = offset_oneRNA))
#'     normalized[is.na(normalized)] <- 0
#'     if (!missing(add_to)) {
#'         add_to_df <- GR2norm_df(add_to)
#'         normalized <- merge(add_to_df, normalized, by = c("RNAid", 
#'             "Pos", "nt"), suffixes = c(".old", ".new"))
#'     }
#'     normalized <- normalized[order(normalized$RNAid, normalized$Pos), 
#'         ]
#'     normalized_GR <- norm_df2GR(normalized)
#'     return(normalized_GR)
#'   }
#' 
#' @import GenomicRanges
#' @export compdata
compdata <- function(Comp_GR, nt_offset=1, add_to){

###Check conditions:
	if(nt_offset < 0){
		print("error: nt_offset must be >= 0")
		stop()
	}
###Define functions:
offset_oneRNA <- function(oneRNA_comp_df){
	if(prod(diff(oneRNA_comp_df$Pos)==rep(1, nrow(oneRNA_comp_df)-1))==1){ #Checks if data is properly sorted.
	oneRNA_comp_df[,columns_to_extract] <- oneRNA_comp_df[(1+nt_offset):(nrow(oneRNA_comp_df)+nt_offset), columns_to_extract]
	return(oneRNA_comp_df)}else{
	stop(paste("Check if data was properly sorted by comp() function. Problem with",oneRNA_comp_df$RNAid[1]))
	}
}
###Function body:
Comp_df <- GR2norm_df(Comp_GR)

columns_to_extract <- names(Comp_df) %in% c("TC","TCR","Cover","PC")

Comp_df[is.na(Comp_df)] <- 0

Comp_df <- Comp_df[order(Comp_df$RNAid, Comp_df$Pos),]

Comp_by_RNA <- split(Comp_df, f=Comp_df$RNAid, drop=T)

normalized <- do.call(rbind, lapply(Comp_by_RNA, FUN=offset_oneRNA))

normalized[is.na(normalized)] <- 0

#If add_to specified, merge with existing normalized data frame:
if(!missing(add_to)){
add_to_df <- GR2norm_df(add_to)
normalized <- merge(add_to_df, normalized, by=c("RNAid", "Pos", "nt"), suffixes=c(".old",".new"))
}
###

normalized <- normalized[order(normalized$RNAid, normalized$Pos),]
normalized_GR <- norm_df2GR(normalized)

return(normalized_GR)
}
