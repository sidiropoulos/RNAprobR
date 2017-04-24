#' Create or extend norm_GR GRanges using Comp_GR GRanges
#'
#' Add metadata present in GRanges made by comp() function
#' (termination count (TC), termination-coverage ratio (TCR), coverage (Cover)
#' and priming count (PC)) to GRanges made by normalizing functions (dtcr(),
#' slograt(), swinsor(), compdata()).
#'
#' @param Comp_GR GRanges object made by comp() function.
#' @param nt_offset how many nucleotides before modification the reverse
#' transcription terminates (default: 1)
#' @param add_to normalized data frame with already performed normalization of
#' another kind. Results will be merged
#' @return \item{norm_GR}{norm_GR GRanges extended by metadata from Comp_GR}
#' @author Lukasz Jan Kielpinski, Nikos Sidiropoulos
#' @seealso \code{\link{comp}}, \code{\link{dtcr}}, \code{\link{slograt}},
#' \code{\link{swinsor}}, \code{\link{GR2norm_df}}, \code{\link{plotRNA}},
#' \code{\link{norm2bedgraph}}
#' @examples
#'
#' dummy_euc_GR_treated <- GRanges(seqnames="DummyRNA",
#'                                 IRanges(start=round(runif(100)*100),
#'                                 width=round(runif(100)*100+1)), strand="+",
#'                                 EUC=round(runif(100)*100))
#' dummy_comp_GR_treated <- comp(dummy_euc_GR_treated)
#' dummy_swinsor <- swinsor(dummy_comp_GR_treated)
#' dummy_swinsor <- compdata(Comp_GR=dummy_comp_GR_treated,
#'                           add_to=dummy_swinsor)
#' dummy_swinsor
#'
#' @import GenomicRanges
#' @export compdata
compdata <- function(Comp_GR, nt_offset=1, add_to){

    ###Check conditions:
    if(nt_offset < 0)
        stop("nt_offset must be >= 0")

    ###Function body:
    Comp_df <- GR2norm_df(Comp_GR)

    columns_to_extract <- names(Comp_df) %in% c("TC","TCR","Cover","PC")

    Comp_df[is.na(Comp_df)] <- 0

    Comp_df <- Comp_df[order(Comp_df$RNAid, Comp_df$Pos),]

    Comp_by_RNA <- split(Comp_df, f=Comp_df$RNAid, drop=TRUE)

    normalized <- do.call(rbind, lapply(Comp_by_RNA, FUN=.offset_oneRNA,
                                        nt_offset = nt_offset,
                                        cols = columns_to_extract))

    normalized[is.na(normalized)] <- 0

    #If add_to specified, merge with existing normalized data frame:
    if(!missing(add_to)){
        add_to_df <- GR2norm_df(add_to)
        normalized <- merge(add_to_df, normalized, by=c("RNAid", "Pos", "nt"),
                            suffixes=c(".old",".new"))
    }

    normalized <- normalized[order(normalized$RNAid, normalized$Pos),]
    norm_df2GR(normalized)
}

###Auxiliary functions:

.offset_oneRNA <- function(oneRNA_comp_df, nt_offset, cols){

    #Check if data is properly sorted
    if(prod(diff(oneRNA_comp_df$Pos)==rep(1, nrow(oneRNA_comp_df)-1))==1){
        oneRNA_comp_df[,cols] <-
            oneRNA_comp_df[(1+nt_offset):(nrow(oneRNA_comp_df)+nt_offset),
                           cols]
        oneRNA_comp_df
    }
    else {
        Message <- "Check if data was properly sorted by comp() function.
                    Problem with"
        stop(paste(strwrap(Message), oneRNA_comp_df$RNAid[1]))
    }
}
