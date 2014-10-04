# Function creating norm_GR from Comp_GR which contains termination count (TC), termination-coverage ratio (TCR), coverage (Cover) and priming count (PC) metadata
# Options:
# Comp_GR - Comp_GR GRanges input (required; created by comp() function)
# nt_offset - how many nucleotides before modification the reverse transcription terminates (default: 0)
# add_to - normalized data frame with already performed normalization of another kind. Results will be merged
#' @export

compdata <- function(Comp_GR, nt_offset=1, add_to){
require(GenomicRanges)
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
