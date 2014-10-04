###Function exporting data in norm_df data frame (product of dtcr, slograt and swinsor) to bedgraph format compatible with UCSC Genome Browser 
# Description of options:
# norm_GR - norm_GR GRanges, required
# txDb - TranscriptDb object with transcript definitions. Names must match those in norm_df
# bed_file - character containing file path to BED file with transcript definitions. Supply txDb XOR bedfile
# norm_method - character specifying which normalized values that norm_df contains should be processed into bedgraph. If not provided first column matching dtcr, slograt or swinsor is transformed.
# genome_build - character specifying which UCSC genome build should data be displayed in, e.g. "mm9"
# bedgraph_out_file - character specifying prefix of output file. Generated file names are: prefix.bedgraph
# track_name - character specifying track name
# track_description  - character specifying track description
##Function writes bedgraph file containing 2 tracks (one for each strand).
#' @export
#' 
norm2bedgraph <- function(norm_GR, txDb, bed_file, norm_method, genome_build, bedgraph_out_file="out_file",track_name="Track_name", track_description="Track_description"){
require(GenomicFeatures)
###Check conditions:
if(missing(genome_build)){
	print("ERROR: Specify genome_build")
	stop()
}
if(missing(norm_method)){
	norm_method <- names(norm_df)[names(norm_df) %in% c("dtcr","slograt","swinsor")][1]
	print(paste("Warning: normalization method to convert not specified.",norm_method,"chosen."))
}
if(missing(txDb) & missing(bed_file)){
	print("Error: specify gene annotation")
	stop()
}
#Define functions:
#Function which removes transcripts for which there is no genomic info and positions within transcripts that are outside annotation (needs norm_df data frame and GRangesList)
remove_unannotated <- function(norm_df, exons_GR){ 
	require(GenomicRanges)
	max_lengths <- sum(width(exons_GR))
	tx_lengths <- data.frame(RNAid=names(max_lengths), max_lengths)
	norm_max_length <- merge(norm_df, tx_lengths, by="RNAid", sort=F)
	discarded_transcripts <- as.character(unique(norm_df[,1])[!(unique(norm_df[,1]) %in% unique(norm_max_length[,1]))])
	if(length(discarded_transcripts) > 0){print(paste("Warning: transcript",discarded_transcripts,"discarded. Genomic location not provided."))}
	norm_out <- subset(norm_max_length, Pos <= max_lengths & Pos > 0)
	norm_out <- norm_out[-which(names(norm_out)=="max_lengths")]
	norm_out <- norm_out[order(norm_out$RNAid, norm_out$Pos),]
	norm_out <- subset(norm_out, Pos > 0)
	norm_out
}
#Function to compress bedgraph by reducing equi-valued runs.
compress_bedgraph <- function(bedgraph_dataframe){
##Define internal function:
	compress_one_chromosome <- function(one_chrom){
		one_chrom <- one_chrom[order(one_chrom[,2]),] #order by position within chromosome
		diff_matrix <- matrix(c(diff(one_chrom[,2]), diff(one_chrom[,4])), ncol=2) #Calculate differences between consecutive positions. delta(position) and delta(value)
		cons_identical <- diff_matrix[,1]==1 & diff_matrix[,2]==0 #Is repeated? or which consecutive (delta(position)==1) positions have the same value (delta(value)==0)
		seq_start_and_end <- diff(cons_identical) #Is the start (+1) or the end (-1) of the repeated sequence
		
		if(sum(seq_start_and_end)==  1){repeat_info <- data.frame(one_chrom, c(F, cons_identical), my_delta=c(0, diff(cons_identical), -1))} #If sum==1, then there is no repeated sequence end - must be at the end
		if(sum(seq_start_and_end)==  0){
			if(min(which(seq_start_and_end==1)) < min(which(seq_start_and_end==-1))){repeat_info <- data.frame(one_chrom, c(F, cons_identical), my_delta=c(0, diff(cons_identical),  0))}
			if(min(which(seq_start_and_end==1)) > min(which(seq_start_and_end==-1))){repeat_info <- data.frame(one_chrom, c(F, cons_identical), my_delta=c(1, diff(cons_identical),  -1))}
			} 
		if(sum(seq_start_and_end)== -1){repeat_info <- data.frame(one_chrom, c(F, cons_identical), my_delta=c(1, diff(cons_identical),  0))} #If sum==-1, then there is no repeated sequence start - must be at the beginning

		repeat_info <- repeat_info[repeat_info[,5]==F | repeat_info[,6]==-1,] #Remove repeated, except first and last
		repeat_info[repeat_info[,6]==1,3] <- repeat_info[repeat_info[,6]== -1,3] #Move end position (column 3) of last repeated to end position (column 3) of first repeated
		repeat_info <- repeat_info[repeat_info[,6]!= -1,] #Remove last repeated
		repeat_info[,1:4] #Return first 4 columns
	}
##End of defining functions
#Run function:
df_list <- split(bedgraph_dataframe, f=bedgraph_dataframe[,1], drop=T)
return(do.call(rbind, lapply(df_list, FUN=compress_one_chromosome)))
}

#End of defining functions.
#Function body:
norm_df <- GR2norm_df(norm_GR)
normalized_column <- which(names(norm_df)==norm_method) #Which column in the norm_df is data to be displayed

if(missing(txDb)){txDb <- BED2txDb(bed_file)} #if transcriptDb not specified, make it from BED
my_exons <- exonsBy(txDb, by="tx", use.names=T) #make GRangesList from TranscriptDb

###Remove transcripts duplicated in transcriptDb:
dup_txs <- names(my_exons)[duplicated(names(my_exons))]
my_exons <- my_exons[!(names(my_exons) %in% dup_txs)]
#Print warning:
my_dup_discarded <- unique(norm_df$RNAid[norm_df$RNAid %in% dup_txs])
if(length(my_dup_discarded) > 0){print(paste("Warning: transcript",my_dup_discarded,"discarded. Provided more than one genomic location."))}
###

norm_good <- remove_unannotated(norm_df, my_exons) #Remove unannotated transcripts and parts of transcripts guided by GRangesList 
norm_list <- split(norm_good, f=norm_good$RNAid, drop=T) #Split normalized data frame into list of data frames\
Pos_list <- lapply(norm_list, FUN=function(norm_df){as.integer(norm_df$Pos)}) #List with positions within transcripts

exons_of_interest <- my_exons[names(my_exons) %in% names(Pos_list)] #GRangesList with transcript present in Pos_list

##Order GRangesList  'exons_if_interest' according to ordering in Pos_list:
ordering_df_Pos <- data.frame(tx_name=names(Pos_list), i=1:length(Pos_list)) #order in data frame. Pos_list is a reference.
ordering_df_exons <- data.frame(tx_name=names(exons_of_interest), i=1:length(exons_of_interest))
ordering_merged <- merge(ordering_df_Pos, ordering_df_exons, by=1, all=T, suffixes=c(".Pos",".ex"), sort=F)
ordering_merged <- ordering_merged[order(ordering_merged[,2]),]
exons_of_interest <- exons_of_interest[ordering_merged[,3]]
##End of ordering

chromosomes <- as.character(runValue(seqnames(exons_of_interest))) #Vector with chromosome names (mathes to Pos_list)
strands <- as.character(runValue(strand(exons_of_interest))) #Vector with strands (mathes to Pos_list)
Pos_genomic_list <- transcriptLocs2refLocs(tlocs=Pos_list, exonStarts=start(exons_of_interest), exonEnds=end(exons_of_interest), strand=strands) #List with genomic coordinates, matches Pos_list

#Merge chromosome, genomic position and strand (first list of data frames...)
df_list <- list()
	for(i in 1:length(Pos_genomic_list)){
	df_list[[i]] <- data.frame(seqname=chromosomes[i], position_off=Pos_genomic_list[[i]]-1, position=Pos_genomic_list[[i]], value=norm_list[[i]][,normalized_column], strand=strands[i])
}
df_for_bedgraph <- do.call(rbind, df_list) #...then merge data frames.

#Look for duplicates:
dups <- duplicated(df_for_bedgraph[,c(1:3,5)]) | duplicated(df_for_bedgraph[,c(1:3,5)], fromLast=TRUE)
dup_count <- sum(dups)
if(dup_count > 0){
	df_for_bedgraph <- df_for_bedgraph[!dups,]
	print(paste("Warning:",dup_count,"lines were removed from dataset due to duplication. Your transcript annotations overlap."))
}
#Remove NA values:

df_for_bedgraph <- df_for_bedgraph[!is.na(df_for_bedgraph$value),]

#And finally split and export two tracks (one for each strand) to a file:
df_plus <- df_for_bedgraph[df_for_bedgraph[,5]=="+",1:4]
df_minus <- df_for_bedgraph[df_for_bedgraph[,5]=="-",1:4]

#Compress bedgraph:
df_plus <- compress_bedgraph(df_plus)
df_minus <- compress_bedgraph(df_minus)

#End of compression

bedgraph_out_file <- paste(bedgraph_out_file,".bedgraph", sep="")
bedgraph_header <- paste('track type=bedGraph name="',track_name,'(plus)" description="',track_description,'-plus_strand" visibility=full color=0,0,100 altColor=0,0,0 priority=100 autoScale=on alwaysZero=on gridDefault=off maxHeightPixels=128:128:11 graphType=bar yLineMark=0 yLineOnOff=on smoothingWindow=off db=',genome_build, sep="")
write(bedgraph_header, file=bedgraph_out_file)
write.table(df_plus, row.names=F, col.names=F, quote=F, sep="\t", file=bedgraph_out_file, append=T)

bedgraph_header <- paste('track type=bedGraph name="',track_name,'(minus)" description="',track_description,'-minus_strand" visibility=full color=0,0,100 altColor=0,0,0 priority=100 autoScale=on alwaysZero=on gridDefault=off maxHeightPixels=128:128:11 graphType=bar yLineMark=0 yLineOnOff=on smoothingWindow=off db=',genome_build, sep="")
write(bedgraph_header, file=bedgraph_out_file, append=T)
write.table(df_minus, row.names=F, col.names=F, quote=F, sep="\t", file=bedgraph_out_file, append=T)
}
