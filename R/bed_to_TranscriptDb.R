###Function to transform BED format file to Bioconductor TranscriptDb object
# Description of options:
# input_bed_path - path to BED file (12 column BED needed)
# Value: TranscriptDb object

BED2txDb <- function(input_bed_path){
require(GenomicFeatures)

#Define functions:
rank_fun <- function(exon_count, strandness){
  ranks <- 1:exon_count
  if(strandness=="-"){ranks <- rev(ranks)}
  ranks
  }
#End of define functions

##Read in BED file and add unique IDs:
input_bed <- read.table(input_bed_path)

if(is.element("tx_id", names(input_bed))){                   #This condition should never be fulfilled, just in case if someone tampers with BED file.
  print("ERROR: BED file not allowed to contain column named tx_id")
  stop()} 

input_bed$tx_id <- (1:nrow(input_bed)) #Add internal transcript ID (required for making TranscriptDb)

##create first data frame required for makeTranscriptDb function (transcripts):
transcripts <- data.frame(tx_id=input_bed$tx_id, tx_name=input_bed[,4], tx_chrom=input_bed[,1], tx_strand=input_bed[,6], tx_start=input_bed[,2]+1, tx_end=input_bed[,3])

##create second data frame required for makeTranscriptDb function (splicing):
input_bed[,11] <- as.character(input_bed[,11]) #column 11 of BED file - exon sizes;  change factor to character
input_bed[,12] <- as.character(input_bed[,12]) #column 12 of BED file - exon starts; change factor to character

exon_size_tx_list <- lapply(strsplit(input_bed[,11],","), as.integer)   #Split strings into list of integers
exon_start_tx_list <- lapply(strsplit(input_bed[,12],","), as.integer) #Split strings into list of integers

exon_start_genome_list <- mapply(FUN=function(x,y){x+y+1}, exon_start_tx_list, input_bed[,2], SIMPLIFY = F)       #Calculate exon starts in genome coordinates, and change to 1-counted notation.
exon_end_genome_list <- mapply(FUN=function(x,y){x+y-1}, exon_start_genome_list, exon_size_tx_list, SIMPLIFY = F) #Calculate exon ends in genome coordinates, and change to 1-counted notation.

exon_rank_list <- mapply(FUN=rank_fun, input_bed[,10], input_bed[,6], SIMPLIFY = F) #Calcualte exon rank: ordering within transcript

tx_id_list <- mapply(FUN=rep, input_bed$tx_id, input_bed[,10], SIMPLIFY = F) #transcript ID

splicings <- data.frame(tx_id=unlist(tx_id_list), exon_rank=unlist(exon_rank_list), exon_start=unlist(exon_start_genome_list), exon_end=unlist(exon_end_genome_list))

##Make TranscriptDb object:
txDb_from_BED <- makeTranscriptDb(transcripts=transcripts, splicings=splicings)
txDb_from_BED
}

