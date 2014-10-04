# Function comp() takes as input EUC_GR produced by readsamples() and produces Comp_GR GRanges with: 
# 1) RNAid, 
# 2) Pos (Position within RNA), 
# and metadata:
# 3) TCR (termination coverage ratio),
# 4) TC (termination count),
# 5) Cover (coverage) and 
# 6) PC (priming count) 
# for each position within each RNA. 
# Description of options: 
# euc_GR - GRanges generated  by readsamples() function
# cutoff - specifies cutoff length, only inserts of this length or longer will be used for processing (default: 0)
# fasta_file - path to fasta file to which reads were mapped. Used to report nucleotide at each position (not required)
#' @export


comp <- function(euc_GR, cutoff=1, fasta_file){
require(GenomicRanges)
require(Rsamtools)

#Check conditions:
	#Check if cutoff is >= 1. If it is not this leads to erroneous coverage calculation:
	if(cutoff < 1){
		print("Cutoff must be >= 1")
		stop()}
	
	#Check if fasta_file is specified, if the file exists - read it in, if it doesn't - print info:
	if(missing(fasta_file)){
		print("Fasta file not specified.")
		fasta_exists <- F
		}else{
			if(file.exists(fasta_file)){
				indexFa(fasta_file)
				fasta_ref <- FaFile(fasta_file)
				fasta_exists <- T
			}else{
				print("Warning: Fasta file not found.")
				fasta_exists <- F
				}
		}

##############
###Define functions:
	process_oneRNA_euc <- function(oneRNA_euc){
		RNAid <- as.character(seqnames(oneRNA_euc[1])) #Name of analysed RNA
		
		RNA_order <- which(names(TC_all)==RNAid) #which coverage vector in the coverage list corresponds to our gene. It is calculated once and assumed to be the same for all three coverage vectors: TC (termination counts), PC (priming counts) and coverage.
		
		if(length(RNA_order) > 1){stop(paste("ERROR: More than one gene with the same ID provided. Check", RNAid))}
		
		#Below I am extracting coverage vectors from coverage list for TC, PC and coverage:
		gene_TC <- as.numeric(TC_all[[RNA_order]]) #coverage to vector
		TC_oneRNA <- data.frame(Pos=1:length(gene_TC), gene_TC)
		
		gene_PC <- as.numeric(PC_all[[RNA_order]]) #coverage to vector
		PC_oneRNA <- data.frame(Pos=1:length(gene_PC), gene_PC)
		
		gene_coverage <- as.numeric(cover_all[[RNA_order]]) #coverage to vector
		Cover_oneRNA <- data.frame(Pos=1:length(gene_coverage), gene_coverage)
		
		#merging by position:
		merged_oneRNA <- merge(TC_oneRNA, PC_oneRNA, by="Pos", all=T)
		merged_oneRNA <- merge(merged_oneRNA, Cover_oneRNA, by="Pos", all=T)
		merged_oneRNA$TCR <- merged_oneRNA$gene_TC/merged_oneRNA$gene_coverage
		
		#And constructing GRanges object for single RNA:
		GR_oneRNA <- GRanges(seqnames=RNAid, IRanges(start=merged_oneRNA$Pos, width=1), strand="+", TC=merged_oneRNA$gene_TC, TCR=merged_oneRNA$TCR, Cover=merged_oneRNA$gene_coverage, PC=merged_oneRNA$gene_PC)
		GR_oneRNA
	}

###Function body:
#Remove inserts shorter than cutoff (keep removed in removed_GR:
good_length <- (width(euc_GR) >= cutoff)
removed_GR <- euc_GR[!good_length]
euc_GR_good <- euc_GR[good_length]

#Calculate coverage:
euc_forCoverage <- euc_GR_good
end(euc_forCoverage) <- end(euc_forCoverage)-cutoff+1
cover_all <- coverage(euc_forCoverage, weight=euc_forCoverage$EUC)

#Calculate TC using coverage function. If element length is set to 1 at the stop site then it corresponds to termination count.
euc_forTC <- euc_GR_good
end(euc_forTC) <- start(euc_forTC)
TC_all <- coverage(euc_forTC, weight=euc_forTC$EUC)

#Calculate PC, the same way as TC:
euc_forPC <- euc_GR_good
start(euc_forPC) <- end(euc_forPC)
PC_all <- coverage(euc_forPC, weight=euc_forPC$EUC)

###Run the single RNA processing for all the RNAs: 
euc_by_RNA <- split(euc_GR_good, f=seqnames(euc_GR_good), drop=T)
Comp_GR <- unlist(suppressWarnings(endoapply(euc_by_RNA, FUN=process_oneRNA_euc)))

#If fasta_file good, add nucleotide identity:
if(fasta_exists){Comp_GR$nt <- as.character(getSeq(fasta_ref, Comp_GR), use.names=F)}else{Comp_GR$nt <- NA}

#Print info on fraction of removed EUC's:
percent_removed <- sum(removed_GR$EUC)/(sum(removed_GR$EUC) + sum(euc_GR_good$EUC))*100
print(paste(round(percent_removed,2), "% of EUCs removed due to cutoff"))

#Return GRanges:
Comp_GR
}

