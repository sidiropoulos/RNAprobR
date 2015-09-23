#' Function converts bam file to bedgraph by counting number of reads starting at each position (termination counts)  
#' It creates two-track bedgraph file (one track for each strand).
#' 
#' #' @param bam_path path to a bam file to be converted
#' @param allowed_flags integer vector with SAM flags should be kept, see https://broadinstitute.github.io/picard/explain-flags.html for explanation
#' @param maxMemory maxMemory of scanBam function used internally
#' @param genome_build character specifying which UCSC genome build should data
#' be displayed in, e.g. "mm9"
#' @param bedgraph_out_file character specifying prefix of output file.
#' Generated file name is: prefix.bedgraph; if file with such a name already
#' exists new tracks will be appended.
#' @param track_name character specifying track name
#' @param track_description character specifying track description
#' 
#' @author Lukasz Jan Kielpinski
#' @import Rsamtools
#' @export bam2bedgraph

bam2bedgraph <- function(bam_path, allowed_flags=0:4095, maxMemory=8000, genome_build,
                         bedgraph_out_file="out_file",track_name="Track_name",
                         track_description="Track_description"){

  
  #Read in bam file:
  input_bam_ref <- BamFile(bam_path)
  scanned_bam <- scanBam(input_bam_ref, maxMemory=maxMemory, param=ScanBamParam(what=c("flag","rname","strand","pos","cigar")))
  
  #Filter on flags:
  kept_flag <- scanned_bam[[1]]$flag %in% allowed_flags
  
  #Extract termini and compress:
  for(local_strand in c("+","-")){
    kept_strand <- scanned_bam[[1]]$strand == local_strand
    kept <- kept_flag & kept_strand
    
    if(any(kept)){
      if(local_strand == "+"){
        positions <- scanned_bam[[1]]$pos[kept]
      }else{
        positions <- scanned_bam[[1]]$pos[kept] + GenomicAlignments::cigarWidthAlongReferenceSpace(scanned_bam[[1]]$cigar[kept]) - 1
      }
      
      chr_pos <- data.frame(seqname=scanned_bam[[1]]$rname[kept], position=positions)
      chr_pos$position_off <- chr_pos$position - 1 
      counts_per_pos <- aggregate(list(value=rep(1, nrow(chr_pos))), by=chr_pos, length)
    
      if(local_strand == "+"){
        df_plus <- .compress_bedgraph(counts_per_pos)
      }else{
        df_minus <- .compress_bedgraph(counts_per_pos)
      }
    }
  }
  
  #####
  
  #Save bedgraph:
  .save_bedgraph(df_plus = df_plus, df_minus = df_minus, genome_build = genome_build, bedgraph_out_file = bedgraph_out_file, track_name = track_name,  track_description = track_description)
  ####
}






