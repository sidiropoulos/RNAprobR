\name{BED2txDb}
\alias{BED2txDb}
\title{
Bedgraph to TranscriptDb object
}
\description{
Function to transform BED format file to Bioconductor TranscriptDb object
}
\usage{
BED2txDb(input_bed_path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_bed_path}{
Path to 12 column BED file
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lukasz Jan Kielpinski
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

write("chr1\t134212702\t134229870\tENSMUST00000072177\t0\t+\t134212806\t134228958\t0\t8\t347,121,24,152,66,120,133,1973,\t0,8827,10080,11571,12005,13832,14433,15195,", file="dummy.bed")
BED2txDb("dummy.bed")

## The function is currently defined as
function (input_bed_path) 
{
    rank_fun <- function(exon_count, strandness) {
        ranks <- 1:exon_count
        if (strandness == "-") {
            ranks <- rev(ranks)
        }
        ranks
    }
    input_bed <- read.table(input_bed_path)
    if (is.element("tx_id", names(input_bed))) {
        print("ERROR: BED file not allowed to contain column named tx_id")
        stop()
    }
    input_bed$tx_id <- (1:nrow(input_bed))
    transcripts <- data.frame(tx_id = input_bed$tx_id, tx_name = input_bed[, 
        4], tx_chrom = input_bed[, 1], tx_strand = input_bed[, 
        6], tx_start = input_bed[, 2] + 1, tx_end = input_bed[, 
        3])
    input_bed[, 11] <- as.character(input_bed[, 11])
    input_bed[, 12] <- as.character(input_bed[, 12])
    exon_size_tx_list <- lapply(strsplit(input_bed[, 11], ","), 
        as.integer)
    exon_start_tx_list <- lapply(strsplit(input_bed[, 12], ","), 
        as.integer)
    exon_start_genome_list <- mapply(FUN = function(x, y) {
        x + y + 1
    }, exon_start_tx_list, input_bed[, 2], SIMPLIFY = F)
    exon_end_genome_list <- mapply(FUN = function(x, y) {
        x + y - 1
    }, exon_start_genome_list, exon_size_tx_list, SIMPLIFY = F)
    exon_rank_list <- mapply(FUN = rank_fun, input_bed[, 10], 
        input_bed[, 6], SIMPLIFY = F)
    tx_id_list <- mapply(FUN = rep, input_bed$tx_id, input_bed[, 
        10], SIMPLIFY = F)
    splicings <- data.frame(tx_id = unlist(tx_id_list), exon_rank = unlist(exon_rank_list), 
        exon_start = unlist(exon_start_genome_list), exon_end = unlist(exon_end_genome_list))
    txDb_from_BED <- makeTranscriptDb(transcripts = transcripts, 
        splicings = splicings)
    txDb_from_BED
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
