\name{norm2bedgraph}
\alias{norm2bedgraph}
\title{Exporting data in norm_df data frame (product of dtcr, slograt and
swinsor) to bedgraph format compatible with UCSC Genome Browser}
\usage{
norm2bedgraph(norm_GR, txDb, bed_file, norm_method, genome_build,
  bedgraph_out_file = "out_file", track_name = "Track_name",
  track_description = "Track_description")
}
\arguments{
\item{norm_GR}{norm_GR GRanges with data to be exported, required}

\item{txDb}{TranscriptDb object with transcript definitions. Names must
match those in norm_df}

\item{bed_file}{character containing file path to BED file with transcript
definitions. Supply txDb XOR bedfile}

\item{norm_method}{character specifying which normalized column
should be processed into bedgraph. If not provided, the first column
matching dtcr, slograt or swinsor is transformed.}

\item{genome_build}{character specifying which UCSC genome build should data
be displayed in, e.g. "mm9"}

\item{bedgraph_out_file}{character specifying prefix of output file.
Generated file name is: prefix.bedgraph; if file with such a name already
exists new tracks will be appended.}

\item{track_name}{character specifying track name}

\item{track_description}{character specifying track description}
}
\value{
Function writes bedgraph file.
}
\description{
Function converts annotation from transcript to genomic coordinates and
creates two-track bedgraph file (one track for each strand)
}
\examples{
dummy_euc_GR_control <- GRanges(seqnames="DummyRNA",
                                 IRanges(start=round(runif(100)*100),
                                 width=round(runif(100)*100+1)), strand="+",
                                 EUC=round(runif(100)*100))
dummy_euc_GR_treated <- GRanges(seqnames="DummyRNA",
                                IRanges(start=round(runif(100)*100),
                                width=round(runif(100)*100+1)), strand="+",
                                EUC=round(runif(100)*100))
dummy_comp_GR_control <- comp(dummy_euc_GR_control)
dummy_comp_GR_treated <- comp(dummy_euc_GR_treated)
dummy_norm <- dtcr(control_GR=dummy_comp_GR_control,
                   treated_GR=dummy_comp_GR_treated)
write(paste(c("chr1", 134212702, 134229870, "DummyRNA", 0, "+", 134212806,
            134228958, 0, 8, "347,121,24,152,66,120,133,1973,",
            "0,8827,10080,11571,12005,13832,14433,15195,"), collapse = "\\t"),
      file="dummy.bed")
norm2bedgraph(norm_GR = dummy_norm, bed_file = "dummy.bed")

}
\seealso{
\code{\link{bedgraph2norm}}, \code{\link{norm_df2GR}},
\code{\link{dtcr}}, \code{\link{slograt}}, \code{\link{swinsor}},
\code{\link{compdata}}
}
\author{
Lukasz Jan Kielpinski, Nikos Sidiropoulos
}
