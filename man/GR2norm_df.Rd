\name{GR2norm_df}
\alias{GR2norm_df}
\title{Export normalized GRanges object to data frame}
\usage{
GR2norm_df(norm_GR, RNAid = "all", norm_methods = "all")
}
\arguments{
\item{norm_GR}{GRanges object made by other normalization function (dtcr(),
slograt(), swinsor(), compdata()) from which data is to be extracted}

\item{RNAid}{Transcript identifiers of transcripts that are to be extracted}

\item{norm_methods}{Names of the columns to be extracted.}
}
\value{
Data frame object with columns: RNAid, Pos and desired metadata
columns (e.g. nt, dtcr)
}
\description{
Function to make data frame out of GRanges output of normalizing functions
(dtcr(), slograt(), swinsor(), compdata()) for all or a set of chosen
transcripts in the file.
}
\examples{
dummy_euc_GR_treated <- GRanges(seqnames="DummyRNA",
                                IRanges(start=round(runif(100)*100),
                                width=round(runif(100)*100+1)), strand="+",
                                EUC=round(runif(100)*100))
dummy_comp_GR_treated <- comp(dummy_euc_GR_treated)
dummy_swinsor <- swinsor(dummy_comp_GR_treated)
GR2norm_df(dummy_swinsor)
}
\author{
Lukasz Jan Kielpinski, Nikos Sidiropoulos
}
\seealso{
\code{\link{norm_df2GR}}, \code{\link{dtcr}},
\code{\link{swinsor}}, \code{\link{slograt}}, \code{\link{compdata}}
}

