\name{compdata}
\alias{compdata}
\title{Create or extend norm_GR GRanges using Comp_GR GRanges}
\usage{
compdata(Comp_GR, nt_offset = 1, add_to)
}
\arguments{
\item{Comp_GR}{GRanges object made by comp() function.}

\item{nt_offset}{how many nucleotides before modification the reverse
transcription terminates (default: 1)}

\item{add_to}{normalized data frame with already performed normalization of
another kind. Results will be merged}
}
\value{
\item{norm_GR}{norm_GR GRanges extended by metadata from Comp_GR}
}
\description{
Add metadata present in GRanges made by comp() function
(termination count (TC), termination-coverage ratio (TCR), coverage (Cover)
and priming count (PC)) to GRanges made by normalizing functions (dtcr(),
slograt(), swinsor(), compdata()).
}
\examples{

dummy_euc_GR_treated <- GRanges(seqnames="DummyRNA",
                                IRanges(start=round(runif(100)*100),
                                width=round(runif(100)*100+1)), strand="+",
                                EUC=round(runif(100)*100))
dummy_comp_GR_treated <- comp(dummy_euc_GR_treated)
dummy_swinsor <- swinsor(dummy_comp_GR_treated)
dummy_swinsor <- compdata(Comp_GR=dummy_comp_GR_treated,
                          add_to=dummy_swinsor)
dummy_swinsor

}
\seealso{
\code{\link{comp}}, \code{\link{dtcr}}, \code{\link{slograt}},
\code{\link{swinsor}}, \code{\link{GR2norm_df}}, \code{\link{plotRNA}},
\code{\link{norm2bedgraph}}
}
\author{
Lukasz Jan Kielpinski, Nikos Sidiropoulos
}
