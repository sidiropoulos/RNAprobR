\name{RNA2df}
\alias{RNA2df}
\title{
Function to make data frame with normalized results for chosen transcripts out of GRanges
}
\description{
Specific information extraction from GRanges object created by normalization functions (dtcr(), swinsor(), slograt() or compdata())
Similar to GR2norm_df() but prints ony desired subset of GRanges.
}
\usage{
RNA2df(norm_GR, RNAid, norm_methods)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{norm_GR}{
GRanges object made by other normalization function (dtcr(), slograt(), swinsor(), compdata()) from which data is to be extracted
}
  \item{RNAid}{
Transcript identifiers of transcripts that are to be extracted
}
  \item{norm_methods}{
Names of the columns to be extracted.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Data frame object with columns: RNAid, Pos and desired metadata columns (e.g. nt, dtcr)
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

dummy_euc_GR_treated <- GRanges(seqnames="DummyRNA", IRanges(start=round(runif(100)*100), width=round(runif(100)*100+1)), strand="+", EUC=round(runif(100)*100))
dummy_comp_GR_treated <- comp(dummy_euc_GR_treated)
dummy_swinsor <- swinsor(dummy_comp_GR_treated)
RNA2df(dummy_swinsor, RNAid="DummyRNA", norm_methods="swinsor")

## The function is currently defined as
function (norm_GR, RNAid, norm_methods) 
{
    if (missing(norm_methods)) {
        norm_methods <- names(mcols(norm_GR))
    }
    oneRNA_GR <- norm_GR[seqnames(norm_GR) \%in\% RNAid, norm_methods]
    GR2norm_df(oneRNA_GR)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
