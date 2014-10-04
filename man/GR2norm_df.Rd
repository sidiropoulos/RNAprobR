\name{GR2norm_df}
\alias{GR2norm_df}
\title{
Function to make data frame out of GRanges output of normalizing functions (dtcr(), slograt(), swinsor(), compdata())
}
\description{
Simple convenience function. Similar to RNA2df() but prints all the information in the GRanges.
}
\usage{
GR2norm_df(norm_GR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{norm_GR}{
GRanges object made by other normalization function (dtcr(), slograt(), swinsor(), compdata()) from which data is to be extracted
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
GR2norm_df(dummy_swinsor)


## The function is currently defined as
function (norm_GR) 
{
    data.frame(RNAid = as.character(seqnames(norm_GR)), Pos = as.integer(start(norm_GR)), 
        mcols(norm_GR))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
