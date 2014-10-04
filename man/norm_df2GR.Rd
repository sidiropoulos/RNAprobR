\name{norm_df2GR}
\alias{norm_df2GR}
\title{
Function to convert norm_df data frame (product of GR2norm_df()) to GRanges
}
\description{
Function to convert norm_df data frame (product of GR2norm_df()) to GRanges. 
}
\usage{
norm_df2GR(norm_df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{norm_df}{
norm_df data frame needs to have columns: RNAid (equivalent to seqnames in GRanges) and Pos (equivalent to start in GRanges) and metadata
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
GRanges compatible with objects created by normalizing functions (dtcr(), slograt(), swinsor(), compdata())
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
dummy_norm_df <- data.frame(RNAid="dummyRNA", Pos=1:100, my_data1=runif(1:100))
norm_df2GR(dummy_norm_df)

## The function is currently defined as
function (norm_df) 
{
    GRanges(seqnames = norm_df$RNAid, IRanges(start = norm_df$Pos, 
        width = 1), strand = "+", norm_df[names(norm_df) != 
        "RNAid" & names(norm_df) != "Pos"])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
