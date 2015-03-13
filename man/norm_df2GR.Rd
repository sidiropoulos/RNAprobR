\name{norm_df2GR}
\alias{norm_df2GR}
\title{Function to convert norm_df data frame (product of GR2norm_df()) to GRanges.}
\usage{
norm_df2GR(norm_df)
}
\arguments{
\item{norm_df}{norm_df data frame needs to have columns: RNAid (equivalent
to seqnames in GRanges) and Pos (equivalent to start in GRanges) and
metadata}
}
\value{
GRanges compatible with objects created by normalizing functions
(dtcr(), slograt(), swinsor(), compdata())
}
\description{
Function to convert norm_df data frame (product of GR2norm_df()) to GRanges.
}
\examples{
dummy_norm_df <- data.frame(RNAid="dummyRNA", Pos=1:100,
                            my_data1=runif(1:100))
norm_df2GR(dummy_norm_df)
}
\author{
Lukasz Jan Kielpinski
}
\seealso{
\code{\link{dtcr}}, \code{\link{slograt}}, \code{\link{swinsor}},
\code{\link{compdata}}, \code{\link{GR2norm_df}}, \code{\link{norm2bedgraph}}
}

