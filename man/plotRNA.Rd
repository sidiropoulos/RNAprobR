\name{plotRNA}
\alias{plotRNA}
\title{Plot normalized values over transcript positions}
\usage{
plotRNA(norm_GR, RNAid, norm_method, stat_method, stat_cutoff, main, type, ylab,
  xlab, ...)
}
\arguments{
\item{norm_GR}{norm_GR GRanges with data to be exported, required}

\item{RNAid}{Transcript identifier, for which transcript plot should be
generated.}

\item{norm_method}{Which normalization method should be to be used for
plotting (column name).}

\item{stat_method}{Name of a column to be used for adding significance
asterisks. If stat_method not provided, function tries to match with
"norm_method", if no guess - empty vector.}

\item{stat_cutoff}{below what value of statistics (from stat_method, p-value
or standard deviation) report significance. If not provided - minimal value
from stat_method used. To suppress reporting significant sites provide
negative value}

\item{main}{an overall title for the plot: see \code{\link{title}}.}

\item{type}{what type of plot should be drawn. See \code{\link{plot}} for
possible types.}

\item{ylab}{a title for the y axis: see \code{\link{title}}.}

\item{xlab}{a title for the x axis: see \code{\link{title}}.}

\item{...}{Arguments to be passed to methods, such as
\code{\link{graphical parameters}} (see \code{\link{par}}).}
}
\value{
Plotting function.
}
\description{
Function plotting normalized values over transcript positions.
}
\examples{
dummy_euc_GR_treated <- GRanges(seqnames="DummyRNA",
                                IRanges(start=round(runif(100)*100),
                                width=round(runif(100)*100+1)), strand="+",
                                EUC=round(runif(100)*100))
dummy_comp_GR_treated <- comp(dummy_euc_GR_treated)
dummy_swinsor <- swinsor(dummy_comp_GR_treated)
plotRNA(dummy_swinsor, RNAid="DummyRNA")
}
\author{
Lukasz Jan Kielpinski
}
\seealso{
\code{\link{plot}}, \code{\link{plot.default}}, \code{\link{dtcr}},
\code{\link{slograt}}, \code{\link{swinsor}}, \code{\link{compdata}}
}

