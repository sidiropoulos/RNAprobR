\name{plotRNA}
\alias{plotRNA}
\title{
Function plotting normalized values over transcript positions
}
\description{
Function plotting normalized values over transcript positions.
}
\usage{
plotRNA(norm_GR, RNAid, norm_method, stat_method, stat_cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{norm_GR}{
norm_GR GRanges with data to be exported, required
}
  \item{RNAid}{
Transcript identifier - for which transcript plot should be generated.
}
  \item{norm_method}{
Names of a column to be used for plotting (normalization method).
}
  \item{stat_method}{
Names of a column to be used for adding significance asterisks. If stat_method not provided, try to guess it, if no guess - empty vector.
}
  \item{stat_cutoff}{
below what value of statistics (from stat_method) report significance. If not provided - minimal value from stat_method used
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plotting function.
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
plotRNA(dummy_swinsor, RNAid="DummyRNA")


## The function is currently defined as
function (norm_GR, RNAid, norm_method, stat_method, stat_cutoff) 
{
    if (missing(norm_method)) {
        norm_method <- names(mcols(norm_GR))[which(names(mcols(norm_GR)) \%in\% 
            c("dtcr", "slograt", "swinsor"))[1]]
    }
    if (length(norm_method) != 1) {
        stop("Specify 1 method")
    }
    if (missing(stat_method)) {
        stat_method <- switch(norm_method, dtcr = "dtcr.p", slograt = "slograt.p", 
            swinsor = "swinsor.sd", character(0))
    }
    oneRNA_GR <- norm_GR[seqnames(norm_GR) == RNAid, c(norm_method, 
        stat_method)]
    plot(mcols(oneRNA_GR)[, 1] ~ start(oneRNA_GR), xlab = "Position", 
        ylab = norm_method, main = RNAid, type = "l")
    if (length(stat_method) == 1) {
        if (missing(stat_cutoff)) {
            stat_cutoff <- min(mcols(oneRNA_GR)[, 2], na.rm = T)
        }
        stat_values <- rep(NA, length(oneRNA_GR))
        stat_values[mcols(oneRNA_GR)[, 2] <= stat_cutoff] <- max(mcols(oneRNA_GR)[, 
            1], na.rm = T) * 1.01
        points(stat_values ~ start(oneRNA_GR), col = "red", pch = "*")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
