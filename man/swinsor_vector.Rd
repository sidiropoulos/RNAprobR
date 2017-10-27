\name{swinsor_vector}
\alias{swinsor_vector}
\title{Smooth Winsor Normalization}
\usage{
swinsor_vector(input_vector, window_size, winsor_level = 0.9,
  only_top = FALSE)
}
\arguments{
\item{input_vector}{Vector with values to be smooth-Winsorized}

\item{window_size}{Size of a sliding window.}

\item{winsor_level}{Winsorization level. Bottom outliers will be set to
(1-winsor_level)/2 quantile and top outliers to (1+winsor_level)/2 quantile.}

\item{only_top}{If TRUE then bottom values are not Winsorized and are set to
0.}
}
\value{
\item{comp1}{Vector with mean Winsorized values for each
input_vector position} \item{comp2}{Vector with standard deviation of
Winsorized values for each input_vector position}
}
\description{
Function performs Winsor normalization (see winsor() function) of each
window of specified window_size, sliding in a given vector by 1 position,
and reports a list of (1) mean Winsorized values for each vector position
(mean of Winsorized value for a given position as calculated within each
overlapping window) and (2) standard deviation of those Winsorized values.
}
\examples{

data_set <- runif(1:100)*100
plot(swinsor_vector(data_set, window_size=71,
                    winsor_level=0.8)[[1]] ~ data_set)

}
\references{
"Analysis of sequencing based RNA structure probing data"
Kielpinski, Sidiropoulos, Vinther. Chapter in "Methods in Enzymology"
(in preparation)
}
\author{
Lukasz Jan Kielpinski
}
